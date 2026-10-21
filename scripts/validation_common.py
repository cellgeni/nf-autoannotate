#!/usr/bin/env python3

"""Shared validation helpers used by validate_reference.py and validate_query.py."""

import re
import runpy

import h5py
import numpy as np

import anndata as ad


MATRIX_CHUNK_SIZE = 1_000_000
INTEGER_TOLERANCE = 1e-6

TOOL_CHOICES = ("celltypist", "scanvi", "panhumanpy")
REFERENCE_TOOLS = frozenset({"celltypist", "scanvi"})


def parse_tools(raw_tools: str | None) -> set[str]:
    tools = {value.strip().lower() for value in (raw_tools or "").split(",") if value.strip()}
    unknown = sorted(tools - set(TOOL_CHOICES))
    if unknown:
        raise ValueError(
            f"Unknown --tools entries: {', '.join(unknown)}. Expected one of: {', '.join(TOOL_CHOICES)}"
        )
    if not tools:
        raise ValueError(f"--tools must name at least one of: {', '.join(TOOL_CHOICES)}")
    return tools


def matrix_value_dataset(handle: h5py.File, label: str) -> h5py.Dataset:
    if "X" not in handle:
        raise ValueError(f"Missing required /X matrix in {label}")

    matrix = handle["X"]
    if isinstance(matrix, h5py.Dataset):
        return matrix
    if isinstance(matrix, h5py.Group) and "data" in matrix and isinstance(matrix["data"], h5py.Dataset):
        return matrix["data"]
    raise ValueError(f"Unsupported /X matrix encoding in {label}")


def iter_dataset_chunks(dataset: h5py.Dataset, chunk_size: int = MATRIX_CHUNK_SIZE):
    if dataset.size == 0:
        return

    if len(dataset.shape) == 1:
        for start in range(0, dataset.shape[0], chunk_size):
            yield dataset[start : start + chunk_size]
        return

    if len(dataset.shape) == 2:
        rows_per_chunk = max(1, chunk_size // max(1, dataset.shape[1]))
        for start in range(0, dataset.shape[0], rows_per_chunk):
            yield dataset[start : start + rows_per_chunk, :]
        return

    raise ValueError("Unsupported matrix dimensionality")


def require_raw_count_matrix(handle: h5py.File, label: str) -> None:
    dataset = matrix_value_dataset(handle, label)
    if not np.issubdtype(dataset.dtype, np.number):
        raise ValueError(f"{label} /X must contain numeric raw count values")
    if np.issubdtype(dataset.dtype, np.unsignedinteger):
        return

    for chunk in iter_dataset_chunks(dataset):
        values = np.asarray(chunk)
        if values.size == 0:
            continue
        if not np.isfinite(values).all():
            raise ValueError(f"{label} /X contains non-finite values; expected raw non-negative integer counts")
        if (values < 0).any():
            raise ValueError(f"{label} /X contains negative values; expected raw non-negative integer counts")
        if np.issubdtype(values.dtype, np.integer):
            continue
        if (np.abs(values - np.rint(values)) > INTEGER_TOLERANCE).any():
            raise ValueError(
                f"{label} /X does not look like raw counts: found non-integer values. "
                "Provide raw non-negative integer counts in adata.X, not normalized or log-transformed expression."
            )


try:
    from anndata.io import sparse_dataset
except ImportError:  # anndata < 0.11
    from anndata.experimental import sparse_dataset


CELLTYPIST_TARGET_SUM = 1e4
CELLTYPIST_SUM_TOLERANCE = 1e-2
CELLTYPIST_SAMPLE_BLOCKS = 5
CELLTYPIST_SAMPLE_BLOCK_ROWS = 200
# Locations checked, in order, for the matrix CellTypist is given.
CELLTYPIST_MATRIX_LOCATIONS = ("X", "layers/counts", "raw/X")


def sample_row_indices(n_rows: int) -> np.ndarray:
    starts = np.linspace(0, max(0, n_rows - CELLTYPIST_SAMPLE_BLOCK_ROWS), CELLTYPIST_SAMPLE_BLOCKS)
    rows = {
        row
        for start in starts
        for row in range(int(start), min(n_rows, int(start) + CELLTYPIST_SAMPLE_BLOCK_ROWS))
    }
    return np.array(sorted(rows), dtype=np.int64)


def sample_csc_values(node: h5py.Group) -> np.ndarray:
    """A few contiguous runs of stored values from a CSC matrix.

    Per-cell totals need every gene, and a CSC matrix stores genes contiguously, so cells
    cannot be sampled without reading the whole file. The stored values are enough to tell
    counts from normalised expression, so only those are read.
    """
    data = node["data"]
    run = CELLTYPIST_SAMPLE_BLOCKS * CELLTYPIST_SAMPLE_BLOCK_ROWS * 100
    starts = np.linspace(0, max(0, data.shape[0] - run), CELLTYPIST_SAMPLE_BLOCKS)
    return np.concatenate([data[int(start) : int(start) + run] for start in starts])


def is_csc(node) -> bool:
    encoding = node.attrs.get("encoding-type", "csr_matrix")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    return encoding == "csc_matrix"


def sample_matrix_rows(handle: h5py.File, location: str) -> np.ndarray | None:
    """Return a dense sample of rows from a matrix in an h5ad, or None if it is absent."""
    if location not in handle:
        return None

    node = handle[location]
    if isinstance(node, h5py.Dataset):
        rows = sample_row_indices(node.shape[0])
        return np.vstack([np.asarray(node[start : stop]) for start, stop in _contiguous_runs(rows)])
    if isinstance(node, h5py.Group) and "indptr" in node and not is_csc(node):
        matrix = sparse_dataset(node)
        rows = sample_row_indices(matrix.shape[0])
        return np.vstack([matrix[start:stop].toarray() for start, stop in _contiguous_runs(rows)])
    return None


def _contiguous_runs(rows: np.ndarray):
    breaks = np.flatnonzero(np.diff(rows) != 1) + 1
    for run in np.split(rows, breaks):
        yield int(run[0]), int(run[-1]) + 1


def looks_like_log1p_normalised(values: np.ndarray) -> bool:
    """True when expm1 of every sampled cell sums to the CellTypist target of 10,000."""
    if values.size == 0 or not np.isfinite(values).all() or (values < 0).any():
        return False
    totals = np.expm1(values.astype(np.float64)).sum(axis=1)
    return bool(np.allclose(totals, CELLTYPIST_TARGET_SUM, rtol=CELLTYPIST_SUM_TOLERANCE))


def looks_like_raw_counts(values: np.ndarray) -> bool:
    if values.size == 0 or not np.isfinite(values).all() or (values < 0).any():
        return False
    return bool((np.abs(values - np.rint(values)) <= INTEGER_TOLERANCE).all())


def classify_celltypist_matrix(handle: h5py.File, location: str) -> str | None:
    """Return 'normalised', 'counts', or None for the matrix at a location.

    Dense and CSR matrices are classified from a sample of cells, with the exact check that
    expm1 sums to 10,000 per cell. CSC matrices are classified from stored values alone,
    where normalised means non-integer and no larger than log1p(10,000); the exact per-cell
    check for those is left to the annotation step, which has the matrix in memory.
    """
    if location not in handle:
        return None
    node = handle[location]

    if isinstance(node, h5py.Group) and "indptr" in node and is_csc(node):
        values = sample_csc_values(node)
        if values.size == 0 or not np.isfinite(values).all() or (values < 0).any():
            return None
        if looks_like_raw_counts(values):
            return "counts"
        if values.max() <= np.log1p(CELLTYPIST_TARGET_SUM) + 1e-3:
            return "normalised"
        return None

    sample = sample_matrix_rows(handle, location)
    if sample is None:
        return None
    if looks_like_log1p_normalised(sample):
        return "normalised"
    if looks_like_raw_counts(sample):
        return "counts"
    return None


def select_celltypist_matrix(handle: h5py.File, label: str) -> tuple[str, bool]:
    """Pick the matrix CellTypist should use and say whether it is already normalised.

    CellTypist needs log1p expression normalised to 10,000 counts per cell. X is used as
    it is when it already has that form. Otherwise raw counts are taken from X,
    layers['counts'], or raw.X, in that order, and normalised by the caller. Returns the
    HDF5 location of the matrix and whether it needs no further normalisation.
    """
    x_kind = classify_celltypist_matrix(handle, "X")
    if x_kind == "normalised":
        return "X", True

    for location in CELLTYPIST_MATRIX_LOCATIONS:
        kind = x_kind if location == "X" else classify_celltypist_matrix(handle, location)
        if kind == "counts":
            return location, False

    raise ValueError(
        f"{label} has no usable expression matrix for CellTypist. CellTypist needs log1p "
        "expression normalised to 10,000 counts per cell, or raw integer counts to "
        "normalise. adata.X is neither, and no raw counts were found in "
        "adata.layers['counts'] or adata.raw.X."
    )


def log1p_totals_are_normalised(matrix, rows_per_chunk: int = 100_000) -> bool:
    """Exact check, over every cell, that expm1 of a loaded matrix sums to 10,000 per cell."""
    from scipy import sparse

    for start in range(0, matrix.shape[0], rows_per_chunk):
        chunk = matrix[start : start + rows_per_chunk]
        if sparse.issparse(chunk):
            chunk = chunk.copy()
            chunk.data = np.expm1(chunk.data.astype(np.float64))
            totals = np.asarray(chunk.sum(axis=1)).ravel()
        else:
            totals = np.expm1(np.asarray(chunk, dtype=np.float64)).sum(axis=1)
        if not np.allclose(totals, CELLTYPIST_TARGET_SUM, rtol=CELLTYPIST_SUM_TOLERANCE):
            return False
    return True


def require_obs_column(adata: ad.AnnData, column: str, label: str) -> None:
    if column not in adata.obs:
        raise ValueError(f"Missing required obs column '{column}' in {label}")


def require_unique_index(adata: ad.AnnData, label: str) -> None:
    if not adata.obs_names.is_unique:
        raise ValueError(f"{label} cell indices are not unique")


def require_unique_var_index(adata: ad.AnnData, label: str) -> None:
    if not adata.var_names.is_unique:
        raise ValueError(f"{label} gene indices are not unique")


def infer_gene_id_type(var_names) -> str:
    names = [str(name) for name in var_names if str(name)]
    if not names:
        raise ValueError("Gene identifiers are empty")

    sample = names[: min(len(names), 1000)]
    ensembl_pattern = re.compile(r"^ENS[A-Z]*G\d+(?:\.\d+)?$")
    ensembl_hits = sum(bool(ensembl_pattern.match(name)) for name in sample)

    if ensembl_hits >= max(1, int(0.8 * len(sample))):
        return "ensembl"
    return "symbol"


def parse_filter_values(raw_values: str | None) -> list[str]:
    if not raw_values:
        return []
    values = [value.strip() for value in raw_values.split(",")]
    values = [value for value in values if value]
    if not values:
        raise ValueError("No valid reference filter values were provided")
    return values


def parse_comma_separated(raw_values: str | None) -> list[str]:
    if not raw_values:
        return []
    return [value.strip() for value in raw_values.split(",") if value.strip()]


def normalize_covariate_keys(value, label: str) -> list[str]:
    if value is None or value == "":
        return []
    if isinstance(value, str):
        return parse_comma_separated(value)
    if isinstance(value, (list, tuple, set)):
        keys = []
        for item in value:
            if not isinstance(item, str):
                raise ValueError(f"{label} must contain column-name strings, not nested values")
            if item.strip():
                keys.append(item.strip())
        return keys
    raise ValueError(f"{label} must be a string or list of column-name strings")


def single_value_from_config(section: str, key: str, value):
    if isinstance(value, (list, tuple, set)):
        values = list(value)
        if len(values) != 1:
            raise ValueError(
                f"{section}.{key} must contain exactly one value for nf-autoannotate, "
                f"not a sweep with {len(values)} values"
            )
        return values[0]
    if isinstance(value, dict):
        raise ValueError(f"{section}.{key} must be a single value, not a nested dictionary")
    return value


def load_nf_scvi_metrics_config(path: str | None) -> dict[str, dict]:
    if not path:
        return {}

    namespace = runpy.run_path(path)
    if "model_input" not in namespace:
        raise ValueError("--nf-scvi-metrics-config must define a model_input dictionary")

    raw_model_input = namespace.get("model_input") or {}
    raw_param_input = namespace.get("param_input") or {}
    if not isinstance(raw_model_input, dict):
        raise ValueError("model_input in --nf-scvi-metrics-config must be a dictionary")
    if not isinstance(raw_param_input, dict):
        raise ValueError("param_input in --nf-scvi-metrics-config must be a dictionary")

    allowed_model_keys = {
        "adata_path",  # nf-scVI-metrics dataset selector; ignored by nf-autoannotate.
        "counts_layer",
        "n_hidden",
        "n_latent",
        "n_layers",
        "dropout_rate",
        "dispersion",
        "gene_likelihood",
        "latent_distribution",
        "max_epochs",
        "accelerator",
        "devices",
        "train_size",
        "validation_size",
        "batch_size",
        "early_stopping",
        "early_stopping_patience",
    }
    allowed_param_keys = {
        "layer",
        "batch_key",
        "categorical_covariate_keys",
        "continuous_covariate_keys",
    }
    optional_sections = {
        "train_input": {
            "max_epochs",
            "accelerator",
            "devices",
            "train_size",
            "validation_size",
            "batch_size",
            "early_stopping",
            "early_stopping_patience",
        },
        "scanvi_train_input": {
            "max_epochs",
            "accelerator",
            "devices",
            "train_size",
            "validation_size",
            "batch_size",
            "early_stopping",
            "early_stopping_patience",
        },
        "query_train_input": {
            "max_epochs",
            "accelerator",
            "devices",
            "train_size",
            "validation_size",
            "batch_size",
            "early_stopping",
            "early_stopping_patience",
            "check_val_every_n_epoch",
            "plan_weight_decay",
        },
    }

    unknown_model = sorted(set(raw_model_input) - allowed_model_keys)
    if unknown_model:
        raise ValueError(
            f"Unsupported model_input keys in --nf-scvi-metrics-config: {', '.join(unknown_model)}"
        )
    unknown_param = sorted(set(raw_param_input) - allowed_param_keys)
    if unknown_param:
        raise ValueError(
            f"Unsupported param_input keys in --nf-scvi-metrics-config: {', '.join(unknown_param)}"
        )

    model_input = {
        key: single_value_from_config("model_input", key, value)
        for key, value in raw_model_input.items()
        if key != "adata_path"
    }
    param_input = {}
    for key, value in raw_param_input.items():
        if key in {"categorical_covariate_keys", "continuous_covariate_keys"}:
            param_input[key] = normalize_covariate_keys(value, f"param_input.{key}")
        else:
            param_input[key] = single_value_from_config("param_input", key, value)

    if "counts_layer" in model_input and "layer" not in param_input:
        param_input["layer"] = model_input["counts_layer"]

    for section, allowed_keys in optional_sections.items():
        raw_values = namespace.get(section, {}) or {}
        if not isinstance(raw_values, dict):
            raise ValueError(f"{section} in --nf-scvi-metrics-config must be a dictionary")
        unknown = sorted(set(raw_values) - allowed_keys)
        if unknown:
            raise ValueError(
                f"Unsupported {section} keys in --nf-scvi-metrics-config: {', '.join(unknown)}"
            )
        for key, value in raw_values.items():
            single_value_from_config(section, key, value)

    return {"param_input": param_input}


def resolve_setup_keys(args, nf_scvi_metrics_config: dict) -> tuple:
    """Resolve layer, batch key, and covariate keys from either the config file or CLI flags."""
    config_param_input = nf_scvi_metrics_config.get("param_input", {})
    layer = config_param_input.get("layer") if nf_scvi_metrics_config else None
    batch_key = config_param_input.get("batch_key") if nf_scvi_metrics_config else args.batch_key
    categorical_covariate_keys = (
        normalize_covariate_keys(
            config_param_input.get("categorical_covariate_keys"),
            "param_input.categorical_covariate_keys",
        )
        if nf_scvi_metrics_config
        else parse_comma_separated(args.categorical_covariate_keys)
    )
    continuous_covariate_keys = (
        normalize_covariate_keys(
            config_param_input.get("continuous_covariate_keys"),
            "param_input.continuous_covariate_keys",
        )
        if nf_scvi_metrics_config
        else parse_comma_separated(args.continuous_covariate_keys)
    )

    if isinstance(layer, (list, tuple, set, dict)):
        raise ValueError("param_input.layer must be a single layer name")
    if isinstance(batch_key, (list, tuple, set, dict)):
        raise ValueError("param_input.batch_key must be a single column name")
    if batch_key and batch_key in categorical_covariate_keys:
        raise ValueError("categorical covariate keys should contain additional columns, not the batch key")

    return layer, batch_key, categorical_covariate_keys, continuous_covariate_keys


def read_gene_list(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def write_gene_list(path: str, genes: list[str]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(genes))
        handle.write("\n")
