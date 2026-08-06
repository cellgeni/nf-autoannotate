#!/usr/bin/env python3

import argparse
import json
import re
import runpy

import h5py
import numpy as np
from pathlib import Path

import anndata as ad


MATRIX_CHUNK_SIZE = 1_000_000
INTEGER_TOLERANCE = 1e-6

TOOL_CHOICES = ("celltypist", "scanvi", "panhumanpy")
REFERENCE_TOOLS = frozenset({"celltypist", "scanvi"})


def parse_args():
    parser = argparse.ArgumentParser(description="Validate reference/query h5ad inputs.")
    parser.add_argument("--ref-h5ad")
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--ref-label-col")
    parser.add_argument("--tools", default=",".join(TOOL_CHOICES))
    parser.add_argument("--batch-key")
    parser.add_argument("--categorical-covariate-keys")
    parser.add_argument("--continuous-covariate-keys")
    parser.add_argument("--nf-scvi-metrics-config")
    parser.add_argument("--panhumanpy-feature-names-col")
    parser.add_argument("--query-cluster-col", required=True)
    parser.add_argument("--compute-missing-clusters", action="store_true")
    parser.add_argument("--ref-filter-col")
    parser.add_argument("--ref-filter-values")
    parser.add_argument("--output", required=True)
    parser.add_argument("--shared-genes-output", required=True)
    return parser.parse_args()


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


def resolve_reference_filter(args: argparse.Namespace) -> tuple[str | None, list[str]]:
    filter_col = args.ref_filter_col
    filter_values = parse_filter_values(args.ref_filter_values)
    if (filter_col and not filter_values) or (not filter_col and filter_values):
        raise ValueError("Both --ref-filter-col and --ref-filter-values are required together")
    return filter_col, filter_values


def require_query_cluster_column(
    query: ad.AnnData,
    query_cluster_col: str,
    compute_missing_clusters: bool,
) -> None:
    if compute_missing_clusters:
        return
    require_obs_column(query, query_cluster_col, "query h5ad")


def main():
    args = parse_args()

    tools = parse_tools(args.tools)
    needs_reference = bool(REFERENCE_TOOLS & tools)

    query_path = Path(args.query_h5ad)
    if not query_path.exists():
        raise FileNotFoundError(f"Query h5ad not found: {query_path}")

    ref_path = None
    if needs_reference:
        if not args.ref_h5ad:
            raise ValueError(
                "--ref-h5ad is required when the selected tools include "
                f"{', '.join(sorted(REFERENCE_TOOLS & tools))}"
            )
        if not args.ref_label_col:
            raise ValueError(
                "--ref-label-col is required when the selected tools include "
                f"{', '.join(sorted(REFERENCE_TOOLS & tools))}"
            )
        ref_path = Path(args.ref_h5ad)
        if not ref_path.exists():
            raise FileNotFoundError(f"Reference h5ad not found: {ref_path}")

    with h5py.File(query_path, "r") as query_handle:
        require_raw_count_matrix(query_handle, "query h5ad")

    ref = ad.read_h5ad(ref_path, backed="r") if needs_reference else None
    query = ad.read_h5ad(query_path, backed="r")

    if needs_reference:
        require_obs_column(ref, args.ref_label_col, "reference h5ad")

    nf_scvi_metrics_config = load_nf_scvi_metrics_config(args.nf_scvi_metrics_config)
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

    if "scanvi" not in tools and (
        layer or batch_key or categorical_covariate_keys or continuous_covariate_keys
    ):
        raise ValueError("scVI/scANVI setup options were provided but 'scanvi' is not among the selected tools")

    if isinstance(layer, (list, tuple, set, dict)):
        raise ValueError("param_input.layer must be a single layer name")
    if layer:
        if layer not in ref.layers:
            raise ValueError(f"Missing required layer '{layer}' in reference h5ad")
        if layer not in query.layers:
            raise ValueError(f"Missing required layer '{layer}' in query h5ad")

    if isinstance(batch_key, (list, tuple, set, dict)):
        raise ValueError("param_input.batch_key must be a single column name")
    if batch_key:
        require_obs_column(ref, batch_key, "reference h5ad")
        require_obs_column(query, batch_key, "query h5ad")

    if batch_key and batch_key in categorical_covariate_keys:
        raise ValueError("categorical covariate keys should contain additional columns, not the batch key")
    for covariate_key in categorical_covariate_keys:
        require_obs_column(ref, covariate_key, "reference h5ad")
        require_obs_column(query, covariate_key, "query h5ad")

    for covariate_key in continuous_covariate_keys:
        require_obs_column(ref, covariate_key, "reference h5ad")
        require_obs_column(query, covariate_key, "query h5ad")

    if args.panhumanpy_feature_names_col:
        if "panhumanpy" not in tools:
            raise ValueError(
                "--panhumanpy-feature-names-col was provided but 'panhumanpy' is not among the selected tools"
            )
        if args.panhumanpy_feature_names_col not in query.var:
            raise ValueError(
                f"Missing required var column '{args.panhumanpy_feature_names_col}' in query h5ad"
            )
    require_query_cluster_column(query, args.query_cluster_col, args.compute_missing_clusters)

    require_unique_index(query, "Query")
    require_unique_var_index(query, "Query")

    ref_filter_col, ref_filter_values = resolve_reference_filter(args)
    if ref_filter_col and not needs_reference:
        raise ValueError("Reference filtering was requested but no reference-based tool is selected")

    ref_gene_id_type = None
    shared_genes = None
    ref_cell_count = None
    filtered_ref_cell_count = None

    if needs_reference:
        require_unique_index(ref, "Reference")
        require_unique_var_index(ref, "Reference")

        ref_cell_count = int(ref.n_obs)
        filtered_ref_cell_count = ref_cell_count
        if ref_filter_col:
            require_obs_column(ref, ref_filter_col, "reference h5ad")
            available_filter_values = sorted({str(value) for value in ref.obs[ref_filter_col].astype(str)})
            missing_filter_values = sorted(set(ref_filter_values) - set(available_filter_values))
            if missing_filter_values:
                raise ValueError(
                    "Requested reference filter values were not found in "
                    f"'{ref_filter_col}': {', '.join(missing_filter_values)}"
                )
            filtered_ref_cell_count = int(ref.obs[ref_filter_col].astype(str).isin(ref_filter_values).sum())
            if filtered_ref_cell_count == 0:
                raise ValueError("Reference filtering removed all cells")

        ref_gene_id_type = infer_gene_id_type(ref.var_names)
        query_gene_id_type = infer_gene_id_type(query.var_names)
        if ref_gene_id_type != query_gene_id_type:
            raise ValueError(
                "Reference and query use different gene identifier types: "
                f"reference={ref_gene_id_type}, query={query_gene_id_type}. "
                "Make sure both h5ad files use the same identifier scheme."
            )

        query_genes = set(query.var_names)
        shared_genes = [gene for gene in ref.var_names if gene in query_genes]
        if not shared_genes:
            raise ValueError("Reference and query h5ad files do not share any genes")

    manifest = {
        "tools": sorted(tools),
        "ref_h5ad": str(ref_path.resolve()) if ref_path else None,
        "query_h5ad": str(query_path.resolve()),
        "ref_label_col": args.ref_label_col if needs_reference else None,
        "batch_key": batch_key,
        "categorical_covariate_keys": categorical_covariate_keys,
        "continuous_covariate_keys": continuous_covariate_keys,
        "panhumanpy_feature_names_col": args.panhumanpy_feature_names_col,
        "query_cluster_col": args.query_cluster_col,
        "compute_missing_clusters": args.compute_missing_clusters,
        "ref_filter_col": ref_filter_col,
        "ref_filter_values": ref_filter_values,
        "gene_id_type": ref_gene_id_type,
        "shared_gene_count": len(shared_genes) if shared_genes is not None else None,
        "ref_cell_count": ref_cell_count,
        "filtered_ref_cell_count": filtered_ref_cell_count,
        "query_cell_count": int(query.n_obs),
    }

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    # Downstream training/annotation processes only consume the shared-gene list when a
    # reference is involved, so a query-only run leaves the file unwritten.
    if shared_genes is not None:
        with open(args.shared_genes_output, "w", encoding="utf-8") as handle:
            handle.write("\n".join(shared_genes))
            handle.write("\n")


if __name__ == "__main__":
    main()
