#!/usr/bin/env python3

import argparse
import json
import re
import runpy
from pathlib import Path

import anndata as ad


def parse_args():
    parser = argparse.ArgumentParser(description="Validate reference/query h5ad inputs.")
    parser.add_argument("--ref-h5ad", required=True)
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--ref-label-col", required=True)
    parser.add_argument("--batch-key")
    parser.add_argument("--categorical-covariate-keys")
    parser.add_argument("--continuous-covariate-keys")
    parser.add_argument("--nf-scvi-metrics-config")
    parser.add_argument("--panhumanpy-feature-names-col")
    parser.add_argument("--ref-filter-col")
    parser.add_argument("--ref-filter-values")
    parser.add_argument("--output", required=True)
    parser.add_argument("--shared-genes-output", required=True)
    return parser.parse_args()


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


def resolve_reference_filter(args: argparse.Namespace) -> tuple[str | None, list[str]]:
    filter_col = args.ref_filter_col
    filter_values = parse_filter_values(args.ref_filter_values)
    if (filter_col and not filter_values) or (not filter_col and filter_values):
        raise ValueError("Both --ref-filter-col and --ref-filter-values are required together")
    return filter_col, filter_values


def main():
    args = parse_args()

    ref_path = Path(args.ref_h5ad)
    query_path = Path(args.query_h5ad)

    if not ref_path.exists():
        raise FileNotFoundError(f"Reference h5ad not found: {ref_path}")
    if not query_path.exists():
        raise FileNotFoundError(f"Query h5ad not found: {query_path}")

    ref = ad.read_h5ad(ref_path, backed="r")
    query = ad.read_h5ad(query_path, backed="r")

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

    if args.panhumanpy_feature_names_col and args.panhumanpy_feature_names_col not in query.var:
        raise ValueError(
            f"Missing required var column '{args.panhumanpy_feature_names_col}' in query h5ad"
        )

    require_unique_index(ref, "Reference")
    require_unique_index(query, "Query")
    require_unique_var_index(ref, "Reference")
    require_unique_var_index(query, "Query")

    ref_filter_col, ref_filter_values = resolve_reference_filter(args)
    filtered_ref_cell_count = int(ref.n_obs)
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
        "ref_h5ad": str(ref_path.resolve()),
        "query_h5ad": str(query_path.resolve()),
        "ref_label_col": args.ref_label_col,
        "batch_key": batch_key,
        "categorical_covariate_keys": categorical_covariate_keys,
        "continuous_covariate_keys": continuous_covariate_keys,
        "panhumanpy_feature_names_col": args.panhumanpy_feature_names_col,
        "ref_filter_col": ref_filter_col,
        "ref_filter_values": ref_filter_values,
        "gene_id_type": ref_gene_id_type,
        "shared_gene_count": len(shared_genes),
        "ref_cell_count": int(ref.n_obs),
        "filtered_ref_cell_count": filtered_ref_cell_count,
        "query_cell_count": int(query.n_obs),
    }

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    with open(args.shared_genes_output, "w", encoding="utf-8") as handle:
        handle.write("\n".join(shared_genes))
        handle.write("\n")


if __name__ == "__main__":
    main()
