#!/usr/bin/env python3

import argparse
import json
import re
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
    if args.batch_key:
        require_obs_column(ref, args.batch_key, "reference h5ad")
        require_obs_column(query, args.batch_key, "query h5ad")

    categorical_covariate_keys = parse_comma_separated(args.categorical_covariate_keys)
    if args.batch_key and args.batch_key in categorical_covariate_keys:
        raise ValueError("--categorical-covariate-keys should contain additional columns, not the batch key")
    for covariate_key in categorical_covariate_keys:
        require_obs_column(ref, covariate_key, "reference h5ad")
        require_obs_column(query, covariate_key, "query h5ad")

    continuous_covariate_keys = parse_comma_separated(args.continuous_covariate_keys)
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
        "batch_key": args.batch_key,
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
