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
    parser.add_argument("--batch-key", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def require_obs_column(adata: ad.AnnData, column: str, label: str) -> None:
    if column not in adata.obs:
        raise ValueError(f"Missing required obs column '{column}' in {label}")


def require_unique_index(adata: ad.AnnData, label: str) -> None:
    if not adata.obs_names.is_unique:
        raise ValueError(f"{label} cell indices are not unique")


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
    require_obs_column(ref, args.batch_key, "reference h5ad")
    require_obs_column(query, args.batch_key, "query h5ad")
    require_unique_index(ref, "Reference")
    require_unique_index(query, "Query")

    ref.var_names_make_unique()
    query.var_names_make_unique()

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
        "gene_id_type": ref_gene_id_type,
        "shared_gene_count": len(shared_genes),
        "shared_genes": shared_genes,
        "ref_cell_count": int(ref.n_obs),
        "query_cell_count": int(query.n_obs),
    }

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


if __name__ == "__main__":
    main()
