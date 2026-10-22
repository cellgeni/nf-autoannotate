#!/usr/bin/env python3

"""Validate one query h5ad and record how well it covers the reference gene space.

Runs once per query. Models are trained from the reference alone, so a query that
covers the reference poorly is reported here rather than silently degrading the
predictions later.
"""

import argparse
import json
import sys
from pathlib import Path

import h5py

from validation_common import (
    H5adMetadata,
    infer_gene_id_type,
    load_nf_scvi_metrics_config,
    parse_tools,
    read_gene_list,
    read_h5ad_metadata,
    require_obs_column,
    require_raw_count_matrix,
    select_celltypist_matrix,
    require_unique_index,
    require_unique_var_index,
    resolve_setup_keys,
)


LOW_COVERAGE_WARNING_FRACTION = 0.8


def parse_args():
    parser = argparse.ArgumentParser(description="Validate a single query h5ad input.")
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--tools", required=True)
    parser.add_argument("--reference-genes")
    parser.add_argument("--reference-manifest")
    parser.add_argument("--batch-key")
    parser.add_argument("--categorical-covariate-keys")
    parser.add_argument("--continuous-covariate-keys")
    parser.add_argument("--nf-scvi-metrics-config")
    parser.add_argument("--panhumanpy-feature-names-col")
    parser.add_argument("--query-cluster-col", required=True)
    parser.add_argument("--compute-missing-clusters", action="store_true")
    parser.add_argument("--skip-merge", action="store_true")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def require_query_cluster_column(
    query: H5adMetadata,
    query_cluster_col: str,
    compute_missing_clusters: bool,
    merge_annotations: bool,
) -> None:
    if not merge_annotations:
        # The cluster column is only read by the merge step, for the cluster summary and
        # the marker gene table. Demanding it from a run that publishes predictions alone
        # would reject queries over a column nothing is going to look at.
        return
    if compute_missing_clusters:
        return
    require_obs_column(query, query_cluster_col, "query h5ad")


def reference_coverage(query_genes: set[str], reference_genes: list[str]) -> tuple[int, float]:
    if not reference_genes:
        return 0, 0.0
    covered = sum(1 for gene in reference_genes if gene in query_genes)
    return covered, covered / len(reference_genes)


def main():
    args = parse_args()

    tools = parse_tools(args.tools)

    query_path = Path(args.query_h5ad)
    if not query_path.exists():
        raise FileNotFoundError(f"Query h5ad not found: {query_path}")

    with h5py.File(query_path, "r") as query_handle:
        query_label = f"query h5ad '{args.sample_id}'"
        if tools == {"celltypist"}:
            # CellTypist alone also accepts X that is already log1p-normalised, or raw
            # counts held in layers['counts'] or raw.X.
            select_celltypist_matrix(query_handle, query_label)
        else:
            require_raw_count_matrix(query_handle, query_label)
        query = read_h5ad_metadata(query_handle, query_label)

    require_unique_index(query, f"Query '{args.sample_id}'")
    require_unique_var_index(query, f"Query '{args.sample_id}'")

    nf_scvi_metrics_config = load_nf_scvi_metrics_config(args.nf_scvi_metrics_config)
    layer, batch_key, categorical_covariate_keys, continuous_covariate_keys = resolve_setup_keys(
        args, nf_scvi_metrics_config
    )

    if "scanvi" not in tools and (
        layer or batch_key or categorical_covariate_keys or continuous_covariate_keys
    ):
        raise ValueError("scVI/scANVI setup options were provided but 'scanvi' is not among the selected tools")

    if layer and layer not in query.layers:
        raise ValueError(f"Missing required layer '{layer}' in query h5ad '{args.sample_id}'")
    if batch_key:
        require_obs_column(query, batch_key, f"query h5ad '{args.sample_id}'")
    for covariate_key in categorical_covariate_keys + continuous_covariate_keys:
        require_obs_column(query, covariate_key, f"query h5ad '{args.sample_id}'")

    if args.panhumanpy_feature_names_col:
        if "panhumanpy" not in tools:
            raise ValueError(
                "--panhumanpy-feature-names-col was provided but 'panhumanpy' is not among the selected tools"
            )
        if args.panhumanpy_feature_names_col not in query.var_columns:
            raise ValueError(
                f"Missing required var column '{args.panhumanpy_feature_names_col}' "
                f"in query h5ad '{args.sample_id}'"
            )

    require_query_cluster_column(
        query,
        args.query_cluster_col,
        args.compute_missing_clusters,
        not args.skip_merge,
    )

    query_gene_id_type = infer_gene_id_type(query.var_names)
    query_genes = {str(gene) for gene in query.var_names}

    reference_gene_count = None
    covered_gene_count = None
    coverage_fraction = None

    if args.reference_genes:
        reference_genes = read_gene_list(args.reference_genes)
        reference_gene_count = len(reference_genes)
        covered_gene_count, coverage_fraction = reference_coverage(query_genes, reference_genes)

        print(
            f"Query '{args.sample_id}' reference gene coverage: "
            f"{covered_gene_count}/{reference_gene_count} ({coverage_fraction:.1%})",
            file=sys.stderr,
        )
        if coverage_fraction < LOW_COVERAGE_WARNING_FRACTION:
            print(
                f"WARNING: query '{args.sample_id}' covers only {coverage_fraction:.1%} of the reference "
                f"gene space. Genes absent from the query are zero-filled when aligning to trained models, "
                "which degrades predictions and may fail the CellTypist model-gene overlap check.",
                file=sys.stderr,
            )
        if covered_gene_count == 0:
            raise ValueError(
                f"Query '{args.sample_id}' shares no genes with the reference h5ad"
            )

    if args.reference_manifest:
        with open(args.reference_manifest, "r", encoding="utf-8") as handle:
            reference_manifest = json.load(handle)
        reference_gene_id_type = reference_manifest.get("gene_id_type")
        if reference_gene_id_type and reference_gene_id_type != query_gene_id_type:
            raise ValueError(
                "Reference and query use different gene identifier types: "
                f"reference={reference_gene_id_type}, query '{args.sample_id}'={query_gene_id_type}. "
                "Make sure both h5ad files use the same identifier scheme."
            )

    manifest = {
        "sample_id": args.sample_id,
        "tools": sorted(tools),
        "query_h5ad": str(query_path.resolve()),
        "batch_key": batch_key,
        "layer": layer,
        "categorical_covariate_keys": categorical_covariate_keys,
        "continuous_covariate_keys": continuous_covariate_keys,
        "panhumanpy_feature_names_col": args.panhumanpy_feature_names_col,
        "query_cluster_col": args.query_cluster_col,
        "compute_missing_clusters": args.compute_missing_clusters,
        "merge_annotations": not args.skip_merge,
        "gene_id_type": query_gene_id_type,
        "query_gene_count": len(query_genes),
        "query_cell_count": int(query.n_obs),
        "reference_gene_count": reference_gene_count,
        "reference_genes_covered": covered_gene_count,
        "reference_gene_coverage_fraction": coverage_fraction,
    }

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


if __name__ == "__main__":
    main()
