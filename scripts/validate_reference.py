#!/usr/bin/env python3

"""Validate the reference h5ad and publish the gene space used for model training.

Training is deliberately independent of any query object, so this runs once per
pipeline execution no matter how many queries are annotated.
"""

import argparse
import json
from pathlib import Path

import h5py

from validation_common import (
    infer_gene_id_type,
    load_nf_scvi_metrics_config,
    parse_filter_values,
    read_h5ad_metadata,
    read_obs_column,
    require_obs_column,
    require_unique_index,
    require_unique_var_index,
    resolve_setup_keys,
    write_gene_list,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Validate the reference h5ad input.")
    parser.add_argument("--ref-h5ad", required=True)
    parser.add_argument("--ref-label-col", required=True)
    parser.add_argument("--batch-key")
    parser.add_argument("--categorical-covariate-keys")
    parser.add_argument("--continuous-covariate-keys")
    parser.add_argument("--nf-scvi-metrics-config")
    parser.add_argument("--ref-filter-col")
    parser.add_argument("--ref-filter-values")
    parser.add_argument("--output", required=True)
    parser.add_argument("--genes-output", required=True)
    return parser.parse_args()


def resolve_reference_filter(args: argparse.Namespace) -> tuple[str | None, list[str]]:
    filter_col = args.ref_filter_col
    filter_values = parse_filter_values(args.ref_filter_values)
    if (filter_col and not filter_values) or (not filter_col and filter_values):
        raise ValueError("Both --ref-filter-col and --ref-filter-values are required together")
    return filter_col, filter_values


def main():
    args = parse_args()

    ref_path = Path(args.ref_h5ad)
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference h5ad not found: {ref_path}")

    with h5py.File(ref_path, "r") as ref_handle:
        ref = read_h5ad_metadata(ref_handle, "reference h5ad")

    require_obs_column(ref, args.ref_label_col, "reference h5ad")
    require_unique_index(ref, "Reference")
    require_unique_var_index(ref, "Reference")

    nf_scvi_metrics_config = load_nf_scvi_metrics_config(args.nf_scvi_metrics_config)
    layer, batch_key, categorical_covariate_keys, continuous_covariate_keys = resolve_setup_keys(
        args, nf_scvi_metrics_config
    )

    if layer and layer not in ref.layers:
        raise ValueError(f"Missing required layer '{layer}' in reference h5ad")
    if batch_key:
        require_obs_column(ref, batch_key, "reference h5ad")
    for covariate_key in categorical_covariate_keys + continuous_covariate_keys:
        require_obs_column(ref, covariate_key, "reference h5ad")

    ref_filter_col, ref_filter_values = resolve_reference_filter(args)
    filtered_ref_cell_count = int(ref.n_obs)
    if ref_filter_col:
        require_obs_column(ref, ref_filter_col, "reference h5ad")
        with h5py.File(ref_path, "r") as ref_handle:
            ref_filter_column = read_obs_column(ref_handle, ref_filter_col).astype(str)
        available_filter_values = sorted(set(ref_filter_column))
        missing_filter_values = sorted(set(ref_filter_values) - set(available_filter_values))
        if missing_filter_values:
            raise ValueError(
                "Requested reference filter values were not found in "
                f"'{ref_filter_col}': {', '.join(missing_filter_values)}"
            )
        filtered_ref_cell_count = int(ref_filter_column.isin(ref_filter_values).sum())
        if filtered_ref_cell_count == 0:
            raise ValueError("Reference filtering removed all cells")

    reference_genes = [str(gene) for gene in ref.var_names]
    gene_id_type = infer_gene_id_type(ref.var_names)

    manifest = {
        "ref_h5ad": str(ref_path.resolve()),
        "ref_label_col": args.ref_label_col,
        "batch_key": batch_key,
        "layer": layer,
        "categorical_covariate_keys": categorical_covariate_keys,
        "continuous_covariate_keys": continuous_covariate_keys,
        "ref_filter_col": ref_filter_col,
        "ref_filter_values": ref_filter_values,
        "gene_id_type": gene_id_type,
        "reference_gene_count": len(reference_genes),
        "ref_cell_count": int(ref.n_obs),
        "filtered_ref_cell_count": filtered_ref_cell_count,
    }

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    write_gene_list(args.genes_output, reference_genes)


if __name__ == "__main__":
    main()
