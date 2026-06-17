#!/usr/bin/env python3

import argparse
import json

import scanpy as sc
import celltypist
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Train a CellTypist model from a reference h5ad.")
    parser.add_argument("--ref-h5ad", required=True)
    parser.add_argument("--validation-manifest", required=True)
    parser.add_argument("--shared-genes", required=True)
    parser.add_argument("--ref-label-col", required=True)
    parser.add_argument("--output-model", required=True)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--training-mode", choices=["standard", "detailed"], default="detailed")
    parser.add_argument("--balance-cells-per-label", type=int, default=500)
    parser.add_argument("--feature-selection", action="store_true")
    parser.add_argument("--feature-selection-top-genes", type=int, default=100)
    parser.add_argument("--feature-selection-max-iter", type=int, default=5)
    parser.add_argument("--max-iter", type=int, default=100)
    return parser.parse_args()


def balanced_training_subset(adata: sc.AnnData, label_col: str, cells_per_label: int) -> tuple[sc.AnnData, int]:
    if cells_per_label <= 0:
        return adata, adata.n_obs

    sampled_index = celltypist.samples.downsample_adata(
        adata,
        mode="each",
        n_cells=cells_per_label,
        by=label_col,
        return_index=True,
    )
    sampled = adata[sampled_index].copy()
    return sampled, sampled.n_obs


def select_features(
    adata: sc.AnnData,
    label_col: str,
    n_jobs: int,
    top_genes: int,
    max_iter: int,
) -> tuple[sc.AnnData, list[str]]:
    top_genes = max(1, min(top_genes, adata.n_vars))

    rough_model = celltypist.train(
        adata,
        labels=label_col,
        n_jobs=n_jobs,
        max_iter=max_iter,
        use_SGD=True,
    )

    gene_index = np.argpartition(
        np.abs(rough_model.classifier.coef_),
        -top_genes,
        axis=1,
    )[:, -top_genes:]
    selected = np.unique(gene_index)
    selected_genes = adata.var_names[selected].tolist()
    return adata[:, selected_genes].copy(), selected_genes


def apply_reference_filter(adata: sc.AnnData, manifest: dict) -> sc.AnnData:
    filter_col = manifest.get("ref_filter_col") or manifest.get("ref_timepoint_col")
    filter_values = manifest.get("ref_filter_values") or manifest.get("ref_timepoints") or []
    if not filter_col or not filter_values:
        return adata
    mask = adata.obs[filter_col].astype(str).isin(filter_values)
    filtered = adata[mask].copy()
    if filtered.n_obs == 0:
        raise ValueError("Reference filtering removed all cells before CellTypist training")
    return filtered


def read_gene_list(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def main():
    args = parse_args()

    with open(args.validation_manifest, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    ref = sc.read_h5ad(args.ref_h5ad)
    ref.var_names_make_unique()
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)
    ref = apply_reference_filter(ref, manifest)
    shared_genes = read_gene_list(args.shared_genes)
    ref = ref[:, shared_genes].copy()

    training_adata = ref
    selected_genes = ref.var_names.tolist()
    balanced_cell_count = ref.n_obs
    feature_selection_used = False

    if args.training_mode == "detailed":
        balanced_ref, balanced_cell_count = balanced_training_subset(
            ref,
            args.ref_label_col,
            args.balance_cells_per_label,
        )
        training_adata = balanced_ref
        selected_genes = balanced_ref.var_names.tolist()
        feature_selection_used = args.feature_selection

        if args.feature_selection:
            training_adata, selected_genes = select_features(
                balanced_ref,
                args.ref_label_col,
                args.n_jobs,
                args.feature_selection_top_genes,
                args.feature_selection_max_iter,
            )

    model = celltypist.train(
        training_adata,
        labels=args.ref_label_col,
        check_expression=not feature_selection_used,
        n_jobs=args.n_jobs,
        max_iter=args.max_iter,
    )
    model.write(args.output_model)

    metadata = {
        "ref_label_col": args.ref_label_col,
        "training_mode": args.training_mode,
        "ref_filter_col": manifest.get("ref_filter_col") or manifest.get("ref_timepoint_col"),
        "ref_filter_values": manifest.get("ref_filter_values") or manifest.get("ref_timepoints", []),
        "shared_gene_count": len(shared_genes),
        "balanced_cell_count": balanced_cell_count,
        "feature_selection": feature_selection_used,
        "feature_selection_top_genes": args.feature_selection_top_genes if feature_selection_used else None,
        "feature_selection_max_iter": args.feature_selection_max_iter if feature_selection_used else None,
        "final_training_max_iter": args.max_iter,
        "balance_cells_per_label": args.balance_cells_per_label if args.training_mode == "detailed" else None,
        "genes": selected_genes,
    }
    with open(args.output_metadata, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


if __name__ == "__main__":
    main()
