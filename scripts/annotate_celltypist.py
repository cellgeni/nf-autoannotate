#!/usr/bin/env python3

import argparse
import sys

import scanpy as sc
import celltypist
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from validation_common import log1p_totals_are_normalised, select_celltypist_matrix

DEFAULT_MIN_MODEL_GENE_OVERLAP_FRACTION = 0.8


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate query cells using a CellTypist model.")
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--model-path", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--min-model-gene-overlap",
        type=float,
        default=DEFAULT_MIN_MODEL_GENE_OVERLAP_FRACTION,
        help="Minimum fraction of model genes that must be present in the query.",
    )
    args = parser.parse_args()
    if not 0 <= args.min_model_gene_overlap <= 1:
        parser.error("--min-model-gene-overlap must be between 0 and 1")
    return args


def align_adata_to_gene_list(adata: sc.AnnData, genes: list[str]) -> sc.AnnData:
    gene_to_index = {gene: idx for idx, gene in enumerate(adata.var_names)}
    matrix = adata.X.tocsc() if sparse.issparse(adata.X) else sparse.csc_matrix(adata.X)

    target_positions = []
    source_columns = []
    for target_idx, gene in enumerate(genes):
        source_idx = gene_to_index.get(gene)
        if source_idx is not None:
            target_positions.append(target_idx)
            source_columns.append(source_idx)

    selected = matrix[:, source_columns].tocoo()
    remapped_columns = np.asarray(target_positions, dtype=int)[selected.col]
    aligned_matrix = sparse.csr_matrix(
        (selected.data, (selected.row, remapped_columns)),
        shape=(adata.n_obs, len(genes)),
        dtype=matrix.dtype,
    )
    return sc.AnnData(X=aligned_matrix, obs=adata.obs.copy(), var=pd.DataFrame(index=genes))


def confidence_for_labels(probability_matrix: pd.DataFrame, labels: pd.Series) -> pd.Series:
    if probability_matrix.empty:
        return pd.Series(index=labels.index, dtype=float)

    aligned_labels = labels.reindex(probability_matrix.index).to_numpy()
    label_to_column = {label: idx for idx, label in enumerate(probability_matrix.columns)}
    column_indices = np.array([label_to_column.get(label, -1) for label in aligned_labels])

    values = probability_matrix.to_numpy()
    row_indices = np.arange(values.shape[0])
    label_scores = values[row_indices, np.clip(column_indices, 0, values.shape[1] - 1)]
    scores = np.where(column_indices >= 0, label_scores, values.max(axis=1))
    return pd.Series(scores, index=probability_matrix.index, dtype=float)


def main():
    args = parse_args()

    with h5py.File(args.query_h5ad, "r") as handle:
        location, already_normalised = select_celltypist_matrix(handle, f"query h5ad '{args.query_h5ad}'")

    query = sc.read_h5ad(args.query_h5ad)
    if location == "layers/counts":
        query.X = query.layers["counts"]
    elif location == "raw/X":
        query = query.raw.to_adata()
    query.var_names_make_unique()
    if already_normalised and not log1p_totals_are_normalised(query.X):
        raise ValueError(
            "adata.X looks log1p-normalised, but expm1 of the cells does not sum to 10,000. "
            "Provide X normalised to 10,000 counts per cell, or raw counts in "
            "adata.layers['counts'] or adata.raw.X."
        )
    if already_normalised:
        print("CellTypist input: adata.X is already log1p-normalised to 10,000 counts per cell", file=sys.stderr)
    else:
        print(f"CellTypist input: normalising raw counts from {location}", file=sys.stderr)
        sc.pp.normalize_total(query, target_sum=1e4)
        sc.pp.log1p(query)

    model = celltypist.models.Model.load(args.model_path)
    model_genes = list(model.classifier.features)
    query_gene_set = set(query.var_names)
    overlapping_genes = [gene for gene in model_genes if gene in query_gene_set]
    overlap_fraction = len(overlapping_genes) / len(model_genes) if model_genes else 0.0
    if not overlapping_genes:
        raise ValueError("No genes overlap between the query and the CellTypist model")
    if overlap_fraction < args.min_model_gene_overlap:
        raise ValueError(
            "Only "
            f"{len(overlapping_genes)}/{len(model_genes)} CellTypist model genes "
            f"({overlap_fraction:.1%}) are present in the query; expected at least "
            f"{args.min_model_gene_overlap:.0%}. Check that both use the same gene identifiers, "
            "and whether genes were filtered out of the query. To annotate anyway, lower "
            "--celltypist_min_model_gene_overlap."
        )
    print(
        "CellTypist model gene overlap: "
        f"{len(overlapping_genes)}/{len(model_genes)} ({overlap_fraction:.1%})",
        file=sys.stderr,
    )
    if overlap_fraction < DEFAULT_MIN_MODEL_GENE_OVERLAP_FRACTION:
        print(
            f"WARNING: annotating with {overlap_fraction:.1%} of CellTypist model genes present. "
            "The missing genes are zero-filled, which degrades predictions.",
            file=sys.stderr,
        )

    aligned_query = align_adata_to_gene_list(query, model_genes)
    # Runs on CPU. CellTypist's use_GPU only affects majority-voting over-clustering, and its GPU
    # neighbour-graph construction calls scanpy's filter_genes on GPU arrays, which fails.
    predictions = celltypist.annotate(aligned_query, model=args.model_path, majority_voting=True)

    label_frame = predictions.predicted_labels
    if "predicted_labels" not in label_frame.columns:
        raise ValueError("CellTypist output is missing the 'predicted_labels' column")

    predicted_labels = label_frame["predicted_labels"].astype(str)
    if "majority_voting" in label_frame.columns:
        majority_labels = label_frame["majority_voting"].astype(str)
    else:
        majority_labels = predicted_labels

    probabilities = predictions.probability_matrix
    scores = confidence_for_labels(probabilities, majority_labels)

    output = pd.DataFrame(
        {
            "cell_id": aligned_query.obs_names.to_list(),
            "celltypist_predicted_label": predicted_labels.to_numpy(),
            "celltypist_majority_voting": majority_labels.to_numpy(),
            "celltypist_confidence": scores.reindex(aligned_query.obs_names).to_numpy(),
        }
    )
    output.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
