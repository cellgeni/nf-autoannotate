#!/usr/bin/env python3

import argparse
import json

import scanpy as sc
import celltypist
import pandas as pd
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate query cells using a CellTypist model.")
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--validation-manifest", required=True)
    parser.add_argument("--model-path", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def align_adata_to_gene_list(adata: sc.AnnData, genes: list[str]) -> sc.AnnData:
    gene_to_index = {gene: idx for idx, gene in enumerate(adata.var_names)}
    matrix = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
    aligned_blocks = []

    for gene in genes:
        if gene in gene_to_index:
            aligned_blocks.append(matrix[:, gene_to_index[gene]])
        else:
            aligned_blocks.append(sparse.csr_matrix((adata.n_obs, 1), dtype=matrix.dtype))

    aligned_matrix = sparse.hstack(aligned_blocks, format="csr")
    return sc.AnnData(X=aligned_matrix, obs=adata.obs.copy(), var=pd.DataFrame(index=genes))


def confidence_for_labels(probability_matrix: pd.DataFrame, labels: pd.Series) -> pd.Series:
    if probability_matrix.empty:
        return pd.Series(index=labels.index, dtype=float)

    scores = []
    for cell_id, row in probability_matrix.iterrows():
        label = labels.loc[cell_id]
        scores.append(row[label] if label in row.index else row.max())
    return pd.Series(scores, index=probability_matrix.index, dtype=float)


def main():
    args = parse_args()

    with open(args.validation_manifest, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    query = sc.read_h5ad(args.query_h5ad)
    query.var_names_make_unique()
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)

    model = celltypist.models.Model.load(args.model_path)
    model_genes = list(model.classifier.features)
    shared_set = set(manifest["shared_genes"])
    if not any(gene in shared_set for gene in model_genes):
        raise ValueError("No genes overlap between query input and CellTypist model features")

    aligned_query = align_adata_to_gene_list(query, model_genes)
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
