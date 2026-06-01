#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scvi.model import SCANVI


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate query cells using a trained scANVI model.")
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--shared-genes", required=True)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--label-column", required=True)
    parser.add_argument("--score-column", required=True)
    parser.add_argument("--query-max-epochs", type=int, default=100)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def align_adata_to_gene_list(adata: ad.AnnData, genes: list[str]) -> ad.AnnData:
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
    return ad.AnnData(X=aligned_matrix, obs=adata.obs.copy(), var=pd.DataFrame(index=genes))


def read_gene_list(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def load_gene_list(model_dir: str, shared_genes_path: str) -> list[str]:
    metadata_path = Path(model_dir) / "scanvi_model_metadata.json"
    if metadata_path.exists():
        with open(metadata_path, "r", encoding="utf-8") as handle:
            return json.load(handle)["genes"]
    return read_gene_list(shared_genes_path)


def ensure_scanvi_label_column(adata: ad.AnnData, metadata: dict) -> ad.AnnData:
    unknown_label = metadata.get("unknown_label")
    if not unknown_label:
        return adata

    prepared = adata.copy()
    prepared.obs["_scanvi_label"] = pd.Categorical([unknown_label] * prepared.n_obs)
    return prepared


def main():
    args = parse_args()

    query = ad.read_h5ad(args.query_h5ad)
    query.var_names_make_unique()

    metadata_path = Path(args.model_dir) / "scanvi_model_metadata.json"
    metadata = {}
    if metadata_path.exists():
        with open(metadata_path, "r", encoding="utf-8") as handle:
            metadata = json.load(handle)

    if metadata.get("inference_mode") == "query_mapping":
        SCANVI.prepare_query_anndata(query, args.model_dir)
        scanvi = SCANVI.load_query_data(query, args.model_dir)
        scanvi.train(
            max_epochs=args.query_max_epochs,
            plan_kwargs=dict(weight_decay=0.0),
            check_val_every_n_epoch=10,
        )
        active_query = query
        labels = scanvi.predict()
        soft = scanvi.predict(soft=True)
    else:
        model_genes = load_gene_list(args.model_dir, args.shared_genes)
        active_query = align_adata_to_gene_list(query, model_genes)
        active_query = ensure_scanvi_label_column(active_query, metadata)
        scanvi = SCANVI.load(args.model_dir, adata=active_query)
        labels = scanvi.predict(active_query)
        soft = scanvi.predict(active_query, soft=True)

    if isinstance(soft, pd.DataFrame):
        scores = soft.max(axis=1).to_numpy()
    else:
        scores = np.asarray(soft).max(axis=1)

    output = pd.DataFrame(
        {
            "cell_id": active_query.obs_names.to_list(),
            args.label_column: pd.Series(labels, index=active_query.obs_names).astype(str).to_numpy(),
            args.score_column: scores,
        }
    )
    output.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
