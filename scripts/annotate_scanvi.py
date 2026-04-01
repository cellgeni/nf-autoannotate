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
    parser.add_argument("--validation-manifest", required=True)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--label-column", required=True)
    parser.add_argument("--score-column", required=True)
    parser.add_argument("--query-max-epochs", type=int, default=100)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def align_adata_to_gene_list(adata: ad.AnnData, genes: list[str]) -> ad.AnnData:
    gene_to_index = {gene: idx for idx, gene in enumerate(adata.var_names)}
    matrix = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
    aligned_blocks = []

    for gene in genes:
        if gene in gene_to_index:
            aligned_blocks.append(matrix[:, gene_to_index[gene]])
        else:
            aligned_blocks.append(sparse.csr_matrix((adata.n_obs, 1), dtype=matrix.dtype))

    aligned_matrix = sparse.hstack(aligned_blocks, format="csr")
    return ad.AnnData(X=aligned_matrix, obs=adata.obs.copy(), var=pd.DataFrame(index=genes))


def load_gene_list(model_dir: str, manifest_genes: list[str]) -> list[str]:
    metadata_path = Path(model_dir) / "scanvi_model_metadata.json"
    if metadata_path.exists():
        with open(metadata_path, "r", encoding="utf-8") as handle:
            return json.load(handle)["genes"]
    return manifest_genes


def main():
    args = parse_args()

    with open(args.validation_manifest, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)

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
        model_genes = load_gene_list(args.model_dir, manifest["shared_genes"])
        active_query = align_adata_to_gene_list(query, model_genes)
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
