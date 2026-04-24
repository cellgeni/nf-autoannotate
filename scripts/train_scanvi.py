#!/usr/bin/env python3

import argparse
import json

import anndata as ad
import pandas as pd
from scvi.model import SCANVI, SCVI


def parse_args():
    parser = argparse.ArgumentParser(description="Train an scANVI model from reference and query h5ad files.")
    parser.add_argument("--ref-h5ad", required=True)
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--validation-manifest", required=True)
    parser.add_argument("--shared-genes", required=True)
    parser.add_argument("--ref-label-col", required=True)
    parser.add_argument("--batch-key", required=True)
    parser.add_argument("--unknown-label", required=True)
    parser.add_argument("--scvi-model-dir")
    parser.add_argument("--output-model-dir", required=True)
    parser.add_argument("--output-metadata", required=True)
    return parser.parse_args()


def apply_reference_timepoint_filter(adata: ad.AnnData, manifest: dict) -> ad.AnnData:
    timepoint_col = manifest.get("ref_timepoint_col")
    timepoints = manifest.get("ref_timepoints") or []
    if not timepoint_col or not timepoints:
        return adata
    mask = adata.obs[timepoint_col].astype(str).isin(timepoints)
    filtered = adata[mask].copy()
    if filtered.n_obs == 0:
        raise ValueError("Reference timepoint filtering removed all cells before scANVI training")
    return filtered


def read_gene_list(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def main():
    args = parse_args()

    with open(args.validation_manifest, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    ref = ad.read_h5ad(args.ref_h5ad)
    query = ad.read_h5ad(args.query_h5ad)
    ref.var_names_make_unique()
    query.var_names_make_unique()
    ref = apply_reference_timepoint_filter(ref, manifest)

    if args.scvi_model_dir:
        ref.obs["_scanvi_label"] = ref.obs[args.ref_label_col].astype(str)
        vae = SCVI.load(args.scvi_model_dir, adata=ref)
        scanvi = SCANVI.from_scvi_model(
            vae,
            adata=ref,
            labels_key="_scanvi_label",
            unlabeled_category=args.unknown_label,
        )
        scanvi.train()
        genes = ref.var_names.tolist()
        inference_mode = "query_mapping"
    else:
        shared_genes = read_gene_list(args.shared_genes)
        ref = ref[:, shared_genes].copy()
        query = query[:, shared_genes].copy()

        ref.obs["_scanvi_label"] = ref.obs[args.ref_label_col].astype(str)
        query.obs["_scanvi_label"] = pd.Categorical([args.unknown_label] * query.n_obs)

        combined = ad.concat(
            {"reference": ref, "query": query},
            join="outer",
            label="_source",
            index_unique=None,
            merge="same",
        )
        combined.obs["_scanvi_label"] = combined.obs["_scanvi_label"].astype(str)

        SCVI.setup_anndata(combined, batch_key=args.batch_key)
        vae = SCVI(combined)
        vae.train()

        scanvi = SCANVI.from_scvi_model(
            vae,
            adata=combined,
            labels_key="_scanvi_label",
            unlabeled_category=args.unknown_label,
        )
        scanvi.train()
        genes = shared_genes
        inference_mode = "direct_load"

    scanvi.save(args.output_model_dir, overwrite=True, save_anndata=True)

    metadata = {
        "ref_label_col": args.ref_label_col,
        "batch_key": args.batch_key,
        "unknown_label": args.unknown_label,
        "ref_timepoint_col": manifest.get("ref_timepoint_col"),
        "ref_timepoints": manifest.get("ref_timepoints", []),
        "scvi_model_dir": args.scvi_model_dir,
        "inference_mode": inference_mode,
        "genes": genes,
    }
    with open(args.output_metadata, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    with open(f"{args.output_model_dir}/scanvi_model_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)


if __name__ == "__main__":
    main()
