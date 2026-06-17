#!/usr/bin/env python3

import argparse
import json
import runpy

import anndata as ad
import pandas as pd
import scanpy as sc
from scvi.model import SCANVI, SCVI


def parse_args():
    parser = argparse.ArgumentParser(description="Train an scANVI model from reference and query h5ad files.")
    parser.add_argument("--ref-h5ad", required=True)
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--validation-manifest", required=True)
    parser.add_argument("--shared-genes", required=True)
    parser.add_argument("--ref-label-col", required=True)
    parser.add_argument("--batch-key")
    parser.add_argument("--categorical-covariate-keys")
    parser.add_argument("--unknown-label", required=True)
    parser.add_argument("--scvi-model-dir")
    parser.add_argument("--hvg-n-top-genes", type=int, default=2500)
    parser.add_argument("--training-config")
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


def parse_comma_separated(raw_values: str | None) -> list[str]:
    if not raw_values:
        return []
    return [value.strip() for value in raw_values.split(",") if value.strip()]


def load_training_config(path: str | None) -> dict[str, dict]:
    if not path:
        return {
            "model_input": {},
            "param_input": {},
            "train_input": {},
            "scanvi_train_input": {},
        }

    namespace = runpy.run_path(path)
    config = {
        "model_input": dict(namespace.get("model_input", {})),
        "param_input": dict(namespace.get("param_input", {})),
        "train_input": dict(namespace.get("train_input", {})),
        "scanvi_train_input": dict(namespace.get("scanvi_train_input", namespace.get("train_input", {}))),
    }

    allowed_keys = {
        "model_input": {
            "n_hidden",
            "n_latent",
            "n_layers",
            "dropout_rate",
            "dispersion",
            "gene_likelihood",
            "latent_distribution",
        },
        "param_input": {
            "layer",
            "batch_key",
            "categorical_covariate_keys",
            "continuous_covariate_keys",
        },
        "train_input": {
            "max_epochs",
            "accelerator",
            "devices",
            "train_size",
            "validation_size",
            "batch_size",
            "early_stopping",
            "early_stopping_patience",
        },
        "scanvi_train_input": {
            "max_epochs",
            "accelerator",
            "devices",
            "train_size",
            "validation_size",
            "batch_size",
            "early_stopping",
            "early_stopping_patience",
        },
    }
    list_allowed = {
        ("param_input", "categorical_covariate_keys"),
        ("param_input", "continuous_covariate_keys"),
    }

    for section, values in config.items():
        unknown = sorted(set(values) - allowed_keys[section])
        if unknown:
            raise ValueError(f"Unsupported {section} keys in --training-config: {', '.join(unknown)}")
        for key, value in values.items():
            if isinstance(value, (list, tuple, set)) and (section, key) not in list_allowed:
                raise ValueError(
                    f"{section}.{key} must be a single value, not a list/tuple/set sweep. "
                    "Use one value per parameter for nf-autoannotate."
                )
            if isinstance(value, dict):
                raise ValueError(f"{section}.{key} must be a single value, not a nested dictionary")

    return config


def merge_setup_params(config_params: dict, batch_key: str | None, categorical_keys: list[str]) -> dict:
    setup_params = dict(config_params)
    if batch_key:
        setup_params["batch_key"] = batch_key
    if categorical_keys:
        setup_params["categorical_covariate_keys"] = categorical_keys
    if setup_params.get("categorical_covariate_keys") == []:
        setup_params["categorical_covariate_keys"] = None
    if setup_params.get("continuous_covariate_keys") == []:
        setup_params["continuous_covariate_keys"] = None
    return setup_params


def subset_to_hvgs(adata: ad.AnnData, n_top_genes: int, batch_key: str | None) -> ad.AnnData:
    if n_top_genes <= 0 or adata.n_vars <= n_top_genes:
        return adata

    hvg_batch_key = batch_key if batch_key and batch_key in adata.obs else None
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key=hvg_batch_key,
        flavor="seurat_v3",
        subset=False,
    )
    hvg_mask = adata.var["highly_variable"].fillna(False).to_numpy()
    if not hvg_mask.any():
        raise ValueError("HVG selection did not identify any genes for scVI training")
    return adata[:, hvg_mask].copy()


def main():
    args = parse_args()
    if args.hvg_n_top_genes < 0:
        raise ValueError("--hvg-n-top-genes must be a non-negative integer")
    training_config = load_training_config(args.training_config)
    effective_batch_key = args.batch_key or training_config["param_input"].get("batch_key")

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
        scanvi.train(**training_config["scanvi_train_input"])
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
        combined = subset_to_hvgs(combined, args.hvg_n_top_genes, effective_batch_key)
        categorical_covariate_keys = parse_comma_separated(args.categorical_covariate_keys)
        setup_params = merge_setup_params(
            training_config["param_input"],
            effective_batch_key,
            categorical_covariate_keys,
        )

        SCVI.setup_anndata(combined, **setup_params)
        vae = SCVI(combined, **training_config["model_input"])
        vae.train(**training_config["train_input"])

        scanvi = SCANVI.from_scvi_model(
            vae,
            adata=combined,
            labels_key="_scanvi_label",
            unlabeled_category=args.unknown_label,
        )
        scanvi.train(**training_config["scanvi_train_input"])
        genes = combined.var_names.tolist()
        inference_mode = "direct_load"

    scanvi.save(args.output_model_dir, overwrite=True, save_anndata=True)

    metadata = {
        "ref_label_col": args.ref_label_col,
        "batch_key": effective_batch_key,
        "unknown_label": args.unknown_label,
        "categorical_covariate_keys": parse_comma_separated(args.categorical_covariate_keys),
        "hvg_n_top_genes": args.hvg_n_top_genes,
        "training_config": training_config,
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
