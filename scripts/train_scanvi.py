#!/usr/bin/env python3

import argparse
import json
import runpy

import anndata as ad
import scanpy as sc
from scvi.model import SCANVI, SCVI


def parse_args():
    parser = argparse.ArgumentParser(description="Train an scANVI reference model from a reference h5ad file.")
    parser.add_argument("--ref-h5ad", required=True)
    parser.add_argument("--reference-manifest", required=True)
    parser.add_argument("--ref-label-col", required=True)
    parser.add_argument("--batch-key")
    parser.add_argument("--categorical-covariate-keys")
    parser.add_argument("--continuous-covariate-keys")
    parser.add_argument("--layer")
    parser.add_argument("--unknown-label", required=True)
    parser.add_argument("--scvi-model-dir")
    parser.add_argument("--hvg-n-top-genes", type=int, default=2500)
    parser.add_argument("--training-config")
    parser.add_argument("--scvi-n-hidden", type=int)
    parser.add_argument("--scvi-n-latent", type=int)
    parser.add_argument("--scvi-n-layers", type=int)
    parser.add_argument("--scvi-dropout-rate", type=float)
    parser.add_argument("--scvi-dispersion")
    parser.add_argument("--scvi-gene-likelihood")
    parser.add_argument("--scvi-latent-distribution")
    parser.add_argument("--scvi-train-max-epochs", type=int)
    parser.add_argument("--scvi-train-batch-size", type=int)
    parser.add_argument("--scvi-train-train-size", type=float)
    parser.add_argument("--scvi-train-validation-size", type=float)
    parser.add_argument("--scvi-train-early-stopping", type=parse_bool)
    parser.add_argument("--scvi-train-early-stopping-patience", type=int)
    parser.add_argument("--scanvi-train-max-epochs", type=int)
    parser.add_argument("--scanvi-train-batch-size", type=int)
    parser.add_argument("--scanvi-train-train-size", type=float)
    parser.add_argument("--scanvi-train-validation-size", type=float)
    parser.add_argument("--scanvi-train-early-stopping", type=parse_bool)
    parser.add_argument("--scanvi-train-early-stopping-patience", type=int)
    parser.add_argument("--output-scvi-model-dir", required=True)
    parser.add_argument("--output-scvi-metadata", required=True)
    parser.add_argument("--output-model-dir", required=True)
    parser.add_argument("--output-metadata", required=True)
    return parser.parse_args()


def parse_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    normalized = value.strip().lower()
    if normalized in {"1", "true", "t", "yes", "y"}:
        return True
    if normalized in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected a boolean value, got {value}")


def apply_reference_filter(adata: ad.AnnData, manifest: dict) -> ad.AnnData:
    filter_col = manifest.get("ref_filter_col")
    filter_values = manifest.get("ref_filter_values") or []
    if not filter_col or not filter_values:
        return adata
    mask = adata.obs[filter_col].astype(str).isin(filter_values)
    filtered = adata[mask].copy()
    if filtered.n_obs == 0:
        raise ValueError("Reference filtering removed all cells before scANVI training")
    return filtered


def parse_comma_separated(raw_values: str | None) -> list[str]:
    if not raw_values:
        return []
    return [value.strip() for value in raw_values.split(",") if value.strip()]


def single_value_from_config(section: str, key: str, value):
    if isinstance(value, (list, tuple, set)):
        values = list(value)
        if len(values) != 1:
            raise ValueError(
                f"{section}.{key} must contain exactly one value for nf-autoannotate, "
                f"not a sweep with {len(values)} values"
            )
        return values[0]
    if isinstance(value, dict):
        raise ValueError(f"{section}.{key} must be a single value, not a nested dictionary")
    return value


def normalize_config_values(section: str, values: dict, allowed_keys: set[str]) -> dict:
    unknown = sorted(set(values) - allowed_keys)
    if unknown:
        raise ValueError(f"Unsupported {section} keys in --training-config: {', '.join(unknown)}")
    return {key: single_value_from_config(section, key, value) for key, value in values.items()}


def normalize_setup_config(values: dict) -> dict:
    allowed_keys = {
        "layer",
        "batch_key",
        "categorical_covariate_keys",
        "continuous_covariate_keys",
    }
    unknown = sorted(set(values) - allowed_keys)
    if unknown:
        raise ValueError(f"Unsupported param_input keys in --training-config: {', '.join(unknown)}")

    normalized = {}
    for key, value in values.items():
        if key in {"categorical_covariate_keys", "continuous_covariate_keys"}:
            normalized[key] = parse_comma_separated(value) if isinstance(value, str) else list(value or [])
        else:
            normalized[key] = single_value_from_config("param_input", key, value)
    return normalized


def load_training_config(path: str | None) -> dict[str, dict]:
    if not path:
        return {
            "model_input": {},
            "param_input": {},
            "train_input": {},
            "scanvi_train_input": {},
        }

    namespace = runpy.run_path(path)
    if "model_input" not in namespace:
        raise ValueError("--training-config must define a model_input dictionary")

    raw_model_input = namespace.get("model_input") or {}
    raw_param_input = namespace.get("param_input") or {}
    if not isinstance(raw_model_input, dict):
        raise ValueError("model_input in --training-config must be a dictionary")
    if not isinstance(raw_param_input, dict):
        raise ValueError("param_input in --training-config must be a dictionary")

    model_keys = {
        "n_hidden",
        "n_latent",
        "n_layers",
        "dropout_rate",
        "dispersion",
        "gene_likelihood",
        "latent_distribution",
    }
    train_keys = {
        "max_epochs",
        "accelerator",
        "devices",
        "train_size",
        "validation_size",
        "batch_size",
        "early_stopping",
        "early_stopping_patience",
    }
    nf_scvi_metrics_model_keys = model_keys | train_keys | {"adata_path", "counts_layer"}
    unknown_model = sorted(set(raw_model_input) - nf_scvi_metrics_model_keys)
    if unknown_model:
        raise ValueError(f"Unsupported model_input keys in --training-config: {', '.join(unknown_model)}")

    normalized_model_input = normalize_config_values(
        "model_input",
        {key: value for key, value in raw_model_input.items() if key not in {"adata_path", "counts_layer"}},
        model_keys | train_keys,
    )
    model_input = {key: value for key, value in normalized_model_input.items() if key in model_keys}
    train_input = {key: value for key, value in normalized_model_input.items() if key in train_keys}
    param_input = normalize_setup_config(raw_param_input)
    if "counts_layer" in raw_model_input and "layer" not in param_input:
        param_input["layer"] = single_value_from_config(
            "model_input", "counts_layer", raw_model_input["counts_layer"]
        )

    explicit_train_input = namespace.get("train_input", {}) or {}
    if explicit_train_input:
        if not isinstance(explicit_train_input, dict):
            raise ValueError("train_input in --training-config must be a dictionary")
        train_input.update(normalize_config_values("train_input", explicit_train_input, train_keys))

    raw_scanvi_train_input = namespace.get("scanvi_train_input", None)
    if raw_scanvi_train_input is None:
        scanvi_train_input = dict(train_input)
    else:
        if not isinstance(raw_scanvi_train_input, dict):
            raise ValueError("scanvi_train_input in --training-config must be a dictionary")
        scanvi_train_input = normalize_config_values(
            "scanvi_train_input", raw_scanvi_train_input, train_keys
        )

    return {
        "model_input": model_input,
        "param_input": param_input,
        "train_input": train_input,
        "scanvi_train_input": scanvi_train_input,
    }


def add_if_set(target: dict, key: str, value):
    if value is not None and value != "":
        target[key] = value


def training_config_from_args(args: argparse.Namespace) -> dict[str, dict]:
    categorical_keys = parse_comma_separated(args.categorical_covariate_keys)
    continuous_keys = parse_comma_separated(args.continuous_covariate_keys)

    param_input = {}
    add_if_set(param_input, "layer", args.layer)
    add_if_set(param_input, "batch_key", args.batch_key)
    if categorical_keys:
        param_input["categorical_covariate_keys"] = categorical_keys
    if continuous_keys:
        param_input["continuous_covariate_keys"] = continuous_keys

    model_input = {}
    add_if_set(model_input, "n_hidden", args.scvi_n_hidden)
    add_if_set(model_input, "n_latent", args.scvi_n_latent)
    add_if_set(model_input, "n_layers", args.scvi_n_layers)
    add_if_set(model_input, "dropout_rate", args.scvi_dropout_rate)
    add_if_set(model_input, "dispersion", args.scvi_dispersion)
    add_if_set(model_input, "gene_likelihood", args.scvi_gene_likelihood)
    add_if_set(model_input, "latent_distribution", args.scvi_latent_distribution)

    train_input = {}
    add_if_set(train_input, "max_epochs", args.scvi_train_max_epochs)
    add_if_set(train_input, "batch_size", args.scvi_train_batch_size)
    add_if_set(train_input, "train_size", args.scvi_train_train_size)
    add_if_set(train_input, "validation_size", args.scvi_train_validation_size)
    add_if_set(train_input, "early_stopping", args.scvi_train_early_stopping)
    add_if_set(train_input, "early_stopping_patience", args.scvi_train_early_stopping_patience)

    scanvi_train_input = {}
    add_if_set(scanvi_train_input, "max_epochs", args.scanvi_train_max_epochs)
    add_if_set(scanvi_train_input, "batch_size", args.scanvi_train_batch_size)
    add_if_set(scanvi_train_input, "train_size", args.scanvi_train_train_size)
    add_if_set(scanvi_train_input, "validation_size", args.scanvi_train_validation_size)
    add_if_set(scanvi_train_input, "early_stopping", args.scanvi_train_early_stopping)
    add_if_set(scanvi_train_input, "early_stopping_patience", args.scanvi_train_early_stopping_patience)

    return {
        "model_input": model_input,
        "param_input": param_input,
        "train_input": train_input,
        "scanvi_train_input": scanvi_train_input or dict(train_input),
    }


def normalize_setup_params(setup_params: dict) -> dict:
    setup_params = dict(setup_params)
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
    if args.training_config:
        training_config = load_training_config(args.training_config)
    else:
        training_config = training_config_from_args(args)
    setup_params = normalize_setup_params(training_config["param_input"])
    effective_batch_key = setup_params.get("batch_key")

    with open(args.reference_manifest, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    ref = ad.read_h5ad(args.ref_h5ad)
    ref.var_names_make_unique()
    ref = apply_reference_filter(ref, manifest)

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
        # HVG selection runs over the full reference gene space so the model stays
        # independent of any particular query.
        ref.obs["_scanvi_label"] = ref.obs[args.ref_label_col].astype(str)
        ref = subset_to_hvgs(ref, args.hvg_n_top_genes, effective_batch_key)
        SCVI.setup_anndata(ref, **setup_params)
        vae = SCVI(ref, **training_config["model_input"])
        vae.train(**training_config["train_input"])

        scanvi = SCANVI.from_scvi_model(
            vae,
            adata=ref,
            labels_key="_scanvi_label",
            unlabeled_category=args.unknown_label,
        )
        scanvi.train(**training_config["scanvi_train_input"])
        genes = ref.var_names.tolist()
        inference_mode = "query_mapping"

    vae.save(args.output_scvi_model_dir, overwrite=True, save_anndata=True)
    scanvi.save(args.output_model_dir, overwrite=True, save_anndata=True)

    metadata = {
        "ref_label_col": args.ref_label_col,
        "batch_key": effective_batch_key,
        "unknown_label": args.unknown_label,
        "categorical_covariate_keys": setup_params.get("categorical_covariate_keys") or [],
        "continuous_covariate_keys": setup_params.get("continuous_covariate_keys") or [],
        "layer": setup_params.get("layer"),
        "hvg_n_top_genes": args.hvg_n_top_genes,
        "training_config": training_config,
        "ref_filter_col": manifest.get("ref_filter_col"),
        "ref_filter_values": manifest.get("ref_filter_values", []),
        "scvi_model_dir": args.scvi_model_dir,
        "inference_mode": inference_mode,
        "genes": genes,
    }
    scvi_metadata = {
        **metadata,
        "model_type": "scvi",
    }
    scanvi_metadata = {
        **metadata,
        "model_type": "scanvi",
        "base_model_dir": args.output_scvi_model_dir,
    }
    with open(args.output_scvi_metadata, "w", encoding="utf-8") as handle:
        json.dump(scvi_metadata, handle, indent=2)
    with open(f"{args.output_scvi_model_dir}/scvi_model_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(scvi_metadata, handle, indent=2)
    with open(args.output_metadata, "w", encoding="utf-8") as handle:
        json.dump(scanvi_metadata, handle, indent=2)
    with open(f"{args.output_model_dir}/scanvi_model_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(scanvi_metadata, handle, indent=2)


if __name__ == "__main__":
    main()
