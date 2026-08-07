#!/usr/bin/env python3

import argparse
import json
import runpy
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scvi.model import SCANVI


MIN_MODEL_GENE_OVERLAP_FRACTION = 0.8


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate query cells using a trained scANVI model.")
    parser.add_argument("--query-h5ad", required=True)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--label-column", required=True)
    parser.add_argument("--score-column", required=True)
    parser.add_argument("--query-max-epochs", type=int, default=100)
    parser.add_argument("--query-batch-size", type=int)
    parser.add_argument("--query-train-size", type=float)
    parser.add_argument("--query-validation-size", type=float)
    parser.add_argument("--query-early-stopping", type=parse_bool)
    parser.add_argument("--query-early-stopping-patience", type=int)
    parser.add_argument("--query-check-val-every-n-epoch", type=int, default=10)
    parser.add_argument("--query-plan-weight-decay", type=float, default=0.0)
    parser.add_argument("--training-config")
    parser.add_argument("--output", required=True)
    parser.add_argument("--latent-output", required=True)
    return parser.parse_args()


def parse_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    normalized = value.strip().lower()
    if normalized in {"1", "true", "t", "yes", "y"}:
        return True
    if normalized in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected a boolean value, got {value!r}")


def add_if_set(target: dict, key: str, value):
    if value is not None and value != "":
        target[key] = value


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


def normalize_single_value_config(section: str, values: dict, allowed_keys: set[str]) -> dict:
    unknown = sorted(set(values) - allowed_keys)
    if unknown:
        raise ValueError(f"Unsupported {section} keys in --training-config: {', '.join(unknown)}")
    return {key: single_value_from_config(section, key, value) for key, value in values.items()}


def query_train_input_from_args(args: argparse.Namespace) -> dict:
    query_train_input = {
        "max_epochs": args.query_max_epochs,
        "plan_kwargs": {"weight_decay": args.query_plan_weight_decay},
        "check_val_every_n_epoch": args.query_check_val_every_n_epoch,
    }
    add_if_set(query_train_input, "batch_size", args.query_batch_size)
    add_if_set(query_train_input, "train_size", args.query_train_size)
    add_if_set(query_train_input, "validation_size", args.query_validation_size)
    add_if_set(query_train_input, "early_stopping", args.query_early_stopping)
    add_if_set(query_train_input, "early_stopping_patience", args.query_early_stopping_patience)
    return query_train_input


def load_query_train_input(path: str | None, args: argparse.Namespace) -> dict:
    if not path:
        return query_train_input_from_args(args)

    namespace = runpy.run_path(path)
    raw_query_train_input = namespace.get("query_train_input", {}) or {}
    if not isinstance(raw_query_train_input, dict):
        raise ValueError("query_train_input in --training-config must be a dictionary")
    allowed_keys = {
        "max_epochs",
        "accelerator",
        "devices",
        "train_size",
        "validation_size",
        "batch_size",
        "early_stopping",
        "early_stopping_patience",
        "check_val_every_n_epoch",
        "plan_weight_decay",
    }
    query_train_input = normalize_single_value_config(
        "query_train_input", raw_query_train_input, allowed_keys
    )
    if "plan_weight_decay" in query_train_input:
        query_train_input["plan_kwargs"] = {"weight_decay": query_train_input.pop("plan_weight_decay")}
    if query_train_input:
        return query_train_input
    return query_train_input_from_args(
        argparse.Namespace(
            query_max_epochs=100,
            query_batch_size=None,
            query_train_size=None,
            query_validation_size=None,
            query_early_stopping=None,
            query_early_stopping_patience=None,
            query_check_val_every_n_epoch=10,
            query_plan_weight_decay=0.0,
        )
    )


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


def load_gene_list(model_dir: str, query: ad.AnnData) -> list[str]:
    """Return the gene space the model was trained on.

    Only the model itself is authoritative. A reference gene list is not a valid
    substitute: scVI keeps `hvg_n_top_genes` highly variable genes and CellTypist keeps a
    per-label subset, so a model's registry is a strict subset of the reference, and
    aligning a query to the wrong list silently scrambles the feature columns.

    Metadata written by this pipeline is preferred because it avoids loading the model,
    and scvi-tools' own registry is the fallback for externally trained directories.
    """
    metadata_path = Path(model_dir) / "scanvi_model_metadata.json"
    if metadata_path.exists():
        with open(metadata_path, "r", encoding="utf-8") as handle:
            metadata = json.load(handle)
        genes = metadata.get("genes")
        if genes:
            return list(genes)

    try:
        var_names = SCANVI.prepare_query_anndata(
            query, model_dir, return_reference_var_names=True
        )
    except Exception as error:  # noqa: BLE001 - surfaced with actionable context below
        raise ValueError(
            f"Could not determine the gene space of the scANVI model in {model_dir}. It has no "
            "scanvi_model_metadata.json from this pipeline, and reading the scvi-tools registry "
            f"failed: {error}"
        ) from error
    return [str(name) for name in var_names]


def report_gene_coverage(query: ad.AnnData, model_genes: list[str]) -> None:
    if not model_genes:
        return

    query_gene_set = set(query.var_names)
    covered = sum(1 for gene in model_genes if gene in query_gene_set)
    fraction = covered / len(model_genes)
    print(
        f"scANVI model gene overlap: {covered}/{len(model_genes)} ({fraction:.1%})",
        file=sys.stderr,
    )
    if covered == 0:
        raise ValueError("No genes overlap between the query and the scANVI model gene space")
    if fraction < MIN_MODEL_GENE_OVERLAP_FRACTION:
        # Missing genes are zero-filled, which reads as "not detected" rather than
        # "not measured" and degrades the latent representation.
        print(
            f"WARNING: only {fraction:.1%} of scANVI model genes are present in the query; "
            "the remainder are zero-filled and predictions may be unreliable.",
            file=sys.stderr,
        )


def ensure_scanvi_label_column(adata: ad.AnnData, metadata: dict) -> ad.AnnData:
    unknown_label = metadata.get("unknown_label")
    if not unknown_label:
        return adata

    prepared = adata.copy()
    prepared.obs["_scanvi_label"] = pd.Categorical([unknown_label] * prepared.n_obs)
    return prepared


def write_scanvi_latent_output(
    path: str,
    obs_names: pd.Index,
    latent: np.ndarray,
    metadata: dict,
) -> None:
    latent = np.asarray(latent)
    if latent.ndim != 2:
        raise ValueError(f"Expected two-dimensional latent representation, got shape {latent.shape}")
    if latent.shape[0] != len(obs_names):
        raise ValueError(
            "Latent representation row count does not match query cells: "
            f"{latent.shape[0]} latent rows for {len(obs_names)} query cells"
        )

    output = ad.AnnData(obs=pd.DataFrame(index=obs_names.copy()))
    output.obsm["X_scanvi"] = latent
    # scANVI uses the same latent-variable model family as scVI. Keep an
    # explicit compatibility key for downstream tools that look for X_scvi.
    output.obsm["X_scvi"] = latent.copy()
    output.uns["scanvi_latent"] = {
        "obsm_keys": ["X_scanvi", "X_scvi"],
        "source": "scvi-tools SCANVI.get_latent_representation",
        "inference_mode": metadata.get("inference_mode"),
        "x_scvi_note": "Compatibility alias of X_scanvi from the annotated scANVI model.",
    }
    output.write_h5ad(path)


def main():
    args = parse_args()

    query = ad.read_h5ad(args.query_h5ad)
    query.var_names_make_unique()

    metadata_path = Path(args.model_dir) / "scanvi_model_metadata.json"
    metadata = {}
    if metadata_path.exists():
        with open(metadata_path, "r", encoding="utf-8") as handle:
            metadata = json.load(handle)

    model_genes = load_gene_list(args.model_dir, query)
    report_gene_coverage(query, model_genes)

    if metadata.get("inference_mode") == "query_mapping":
        SCANVI.prepare_query_anndata(query, args.model_dir)
        scanvi = SCANVI.load_query_data(query, args.model_dir)
        scanvi.train(**load_query_train_input(args.training_config, args))
        active_query = query
        labels = scanvi.predict()
        soft = scanvi.predict(soft=True)
        latent = scanvi.get_latent_representation()
    else:
        active_query = align_adata_to_gene_list(query, model_genes)
        active_query = ensure_scanvi_label_column(active_query, metadata)
        scanvi = SCANVI.load(args.model_dir, adata=active_query)
        labels = scanvi.predict(active_query)
        soft = scanvi.predict(active_query, soft=True)
        latent = scanvi.get_latent_representation(active_query)

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
    write_scanvi_latent_output(args.latent_output, active_query.obs_names, latent, metadata)


if __name__ == "__main__":
    main()
