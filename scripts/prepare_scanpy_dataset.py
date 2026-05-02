#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen

import numpy as np
import scanpy as sc
from sklearn.model_selection import train_test_split


PBMC68K_URL = "https://ndownloader.figshare.com/files/27686886"


DATASET_LOADERS = {
    "pbmc68k_reduced": sc.datasets.pbmc68k_reduced,
    "pbmc3k": sc.datasets.pbmc3k,
    "pbmc68k": lambda: sc.read("pbmc68k.h5ad", backup_url=PBMC68K_URL),
}

DEFAULT_STRATIFY_CANDIDATES = [
    "cell_type",
    "celltype",
    "bulk_labels",
    "louvain",
    "leiden",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load an AnnData dataset and split it into reference/query h5ad files. "
            "Default split is 80%% reference and 20%% query."
        )
    )
    source_group = parser.add_mutually_exclusive_group()
    source_group.add_argument("--dataset", choices=sorted(DATASET_LOADERS.keys()))
    source_group.add_argument("--input-h5ad")
    source_group.add_argument("--input-url")
    parser.add_argument("--query-fraction", type=float, default=0.2)
    parser.add_argument("--stratify-col")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    if not any([args.dataset, args.input_h5ad, args.input_url]):
        args.dataset = "pbmc68k"
    return args


def infer_download_filename(url: str, fallback: str) -> str:
    parsed = urlparse(url)
    candidate = Path(parsed.path).name
    if candidate and "." in candidate:
        return candidate
    return fallback


def download_to_outdir(url: str, outdir: Path, filename_override: str | None = None) -> Path:
    raw_dir = outdir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    filename = filename_override or infer_download_filename(url, fallback="downloaded_input.h5ad")
    destination = raw_dir / filename

    print(f"Downloading dataset from: {url}")
    with urlopen(url) as response, open(destination, "wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)

    print(f"Downloaded file to: {destination}")
    return destination


def load_dataset(args: argparse.Namespace, outdir: Path):
    if args.dataset:
        print(f"Loading built-in Scanpy dataset: {args.dataset}")
        return DATASET_LOADERS[args.dataset]().copy(), "scanpy_dataset", args.dataset

    if args.input_h5ad:
        input_path = Path(args.input_h5ad)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")
        print(f"Loading local AnnData file: {input_path}")
        return sc.read(input_path).copy(), "local_file", str(input_path)

    downloaded = download_to_outdir(args.input_url, outdir)
    print(f"Loading downloaded AnnData file: {downloaded}")
    return sc.read(downloaded).copy(), "remote_url", args.input_url


def choose_stratify_column(adata: sc.AnnData, requested: str | None) -> str | None:
    if requested:
        if requested not in adata.obs.columns:
            raise ValueError(f"Requested --stratify-col '{requested}' not found in adata.obs")
        return requested

    for candidate in DEFAULT_STRATIFY_CANDIDATES:
        if candidate in adata.obs.columns and adata.obs[candidate].nunique() > 1:
            return candidate

    return None


def maybe_stratify_labels(adata: sc.AnnData, column: str | None):
    if not column:
        return None

    labels = adata.obs[column].astype(str)
    label_counts = labels.value_counts()
    if (label_counts < 2).any():
        print(
            "Stratified splitting skipped because at least one label occurs only once. "
            "Falling back to random split."
        )
        return None

    return labels


def main() -> None:
    args = parse_args()

    if not 0.0 < args.query_fraction < 1.0:
        raise ValueError("--query-fraction must be strictly between 0 and 1")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata, source_type, source_value = load_dataset(args, outdir)
    adata.var_names_make_unique()

    stratify_col = choose_stratify_column(adata, args.stratify_col)
    stratify_labels = maybe_stratify_labels(adata, stratify_col)

    all_indices = np.arange(adata.n_obs)
    ref_indices, query_indices = train_test_split(
        all_indices,
        test_size=args.query_fraction,
        random_state=args.seed,
        shuffle=True,
        stratify=stratify_labels,
    )

    ref_adata = adata[ref_indices].copy()
    query_adata = adata[query_indices].copy()

    ref_path = outdir / "reference.h5ad"
    query_path = outdir / "query.h5ad"
    metadata_path = outdir / "dataset_split_metadata.json"

    ref_adata.write_h5ad(ref_path)
    query_adata.write_h5ad(query_path)

    metadata = {
        "source_type": source_type,
        "source": source_value,
        "dataset": args.dataset,
        "n_cells_total": int(adata.n_obs),
        "n_genes_total": int(adata.n_vars),
        "n_cells_reference": int(ref_adata.n_obs),
        "n_cells_query": int(query_adata.n_obs),
        "query_fraction": float(args.query_fraction),
        "reference_fraction": float(1 - args.query_fraction),
        "seed": int(args.seed),
        "stratify_col": stratify_col,
        "stratified": stratify_labels is not None,
        "reference_h5ad": str(ref_path),
        "query_h5ad": str(query_path),
    }
    with open(metadata_path, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    print(f"Reference saved to: {ref_path}")
    print(f"Query saved to: {query_path}")
    print(f"Metadata saved to: {metadata_path}")


if __name__ == "__main__":
    main()
