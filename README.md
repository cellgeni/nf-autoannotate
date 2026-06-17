# nf-autoannotate

`nf-autoannotate` is a Nextflow pipeline for annotating a query `.h5ad` dataset with three independent methods:

- CellTypist
- scANVI
- PanHumanPy / Pan-Human Azimuth

The pipeline writes the method outputs back into the original query AnnData object and produces a cluster-level summary table. It does not create a consensus label.

## Workflow

The workflow is implemented in [main.nf](main.nf), with supporting Python scripts under [scripts/](scripts).

[prepare_scanpy_dataset.py](scripts/prepare_scanpy_dataset.py) is a utility script for creating local demo reference/query `.h5ad` files from a Scanpy dataset or an existing AnnData file.

The pipeline stages are:

1. Validate the reference and query AnnData files, then write `validation_manifest.json` and `shared_genes.txt`.
2. Train or load a CellTypist model.
3. Train or load an scANVI model.
4. Annotate the query with CellTypist.
5. Annotate the query with scANVI.
6. Annotate the query with PanHumanPy.
7. Merge method outputs into the query `.h5ad` and write a cluster annotation summary.

## Required Inputs

- `--ref_h5ad`: reference AnnData file.
- `--query_h5ad`: query AnnData file to annotate.
- `--project_tag`: output filename prefix.
- `--ref_label_col`: reference `obs` column containing cell type labels.

Optional reference filtering:

- `--ref_timepoint_col`: reference `obs` column containing timepoint labels.
- `--ref_timepoints`: comma-separated timepoints to keep, for example `E12,E14`.

`--ref_timepoint_col` and `--ref_timepoints` must be supplied together.

## Tool-Specific Parameters

The pipeline selects model behavior from artifact flags instead of explicit mode flags:

- CellTypist trains from the reference unless `--celltypist_model` is supplied.
- scANVI loads `--scanvi_model_dir` when supplied; otherwise it trains a model, optionally initialized from `--scvi_reference_model_dir`.
- PanHumanPy always annotates the query directly.

Required when using pretrained or external artifacts:

- `--celltypist_model`: required to reuse a pretrained CellTypist model.
- `--scanvi_model_dir`: required to reuse a pretrained scANVI model directory.
- `--scvi_reference_model_dir`: required only when initializing scANVI training from a pretrained scVI reference model. Do not combine this with `--scanvi_model_dir`.

CellTypist optional controls:

- `--celltypist_training_mode`: `standard` or `detailed`. Default: `detailed`.
- `--celltypist_balance_cells_per_label`: maximum reference cells per label in detailed mode. Default: `500`.
- `--celltypist_feature_selection`: run tutorial-style feature selection in detailed mode. Default: `true`.
- `--celltypist_feature_selection_top_genes`: top genes per label for feature selection. Default: `100`.
- `--celltypist_feature_selection_max_iter`: quick SGD feature-selection iterations. Default: `5`.
- `--celltypist_max_iter`: final CellTypist fit iterations. Default: `100`.

CellTypist model reuse checks that at least 80% of model genes are present in the query. This catches gene identifier mismatches before silently filling most model features with zeros.

scANVI optional controls:

- `--batch_key`: optional column present in both reference and query `obs`; omitted values stay as real `None`/`null`. Use this when batches should be modeled.
- `--scanvi_categorical_covariate_keys`: optional comma-separated additional covariate columns present in both objects.
- `--scanvi_hvg_n_top_genes`: HVGs for new scVI/scANVI training. Set `0` to keep all shared genes. Default: `2500`.
- `--scanvi_query_max_epochs`: max epochs for query mapping from a pretrained scVI reference model. Default: `100`.
- `--unknown_label`: unlabeled query category for scANVI training. Default: `Unknown`.

PanHumanPy optional controls:

- `--panhumanpy_feature_names_col`: optional query `var` column containing gene symbols when `query.var_names` are not gene symbols.

Output naming and merge behavior:

- `--outdir`: output directory. Default: `autoannotate-results-<project_tag>`.
- `--scanvi_obs_col`: scANVI label column name in output. Default: `scanvi_label`.
- `--scanvi_score_col`: scANVI confidence column name in output. Default: `scanvi_confidence`.
- `--query_cluster_col`: existing query `obs` cluster column used for cluster summaries. Default: `leiden`.
- `--compute_missing_clusters`: compute Leiden clusters in the merge step if `--query_cluster_col` is absent. Default: `false`.
- `--marker_genes_n_top`: number of marker genes to report per cluster. Default: `50`.
- `--marker_genes_method`: differential-expression test for marker ranking, either `wilcoxon` or `t-test`. Default: `wilcoxon`.

When missing clusters are explicitly computed, the merge step reuses existing neighbors or `X_pca` if present. Otherwise it runs PCA/neighbors/Leiden on a normalized, log-transformed copy so the final query matrix is not modified.

Parameter defaults live in [nextflow.config](nextflow.config).

## Validation

Before model steps run, the pipeline checks that:

- input files exist
- `ref_label_col` exists in `ref_h5ad.obs`
- `batch_key`, if supplied, exists in both reference and query `obs`
- all `scanvi_categorical_covariate_keys`, if supplied, exist in both objects
- `panhumanpy_feature_names_col`, if supplied, exists in `query_h5ad.var`
- reference timepoint filters are valid and keep at least one reference cell
- reference and query cell indices are unique
- reference and query gene indices are unique
- reference and query gene identifiers use the same broad scheme, such as symbols or Ensembl IDs
- reference and query share at least one gene

The merge step also validates that every annotation CSV has unique `cell_id` values matching `query.obs_names` exactly. Missing or extra cells now fail the merge instead of becoming silent `NaN` values.

## Outputs

Published files:

- `<project_tag>.annotated.h5ad`
- `<project_tag>.cluster_annotation_summary.csv`
- `<project_tag>.cluster_marker_genes.tsv`

The annotated `.h5ad` is the query object with additional `obs` columns.

Method columns copied into `query.obs`:

- `celltypist_predicted_label`: raw per-cell CellTypist label before majority voting.
- `celltypist_majority_voting`: CellTypist label after majority voting.
- `celltypist_confidence`: probability for the majority-voting label when available; otherwise the row maximum.
- `<scanvi_obs_col>`: scANVI predicted label. Default: `scanvi_label`.
- `<scanvi_score_col>`: maximum scANVI class probability. Default: `scanvi_confidence`.
- `panhumanpy_full_hierarchical_label`: PanHumanPy full hierarchy label.
- `panhumanpy_level_zero_label`: PanHumanPy broadest label.
- `panhumanpy_final_level_label`: PanHumanPy final selected label.
- `panhumanpy_confidence`: PanHumanPy final-level softmax probability.
- `panhumanpy_azimuth_broad`, `panhumanpy_azimuth_medium`, `panhumanpy_azimuth_fine`: refined PanHumanPy labels when returned.

Comparison columns created in `query.obs`:

- `celltypist_scanvi_agree`: per-cell boolean; `true` only when both methods have non-missing labels and `celltypist_majority_voting` equals the scANVI label.
- `celltypist_scanvi_cluster_agreement_fraction`: cluster-level fraction copied onto each cell in the cluster; the fraction of cells where `celltypist_scanvi_agree` is `true`.

The cluster summary is grouped by `--query_cluster_col` and includes:

- `<query_cluster_col>`: cluster identifier.
- `n_cells`: number of cells in the cluster.
- `celltypist_majority_top_label`: most frequent CellTypist majority-voting label.
- `celltypist_majority_top_fraction`: fraction of non-missing CellTypist labels assigned to that top label.
- `scanvi_top_label`: most frequent scANVI label.
- `scanvi_top_fraction`: fraction of non-missing scANVI labels assigned to that top label.
- `panhumanpy_top_label`: most frequent PanHumanPy final-level label, reported independently.
- `panhumanpy_top_fraction`: fraction of non-missing PanHumanPy labels assigned to that top label.
- `celltypist_scanvi_cluster_modal_match`: whether the CellTypist and scANVI modal labels match in the cluster.
- `celltypist_scanvi_cluster_agreement_fraction`: fraction of cells in the cluster where CellTypist and scANVI agree.

`celltypist_scanvi_agree` and `celltypist_scanvi_cluster_modal_match` are intentionally different summaries. A cluster can have matching modal labels while many individual cells disagree. The modal-match boolean is therefore kept in the cluster summary, while the more interpretable agreement fraction is copied to `query.obs`.

The marker-gene table is grouped by `--query_cluster_col` and written to `<project_tag>.cluster_marker_genes.tsv`. Marker genes are ranked with Scanpy `rank_genes_groups` on a normalized, log-transformed copy of the query object, so marker calculation does not modify the matrix written to `<project_tag>.annotated.h5ad`. The table includes:

- `cluster`: cluster identifier.
- `rank`: one-based marker rank within the cluster.
- `gene`: gene identifier from `query.var_names`.
- `score`: Scanpy rank score for the selected method.
- `logfoldchange`: estimated log fold change for the cluster against the remaining cells.
- `p_value`: nominal p-value.
- `p_value_adj`: adjusted p-value.

If the query contains fewer than two clusters, the marker-gene TSV is still written with headers but no marker rows because there is no between-cluster comparison.

Intermediate files remain in Nextflow work directories, including:

- `validation_manifest.json`
- `shared_genes.txt`
- `celltypist_model.pkl`
- `celltypist_model_metadata.json`
- `celltypist_predictions.csv`
- `scanvi_model/`
- `scanvi_model_metadata.json`
- `scanvi_predictions.csv`
- `panhumanpy_predictions.csv`

Nextflow timeline, report, and trace files are written under `autoannotate-reports-<project_tag>/`.

## Demo Data

To create local reference and query inputs, run:

```bash
python scripts/prepare_scanpy_dataset.py --outdir data/demo
```

By default this downloads the full PBMC68K object into `data/demo/raw/`, writes an 80/20 reference/query split, and records split metadata:

- `data/demo/reference.h5ad`
- `data/demo/query.h5ad`
- `data/demo/dataset_split_metadata.json`

Use the generated files for a training-mode demo run:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad data/demo/reference.h5ad \
  --query_h5ad data/demo/query.h5ad \
  --project_tag demo \
  --ref_label_col bulk_labels \
  --query_cluster_col louvain \
  --outdir results/demo
```

The helper also supports other sources:

```bash
python scripts/prepare_scanpy_dataset.py --dataset pbmc68k_reduced --outdir data/pbmc68k_reduced_demo
python scripts/prepare_scanpy_dataset.py --dataset pbmc3k --outdir data/pbmc3k_demo
python scripts/prepare_scanpy_dataset.py --input-h5ad /path/to/input.h5ad --outdir data/custom_demo
```

`--dataset pbmc68k` is the default and downloads the larger Figshare-hosted PBMC68K object into the selected output directory before splitting it. `--dataset pbmc68k_reduced` and `--dataset pbmc3k` use datasets provided directly by Scanpy.

## Examples

Train CellTypist and scANVI from the reference:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag test_run \
  --ref_label_col cell_type \
  --batch_key donor \
  --query_cluster_col leiden
```

Reuse a pretrained CellTypist model while training scANVI:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag test_run \
  --ref_label_col cell_type \
  --batch_key donor \
  --celltypist_model /path/to/celltypist_model.pkl
```

Reuse pretrained CellTypist and scANVI artifacts:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag test_run \
  --ref_label_col cell_type \
  --celltypist_model /path/to/celltypist_model.pkl \
  --scanvi_model_dir /path/to/scanvi_model_dir
```

Initialize scANVI training from a pretrained scVI reference model:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag test_run \
  --ref_label_col cell_type \
  --batch_key donor \
  --scvi_reference_model_dir /path/to/scvi_model_dir
```

Allow the merge step to compute missing Leiden clusters:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag test_run \
  --ref_label_col cell_type \
  --compute_missing_clusters true
```

## Assumptions

- Input `.h5ad` files should already be normalized/preprocessed appropriately for the chosen annotation methods.
- Gene alignment for CellTypist and scANVI uses the reference-query intersection from validation.
- PanHumanPy receives the original query feature space, optionally using `--panhumanpy_feature_names_col` for symbols.
- CellTypist, scANVI, and PanHumanPy outputs are kept separate.
- PanHumanPy is not compared with CellTypist or scANVI because it uses a different reference.
