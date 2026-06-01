# nf-autoannotate

`nf-autoannotate` is a Nextflow pipeline for annotating a query `.h5ad` dataset with three independent methods:

- CellTypist
- scANVI
- PanHumanPy / Pan-Human Azimuth

The pipeline writes the method outputs back into the original query AnnData object and produces a cluster-level summary table. It does not create a consensus label.

## Workflow

The workflow is implemented in [main.nf](main.nf) and calls helper scripts in [scripts/](scripts):

- [validate_inputs.py](scripts/validate_inputs.py)
- [train_celltypist.py](scripts/train_celltypist.py)
- [annotate_celltypist.py](scripts/annotate_celltypist.py)
- [train_scanvi.py](scripts/train_scanvi.py)
- [annotate_scanvi.py](scripts/annotate_scanvi.py)
- [annotate_panhumanpy.py](scripts/annotate_panhumanpy.py)
- [merge_annotations.py](scripts/merge_annotations.py)

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

## Model And Annotation Options

CellTypist:

- `--celltypist_mode`: `train` or `pretrained`. If unset, the pipeline uses `pretrained` when `--celltypist_model` is supplied, otherwise `train`.
- `--celltypist_model`: pretrained CellTypist model file.
- `--celltypist_training_mode`: `standard` or `detailed`. Default: `detailed`.
- `--celltypist_balance_cells_per_label`: maximum reference cells per label in detailed mode. Default: `500`.
- `--celltypist_feature_selection`: run tutorial-style feature selection in detailed mode. Default: `true`.
- `--celltypist_feature_selection_top_genes`: top genes per label for feature selection. Default: `100`.
- `--celltypist_feature_selection_max_iter`: quick SGD feature-selection iterations. Default: `5`.
- `--celltypist_max_iter`: final CellTypist fit iterations. Default: `100`.

CellTypist model reuse now checks that at least 80% of model genes are present in the query. This catches gene identifier mismatches before silently filling most model features with zeros.

scANVI:

- `--scanvi_mode`: `train` or `pretrained`. If unset, the pipeline uses `pretrained` when `--scanvi_model_dir` is supplied, otherwise `train`.
- `--scanvi_model_dir`: pretrained scANVI model directory.
- `--scvi_reference_model_dir`: pretrained scVI reference model directory used to initialize scANVI training.
- `--batch_key`: optional column present in both reference and query `obs`; omitted values stay as real `None`/`null`.
- `--scanvi_categorical_covariate_keys`: optional comma-separated additional covariate columns present in both objects.
- `--scanvi_hvg_n_top_genes`: HVGs for new scVI/scANVI training. Set `0` to keep all shared genes. Default: `2500`.
- `--scanvi_query_max_epochs`: max epochs for query mapping from a pretrained scVI reference model. Default: `100`.
- `--unknown_label`: unlabeled query category for scANVI training. Default: `Unknown`.

PanHumanPy:

- `--panhumanpy_feature_names_col`: optional query `var` column containing gene symbols when `query.var_names` are not gene symbols.

Output naming and merge behavior:

- `--outdir`: output directory. Default: `autoannotate-results-<project_tag>`.
- `--scanvi_obs_col`: scANVI label column name in output. Default: `scanvi_label`.
- `--scanvi_score_col`: scANVI confidence column name in output. Default: `scanvi_confidence`.
- `--query_cluster_col`: existing query `obs` cluster column used for cluster summaries. Default: `leiden`.
- `--compute_missing_clusters`: compute Leiden clusters in the merge step if `--query_cluster_col` is absent. Default: `false`.

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

## Examples

Train both reference-driven models:

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

Reuse a pretrained CellTypist model:

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

Reuse CellTypist and scANVI models:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag test_run \
  --ref_label_col cell_type \
  --celltypist_mode pretrained \
  --celltypist_model /path/to/celltypist_model.pkl \
  --scanvi_mode pretrained \
  --scanvi_model_dir /path/to/scanvi_model_dir
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

## Resources And Containers

The configured process labels in [nextflow.config](nextflow.config) are:

- `validation`: lightweight input checks.
- `normal`: single-core CPU tasks.
- `cpu4`: four-core CellTypist training.
- `gpu`: scANVI training and annotation on `gpu-normal`.

`SCANVI_ANNOTATE` uses the shared `gpu` label but overrides memory to 128 GB. The unused `long` label and duplicate GPU labels were removed.

Two fixed Singularity images are configured:

- `/nfs/cellgeni/singularity/images/autoannotate-0.1.0.sif` for CellTypist, scVI/scANVI, validation, and merge steps.
- `/nfs/cellgeni/singularity/images/autoannotate-0.1.0-phpy.sif` for `PANHUMANPY_ANNOTATE` only.

PanHumanPy uses a separate image because its TensorFlow dependency stack is isolated from the main scverse/RAPIDS image.

## Assumptions

- Input `.h5ad` files should already be normalized/preprocessed appropriately for the chosen annotation methods.
- Gene alignment for CellTypist and scANVI uses the reference-query intersection from validation.
- PanHumanPy receives the original query feature space, optionally using `--panhumanpy_feature_names_col` for symbols.
- CellTypist, scANVI, and PanHumanPy outputs are kept separate.
- PanHumanPy is not compared with CellTypist or scANVI because it uses a different reference.

## Development

Core files:

- [main.nf](main.nf)
- [nextflow.config](nextflow.config)
- [Dockerfile](Dockerfile)
- [Dockerfile.panhumanpy](Dockerfile.panhumanpy)
- [scripts/](scripts)
- [tests/](tests)
