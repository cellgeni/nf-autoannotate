# nf-autoannotate

`nf-autoannotate` is a Nextflow pipeline that annotates a query `.h5ad` dataset with three methods in parallel:

- CellTypist
- scANVI
- PanHumanPy / Pan-Human Azimuth

It then merges the per-cell predictions back into the original query AnnData object and writes a final annotated `.h5ad`.

## What The Pipeline Does

The pipeline takes:

- a labeled reference `.h5ad`
- an unlabeled or partially labeled query `.h5ad`

It performs these stages:

1. Validate the inputs and compute a shared feature set.
2. Prepare a CellTypist model.
   This can either train from the reference dataset or reuse a pretrained model.
3. Prepare an scANVI model.
   This can either train from the reference/query pair or reuse a pretrained model directory.
4. Run CellTypist annotation on the query.
5. Run scANVI annotation on the query.
6. Run PanHumanPy annotation on the query.
7. Merge all outputs into the query `.h5ad`.

The final output contains all methods independently. This version does not compute a consensus label.

CellTypist training supports two modes:

- `standard`
  Train directly on the validated reference dataset restricted to shared genes.
- `detailed`
  Follow the official large-scale cross-dataset tutorial approach: cap the number of cells per label, optionally run a quick SGD feature-selection pass, and train the final model on the balanced subset and selected genes. Source: [CellTypist cross-dataset tutorial](https://celltypist.readthedocs.io/en/latest/notebook/celltypist_tutorial_cv.html).

## Pipeline Structure

The workflow is implemented in [`main.nf`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/main.nf) and calls helper scripts in [`scripts/`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/scripts):

- `validate_inputs.py`
- `train_celltypist.py`
- `annotate_celltypist.py`
- `train_scanvi.py`
- `annotate_scanvi.py`
- `annotate_panhumanpy.py`
- `merge_annotations.py`

## Inputs

### Required Parameters

- `--ref_h5ad`
  Path to the reference AnnData file.
- `--query_h5ad`
  Path to the query AnnData file that will receive annotations.
- `--project_tag`
  Used in report paths and in the final output filename.
- `--outdir`
  Directory where final published outputs are copied.
  Default: `autoannotate-results-<project_tag>`
- `--ref_label_col`
  Column in `ref_h5ad.obs` containing reference cell type labels.
- `--ref_timepoint_col`
  Optional reference `obs` column containing timepoint information used to restrict the reference before training.
  Must be used together with `--ref_timepoints`.
- `--ref_timepoints`
  Optional comma-separated list of reference timepoint values to keep for training, for example `E12,E14`.
  Must be used together with `--ref_timepoint_col`.

### Optional Parameters

- `--batch_key`
  Optional column present in both `ref_h5ad.obs` and `query_h5ad.obs` used for scVI/scANVI setup.
  If omitted, scVI/scANVI is set up with `batch_key=None`.
- `--scanvi_categorical_covariate_keys`
  Optional comma-separated `obs` columns present in both `ref_h5ad` and `query_h5ad` to pass to `SCVI.setup_anndata(..., categorical_covariate_keys=...)`.
  These should be additional covariates, not the batch key itself.
  Default: unset
- `--scanvi_hvg_n_top_genes`
  Number of highly variable genes to use for newly trained scVI/scANVI models.
  Set to `0` to keep all shared genes.
  Default: `2500`
- `--panhumanpy_feature_names_col`
  Optional query `var` column containing gene symbols for PanHumanPy when `query.var_names` are not gene symbols.
  This is passed to `panhumanpy.AzimuthNN(..., feature_names_col=...)`.
  Default: unset
- `--panhumanpy_container`
  Singularity image used only for the `PANHUMANPY_ANNOTATE` process.
  Default: `/nfs/cellgeni/singularity/images/autoannotate-panhumanpy-0.1.0.sif`
- `--celltypist_model`
  Path to a pretrained CellTypist model file.
- `--scanvi_model_dir`
  Path to a pretrained scANVI model directory.
- `--scvi_reference_model_dir`
  Path to a pretrained scVI reference model directory to initialize scANVI training.
  This is distinct from `--scanvi_model_dir`, which skips scANVI training entirely.
- `--celltypist_mode`
  One of `train` or `pretrained`.
  If unset, the pipeline auto-selects `pretrained` when `--celltypist_model` is supplied, otherwise `train`.
- `--scanvi_mode`
  One of `train` or `pretrained`.
  If unset, the pipeline auto-selects `pretrained` when `--scanvi_model_dir` is supplied, otherwise `train`.
- `--scanvi_obs_col`
  Output column name for scANVI predicted labels.
  Default: `scanvi_label`
- `--scanvi_score_col`
  Output column name for scANVI confidence.
  Default: `scanvi_confidence`
- `--query_cluster_col`
  Query `obs` column used for CellTypist/scANVI cluster-level comparison and independent per-method cluster summaries during merge.
  If the column is absent, the merge step computes Leiden clusters with this
  column name.
  Default: `leiden`
- `--scanvi_query_max_epochs`
  Max epochs for scANVI query mapping when using a pretrained scVI reference model.
  Default: `100`
- `--unknown_label`
  Label used for unlabeled query cells when training scANVI.
  Default: `Unknown`
- `--celltypist_balance_cells_per_label`
  Maximum number of reference cells to use per label when training CellTypist.
  Default: `500`
- `--celltypist_training_mode`
  One of `standard` or `detailed`.
  Default: `detailed`
- `--celltypist_feature_selection`
  Enable tutorial-style CellTypist feature selection before final training.
  Default: `true`
- `--celltypist_feature_selection_top_genes`
  Number of top genes per label to retain during CellTypist feature selection.
  Default: `100`
- `--celltypist_feature_selection_max_iter`
  Iteration count for the quick SGD pass used for CellTypist feature selection.
  Default: `5`
- `--celltypist_max_iter`
  Iteration count for the final CellTypist logistic-regression fit.
  Default: `100`

The parameter defaults are defined in [`nextflow.config`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/nextflow.config).

## Validation Rules

Before any model runs, the pipeline validates that:

- `ref_h5ad` exists
- `query_h5ad` exists
- `ref_label_col` exists in `ref_h5ad.obs`
- if supplied, `batch_key` exists in both `ref_h5ad.obs` and `query_h5ad.obs`
- if supplied, every `scanvi_categorical_covariate_keys` column exists in both `ref_h5ad.obs` and `query_h5ad.obs`
- if supplied, `panhumanpy_feature_names_col` exists in `query_h5ad.var`
- if timepoint filtering is requested, `ref_timepoint_col` exists in `ref_h5ad.obs`
- if timepoint filtering is requested, all requested `ref_timepoints` exist in the reference
- cell indices are unique in both objects
- the reference and query share at least one gene

If validation passes, the pipeline writes a `validation_manifest.json` file. This manifest contains:

- the resolved input file paths
- `ref_label_col`
- `batch_key`
- `categorical_covariate_keys`
- `panhumanpy_feature_names_col`
- the optional reference timepoint filter column and selected timepoints
- the shared gene count
- reference and query cell counts

The full validated feature intersection is written separately as `shared_genes.txt` and passed to downstream steps that need a consistent gene order.

## Execution Modes

### Train Both Models

If neither pretrained artifact is supplied, the pipeline:

- trains CellTypist from `ref_h5ad`
- trains scANVI from `ref_h5ad` and `query_h5ad`
- annotates the query with CellTypist, scANVI, and PanHumanPy

### Reuse Pretrained CellTypist

If `--celltypist_model` is supplied, or `--celltypist_mode pretrained` is set, the pipeline skips CellTypist training and uses that model for inference.

If CellTypist is trained in `standard` mode, the training logic:

- optionally filters the reference to the requested timepoint subset
- subsets the reference to the validated shared genes
- trains the final model directly on that dataset

If CellTypist is trained in `detailed` mode, the training logic:

- optionally filters the reference to the requested timepoint subset
- subsets the reference to the validated shared genes
- downsamples each label to at most `--celltypist_balance_cells_per_label` cells
- optionally trains a quick SGD model to rank informative genes
- unions the top genes across labels
- fits the final CellTypist model on the balanced subset using the selected genes

When CellTypist annotation runs with majority voting enabled, the pipeline stores:

- the raw CellTypist `predicted_labels`
- the CellTypist `majority_voting` labels
- one CellTypist confidence column derived from the probability assigned to the final majority-voting label when available, with fallback to the row maximum otherwise

These CellTypist output columns are fixed:

- `celltypist_predicted_label`
- `celltypist_majority_voting`
- `celltypist_confidence`

### PanHumanPy / Pan-Human Azimuth Annotation

PanHumanPy runs directly on the query AnnData object using the high-level `panhumanpy.AzimuthNN` interface. The pipeline also calls the built-in refinement step when available.

PanHumanPy output columns are fixed:

- `panhumanpy_full_hierarchical_label`
- `panhumanpy_level_zero_label`
- `panhumanpy_final_level_label`
- `panhumanpy_confidence`
- `panhumanpy_azimuth_broad`
- `panhumanpy_azimuth_medium`
- `panhumanpy_azimuth_fine`

### Reuse Pretrained scANVI

If `--scanvi_model_dir` is supplied, or `--scanvi_mode pretrained` is set, the pipeline skips scANVI training and uses that model directory for inference.

When the pipeline trains scVI/scANVI itself, scVI is fit on highly variable genes only. The default is `--scanvi_hvg_n_top_genes 2500`; set it to `0` to use all validated shared genes. If `--scvi_reference_model_dir` or `--scanvi_model_dir` is supplied, the model's existing feature space is used instead.

### Initialize scANVI From A Pretrained scVI Reference Model

If `--scvi_reference_model_dir` is supplied while `scanvi_mode` remains `train`, the pipeline:

- optionally filters the reference to the requested timepoint subset before scANVI training
- loads the pretrained scVI reference model
- initializes a reference scANVI model from it
- trains scANVI on the labeled reference only
- maps the query during the annotation step using scvi-tools query loading

### Mixed Mode

The methods are independent. CellTypist and scANVI can be trained or reused separately, while PanHumanPy runs directly on the query.

## Outputs

### Final Output

The merge step publishes final outputs with `publishDir` into:

- `--outdir` if provided
- otherwise `autoannotate-results-<project_tag>`

Published files:

- `<project_tag>.annotated.h5ad`
- `<project_tag>.cluster_annotation_summary.csv`

This file is the original query AnnData with additional columns in `obs`. The merge step copies method outputs into `query.obs` and creates only CellTypist-vs-scANVI comparison columns. PanHumanPy is not compared with the other methods because it uses a different reference.

Per-cell annotation columns copied into `query.obs`:

- `celltypist_predicted_label`: raw per-cell CellTypist prediction before majority voting.
- `celltypist_majority_voting`: CellTypist label after neighborhood majority voting; this is the CellTypist label used for CellTypist-vs-scANVI comparisons.
- `celltypist_confidence`: CellTypist probability assigned to the majority-voting label when available, otherwise the maximum probability in that cell.
- `<scanvi_obs_col>`, default `scanvi_label`: scANVI predicted label.
- `<scanvi_score_col>`, default `scanvi_confidence`: maximum scANVI class probability for the predicted label.
- `panhumanpy_full_hierarchical_label`: PanHumanPy full hierarchy label from the Pan-Human Azimuth reference.
- `panhumanpy_level_zero_label`: PanHumanPy broadest level-zero label.
- `panhumanpy_final_level_label`: PanHumanPy final selected label for the cell.
- `panhumanpy_confidence`: PanHumanPy final-level softmax probability.
- `panhumanpy_azimuth_broad`, `panhumanpy_azimuth_medium`, `panhumanpy_azimuth_fine`: refined PanHumanPy labels when the refinement step returns them.

Columns created by `merge_annotations.py` in `query.obs`:

- `celltypist_scanvi_agree`: per-cell boolean; `true` when `celltypist_majority_voting` equals the scANVI label column.
- `celltypist_scanvi_cluster_modal_match`: cluster-level boolean copied onto every cell in the cluster; `true` when the top CellTypist majority-voting label in that cluster equals the top scANVI label in that cluster.
- `celltypist_scanvi_cluster_agreement_fraction`: cluster-level fraction copied onto every cell in the cluster; the fraction of cells in that cluster where `celltypist_scanvi_agree` is `true`.

The cluster summary is grouped by `--query_cluster_col` and includes:

- `<query_cluster_col>`: cluster identifier.
- `n_cells`: number of cells in the cluster.
- `celltypist_majority_top_label`: most frequent CellTypist majority-voting label in the cluster.
- `celltypist_majority_top_fraction`: fraction of non-missing CellTypist labels in the cluster assigned to `celltypist_majority_top_label`.
- `scanvi_top_label`: most frequent scANVI label in the cluster.
- `scanvi_top_fraction`: fraction of non-missing scANVI labels in the cluster assigned to `scanvi_top_label`.
- `panhumanpy_top_label`: most frequent PanHumanPy final-level label in the cluster, reported independently.
- `panhumanpy_top_fraction`: fraction of non-missing PanHumanPy final-level labels in the cluster assigned to `panhumanpy_top_label`.
- `celltypist_scanvi_cluster_modal_match`: whether the cluster-level CellTypist and scANVI modal labels match.
- `celltypist_scanvi_cluster_agreement_fraction`: fraction of cells in the cluster where CellTypist and scANVI agree.

### Intermediate Outputs

The workflow also produces intermediates in task work directories, including:

- `validation_manifest.json`
- `shared_genes.txt`
- `celltypist_model.pkl`
- `celltypist_model_metadata.json`
- `celltypist_predictions.csv`
- `scanvi_model/`
- `scanvi_model_metadata.json`
- `scanvi_predictions.csv`
- `panhumanpy_predictions.csv`

### Reports

Nextflow reporting is enabled and writes:

- timeline report
- execution report
- trace file

under `autoannotate-reports-<project_tag>/`.

## Example Commands

### Train Both Models

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag test_run \
  --ref_label_col cell_type \
  --batch_key donor
```

### Reuse A Pretrained CellTypist Model

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

### Reuse Both Models

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag test_run \
  --ref_label_col cell_type \
  --batch_key donor \
  --celltypist_mode pretrained \
  --celltypist_model /path/to/celltypist_model.pkl \
  --scanvi_mode pretrained \
  --scanvi_model_dir /path/to/scanvi_model_dir
```

## Resource Configuration

Default process labels are configured in [`nextflow.config`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/nextflow.config):

- `normal`
  general lightweight tasks
- `normal4core`
  CellTypist training and annotation
- `normalgpu`
  scANVI training and annotation

The pipeline uses two containers:

- `/nfs/cellgeni/singularity/images/autoannotate-0.1.0.sif` for the main CellTypist/scVI/scANVI steps
- `--panhumanpy_container`, defaulting to `/nfs/cellgeni/singularity/images/autoannotate-panhumanpy-0.1.0.sif`, for `PANHUMANPY_ANNOTATE` only

This keeps PanHumanPy's TensorFlow dependency stack isolated from the main scverse/RAPIDS image. The pipeline runs with LSF executor settings defined in `nextflow.config`.

## Assumptions And Current Behavior

- The input `.h5ad` files are assumed to already be preprocessed appropriately for annotation.
- The pipeline only performs structural validation and feature alignment.
- Gene alignment is based on the reference-query intersection computed during validation.
- PanHumanPy receives the original query feature space; the upstream PanHumanPy guidance is to keep all quantified genes for annotation.
- CellTypist `detailed` training is label-balanced by capped per-label downsampling, which is intended to reduce domination by very common labels.
- CellTypist, scANVI, and PanHumanPy outputs are kept separate.
- Only the top predicted label and one confidence value per method are stored.
- The final merge updates `query.obs` and preserves the rest of the AnnData object.

## Development Notes

Core files:

- [`main.nf`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/main.nf)
- [`nextflow.config`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/nextflow.config)
- [`Dockerfile`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/Dockerfile)
- [`Dockerfile.panhumanpy`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/Dockerfile.panhumanpy)
- [`scripts/validate_inputs.py`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/scripts/validate_inputs.py)
- [`scripts/train_celltypist.py`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/scripts/train_celltypist.py)
- [`scripts/annotate_celltypist.py`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/scripts/annotate_celltypist.py)
- [`scripts/train_scanvi.py`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/scripts/train_scanvi.py)
- [`scripts/annotate_scanvi.py`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/scripts/annotate_scanvi.py)
- [`scripts/annotate_panhumanpy.py`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/scripts/annotate_panhumanpy.py)
- [`scripts/merge_annotations.py`](/lustre/scratch127/cellgen/cellgeni/bc8/_repos/nf-autoannotate/scripts/merge_annotations.py)
