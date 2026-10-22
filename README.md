# nf-autoannotate

`nf-autoannotate` is a Nextflow pipeline for annotating a query `.h5ad` dataset with three independent methods:

- CellTypist
- scANVI
- PanHumanPy / Pan-Human Azimuth

The pipeline writes the method outputs back into the original query AnnData object and produces a cluster-level summary table. It does not create a consensus label.

All three methods run by default. Use `--tools` to run a subset, for example `--tools celltypist,panhumanpy`.

## Workflow

The workflow is implemented in [main.nf](main.nf), with supporting Python scripts under [scripts/](scripts).

[prepare_scanpy_dataset.py](scripts/prepare_scanpy_dataset.py) is a utility script for creating local demo reference/query `.h5ad` files from a Scanpy dataset or an existing AnnData file.

The default workflow stages are:

Once per run, from the reference alone:

1. Validate the reference AnnData file, then write `reference_manifest.json` and `reference_genes.txt`.
2. Train or load a CellTypist model.
3. Train or load scVI/scANVI reference models.

Once per query:

4. Validate the query AnnData file and record its coverage of the reference gene space.
5. Annotate the query with CellTypist.
6. Annotate the query with scANVI query mapping.
7. Annotate the query with PanHumanPy.
8. Merge method outputs into the query `.h5ad` and write a cluster annotation summary.

Stages 2, 3, and 5-7 only run for the tools named in `--tools`. Stage 4 always runs, so a subset run still validates its inputs. Stage 8 runs unless it is turned off with `--merge_annotations false`; see [Skipping the Merge Step](#skipping-the-merge-step).

Models are trained from the reference and never from a reference-query intersection, which is what lets one training run serve any number of queries. Their final feature space is narrower than the reference, because CellTypist feature selection and scVI HVG selection each reduce it. Every query is aligned to the trained model's own feature list at annotation time, with genes it lacks zero-filled.

Two other entry workflows are available:

- `-entry TRAIN_MODELS` runs reference validation and model training, publishes CellTypist plus scVI/scANVI artifacts, and skips query annotation, PanHumanPy, and merge outputs. It requires `--tools` to include at least one of `celltypist` or `scanvi`, and because training never reads a query it does not need `--query_h5ad` at all.
- `-entry ANNOTATE_WITH_MODELS` annotates queries from models that already exist, with the (model, query) pairs listed in a CSV or TSV. See [Annotating From Existing Models](#annotating-from-existing-models).

## Tool Selection

- `--tools`: annotation tools to run. Defaults to `celltypist,scanvi,panhumanpy`.

Accepts any subset of `celltypist`, `scanvi`, and `panhumanpy`, as a comma-separated string on the command line or as either a string or a list in a `-params-file` YAML:

```yaml
tools: celltypist,panhumanpy
# or
tools:
  - celltypist
  - panhumanpy
```

Selecting a subset changes what the run requires and produces:

- A reference is only needed to train CellTypist or scANVI. `--tools panhumanpy` annotates the query on its own, and so does a run that supplies pretrained models for every selected reference-based tool; in both cases `--ref_h5ad` and `--ref_label_col` are optional and the reference validation stage is skipped.
- The cluster annotation summary only contains columns for the tools that ran. The `celltypist_scanvi_*` agreement columns need both CellTypist and scANVI, so they are omitted otherwise.
- Supplying an option for a tool that was not selected is an error rather than a silent no-op. For example `--tools panhumanpy --celltypist_model model.pkl` fails immediately.

## Required Inputs

- `--query_h5ad`: query AnnData file(s) to annotate. A path, a glob, or a list. See [Multiple Queries](#multiple-queries).
- `--project_tag`: output filename prefix.
- `--ref_h5ad`: reference AnnData file. Required when a reference-based tool is trained, that is when `--tools` includes `celltypist` or `scanvi` without a matching pretrained model. Optional otherwise, and still used for the reference-coverage diagnostic when supplied.
- `--ref_label_col`: reference `obs` column containing cell type labels. Required whenever a reference h5ad is used.

## Multiple Queries

`--query_h5ad` accepts one path, a shell glob, or a list, and each query is annotated independently:

```yaml
# any of these
query_h5ad: /data/donor1.h5ad
query_h5ad: /data/queries/*.h5ad
query_h5ad:
  - /data/donor1.h5ad
  - /data/donor2.h5ad
```

```bash
nextflow run main.nf --query_h5ad '/data/queries/*.h5ad' ...
```

Quote globs on the command line so the shell passes the pattern through rather than expanding it into several arguments.

The sample ID for each query is its filename without the final extension, so `/data/donor1.h5ad` becomes `donor1`. Sample IDs must be unique; two files that reduce to the same ID fail before any task is submitted. Sample IDs name the per-query output subdirectory and appear in every per-query filename.

Reference validation and model training run **once** regardless of how many queries are supplied, because no query takes part in training. Each query then validates, annotates, and merges independently, so one slow query does not hold up the others.

CellTypist needs log1p expression normalised to 10,000 counts per cell. The CellTypist step uses `adata.X` as it is when it already has that form. Otherwise it normalises raw counts taken from `adata.X`, then `adata.layers['counts']`, then `adata.raw.X`, whichever holds them first, and fails if none does. For a `--tools celltypist` run, validation applies the same rule up front, so a query with normalised `adata.X` is accepted. Any run that includes `scanvi` or `panhumanpy` still requires raw counts in `adata.X`.

Queries are never combined. There is no cross-query summary, and no query influences another's labels.

Models are trained from the reference rather than from a reference-query intersection, but their final feature space is narrower than the reference: CellTypist keeps a per-label subset after feature selection, and scVI keeps `--scanvi_hvg_n_top_genes` highly variable genes. A query that lacks some of those features has them zero-filled when it is aligned to the model.

Two separate numbers are reported as a result:

- **Reference coverage**, recorded per query in `<sample_id>.query_manifest.json` and warned on below 80%. This is a broad diagnostic against the whole reference gene space, not a statement about model compatibility.
- **Model gene overlap**, printed by each annotation step against that model's actual feature list. This is the authoritative check: CellTypist fails outright below 80%, and scANVI warns.

A query can show good reference coverage and still overlap poorly with a model's selected features, so treat the per-model number as the one that matters.

Optional reference filtering:

- `--ref_filter_col`: reference `obs` column used to filter the reference before model training.
- `--ref_filter_values`: comma-separated values to keep from `--ref_filter_col`, for example `E12,E14` or `treated,control`.

`--ref_filter_col` and `--ref_filter_values` must be supplied together, and only affect model training, so they require a reference-based tool that is actually being trained.

Run parameters are normally defined in a YAML file passed with `-params-file`. Start from [examples/autoannotate.params.yml](examples/autoannotate.params.yml), edit the values for your run, and launch the workflow with that file. Parameter defaults and runtime settings live in [nextflow.config](nextflow.config), and user-facing parameters are described and validated by [nextflow_schema.json](nextflow_schema.json).

## Annotating From Existing Models

The default workflow applies one set of models to every query. When different queries need
different models, use `-entry ANNOTATE_WITH_MODELS` and list the pairs in a sheet:

```bash
nextflow run main.nf -entry ANNOTATE_WITH_MODELS \
  --model_sheet /path/to/model_sheet.csv \
  --project_tag reuse_run \
  --outdir results/reuse_run
```

```csv
tool,model,model_id,query_h5ad
celltypist,/runs/gut/celltypist_model.pkl,gut_ct,/data/obj1.h5ad
celltypist,/runs/gut/celltypist_model.pkl,gut_ct,/data/obj2.h5ad
scanvi,/runs/lung/scanvi_model,lung_scanvi,/data/obj3.h5ad
```

That sheet runs the gut CellTypist model over `obj1` and `obj2`, and the lung scANVI model over `obj3`, in one pipeline run. Start from [examples/model_sheet.csv](examples/model_sheet.csv) and [examples/annotate-with-models.params.yml](examples/annotate-with-models.params.yml).

| column | required | meaning |
| --- | --- | --- |
| `tool` | yes | `celltypist`, `scanvi`, or `panhumanpy` |
| `model` | yes, except for `panhumanpy` | the `celltypist_model.pkl` file, or the `scanvi_model` directory. Left empty for PanHumanPy, which carries its own model. |
| `query_h5ad` | yes | one query path, or a glob that expands to one pair per matching file |
| `model_id` | no | names the pair's output directory. Defaults to the model file or directory name. |

Sheet rules:

- A tab anywhere in the header line makes the file a TSV; anything else is read as a CSV. Blank lines and lines starting with `#` are ignored, and fields are trimmed. There are no quoting rules, so a path cannot contain the delimiter.
- Each row is independent. One row means one tool applied to one query, so a query listed twice is annotated twice and its two results never mix.
- **Set `model_id` whenever two models appear in one sheet.** Every model this pipeline trains is called `celltypist_model.pkl` or `scanvi_model`, so the default ID collides as soon as two training runs are combined, and the run stops with an error rather than publishing one pair over another.
- The whole sheet is parsed and checked before any job is submitted: unknown tools, missing or unexpected columns, a model path that does not exist, a `scanvi` row pointing at a file instead of a directory, a glob that matched nothing, a duplicated pair, and IDs that cannot be used in a path all fail immediately.
- The sheet is the only place tools, models, and queries are named for this entry point, so `--tools`, `--query_h5ad`, `--ref_h5ad`, `--ref_label_col`, `--ref_filter_*`, `--celltypist_model`, `--scanvi_model_dir`, and `--scvi_reference_model_dir` are not used. Supplying one of them is an error rather than a silent no-op.
- No reference takes part. Each model already carries the gene space it was trained on, and each annotation step reports its own overlap against that space, so the reference-coverage diagnostic is not available here.

Each pair publishes to `<outdir>/<sample_id>.<model_id>/`, so two models over one query keep separate results. Because each pair runs a single tool, merging is **off by default** for this entry point: a run publishes the per-tool prediction tables and nothing else. Pass `--merge_annotations true` to also write the annotated `.h5ad`, cluster summary, and marker genes for every pair.

## Skipping the Merge Step

`--merge_annotations` controls whether per-tool predictions are written back into the query object:

- `true` merges them into `<project_tag>.<id>.annotated.h5ad` and writes the cluster annotation summary and the marker-gene outputs.
- `false` skips that step. The per-tool prediction tables are still published, and `--query_cluster_col` is no longer required in the query, since only the merge step reads it.

It can be set in the `-params-file` YAML as well as on the command line:

```yaml
merge_annotations: false
```

Left unset it defaults to `true` for the default workflow and `false` for `ANNOTATE_WITH_MODELS`. Every run logs which mode it resolved to.

Turning it off is worth considering for large queries. The merge step is by far the most expensive part of a run: it holds a second copy of the query matrix in memory to rank marker genes, single-threaded, and on a 66 GB query it has taken over 12 hours and ~190 GB of RAM. If you do want merged outputs for an object that size, give the process a longer queue:

```groovy
// merge.config
process {
  withName: MERGE_ANNOTATIONS {
    queue = 'long'
  }
}
```

```bash
nextflow run main.nf -c merge.config -params-file your-run.params.yml
```

The parameters that only affect merging - `--query_cluster_col`, `--compute_missing_clusters`, `--marker_genes_*`, and `--marker_dotplot*` - are ignored when the step is skipped.

## Tool-Specific Parameters

The pipeline selects model behavior from artifact flags instead of explicit mode flags:

- CellTypist trains from the reference unless `--celltypist_model` is supplied.
- scANVI loads `--scanvi_model_dir` when supplied; otherwise it trains a reference scANVI model from labeled reference cells, optionally initialized from `--scvi_reference_model_dir`, then maps the query with scANVI query mapping. This query-mapping step uses the scArches-style transfer-learning support in scvi-tools.
- PanHumanPy always annotates the query directly.

Required when using pretrained or external artifacts:

- `--celltypist_model`: required to reuse a pretrained CellTypist model.
- `--scanvi_model_dir`: required to reuse a pretrained scANVI model directory. It must contain the `scanvi_model_metadata.json` written by this pipeline, which records the model's gene space; a directory saved by scvi-tools alone does not carry that information and is rejected.
- `--scvi_reference_model_dir`: required only when initializing scANVI training from a pretrained scVI reference model. Do not combine this with `--scanvi_model_dir`.

CellTypist optional controls:

- `--celltypist_training_mode`: `standard` or `detailed`. Default: `detailed`.
- `--celltypist_balance_cells_per_label`: maximum reference cells per label in detailed mode. Default: `500`.
- `--celltypist_feature_selection`: run tutorial-style feature selection in detailed mode. Default: `true`.
- `--celltypist_feature_selection_top_genes`: top genes per label for feature selection. Default: `100`.
- `--celltypist_feature_selection_max_iter`: quick SGD feature-selection iterations. Default: `5`.
- `--celltypist_max_iter`: final CellTypist fit iterations. Default: `100`.
- `--celltypist_use_gpu`: train CellTypist models on the GPU queue. Default: `true`. Set `false` to train on CPU queues. On GPU, the final model fit uses cuML logistic regression; feature selection still runs on CPU, because CellTypist fits it with SGD. Training fails if the container lacks cuML, rather than falling back to CPU while holding a GPU. GPU and CPU training can give slightly different models, because the logistic regression solvers differ. Annotation always runs on CPU: CellTypist predicts on CPU regardless, and its GPU option only covers the majority-voting over-clustering, whose GPU path fails in CellTypist 1.7.1 when it has to build the neighbour graph itself.

CellTypist model reuse checks that at least 80% of model genes are present in the query. This catches gene identifier mismatches before silently filling most model features with zeros.

scANVI optional controls:

- `--batch_key`: optional column present in both reference and query `obs`; omitted values stay as real `None`/`null`. Use this when batches should be modeled.
- `--scvi_layer`: optional AnnData layer passed to `SCVI.setup_anndata`.
- `--scanvi_categorical_covariate_keys`: optional comma-separated additional categorical covariate columns present in both objects.
- `--scanvi_continuous_covariate_keys`: optional comma-separated continuous covariate columns present in both objects.
- `--scanvi_hvg_n_top_genes`: HVGs for new scVI/scANVI training. Set `0` to keep all shared genes. Default: `2500`.
- `--scvi_n_hidden`, `--scvi_n_latent`, `--scvi_n_layers`, `--scvi_dropout_rate`: optional scVI architecture controls. Unset values use scvi-tools defaults.
- `--scvi_dispersion`, `--scvi_gene_likelihood`, `--scvi_latent_distribution`: optional scVI model distribution controls.
- `--scvi_train_max_epochs`, `--scvi_train_batch_size`, `--scvi_train_train_size`, `--scvi_train_validation_size`, `--scvi_train_early_stopping`, `--scvi_train_early_stopping_patience`: optional `SCVI.train` controls.
- `--scanvi_train_max_epochs`, `--scanvi_train_batch_size`, `--scanvi_train_train_size`, `--scanvi_train_validation_size`, `--scanvi_train_early_stopping`, `--scanvi_train_early_stopping_patience`: optional reference `SCANVI.train` controls.
- `--scanvi_query_max_epochs`, `--scanvi_query_batch_size`, `--scanvi_query_train_size`, `--scanvi_query_validation_size`, `--scanvi_query_early_stopping`, `--scanvi_query_early_stopping_patience`, `--scanvi_query_check_val_every_n_epoch`, `--scanvi_query_plan_weight_decay`: optional query-mapping `SCANVI.train` controls. `--scanvi_query_max_epochs` defaults to `100`; the other values preserve current script defaults when unset.
- `--unknown_label`: unlabeled category used by scANVI during query mapping. Default: `Unknown`.

There are two ways to configure scVI/scANVI training:

1. Use native nf-autoannotate params, usually through `-params-file examples/autoannotate.params.yml`. This is the recommended mode for new runs.
2. Use an nf-scVI-metrics-style Python config with `--nf_scvi_metrics_config config.py`. This compatibility mode is for reusing existing nf-scVI-metrics config files without adding that Python file path to the YAML params file.

`--nf_scvi_metrics_config` accepts the same top-level `model_input` and `param_input` dictionaries used by `nf-scVI-metrics/examples/inputs.py`. Because nf-autoannotate runs one model configuration per pipeline run, each model/training parameter must resolve to one value. Single-item lists such as `n_latent = [30]` are accepted and unwrapped; sweep lists such as `n_layers = [1, 2]` are rejected during input validation.

Supported `model_input` keys are `n_hidden`, `n_latent`, `n_layers`, `dropout_rate`, `dispersion`, `gene_likelihood`, `latent_distribution`, `max_epochs`, `accelerator`, `devices`, `train_size`, `validation_size`, `batch_size`, `early_stopping`, `early_stopping_patience`, and `counts_layer`. `counts_layer` is mapped to the scvi-tools `layer` setup argument. nf-scVI-metrics dataset selectors such as `adata_path` are ignored because nf-autoannotate uses `--ref_h5ad` and `--query_h5ad`.

Supported `param_input` keys are `batch_key`, `categorical_covariate_keys`, and `continuous_covariate_keys`. `layer` is also accepted for native scvi-tools naming. `categorical_covariate_keys` and `continuous_covariate_keys` may be lists because they are single scvi-tools parameter values, not sweep lists.

Two optional nf-autoannotate-only sections are also supported: `scanvi_train_input` for reference scANVI fine-tuning, and `query_train_input` for query-mapping `SCANVI.train`. If `scanvi_train_input` is omitted, the `model_input` training parameters are reused. In `query_train_input`, use `plan_weight_decay` for query-training weight decay.

When `--nf_scvi_metrics_config` is supplied, individual Nextflow params for setup/model/train/query-train settings are intentionally ignored so the Python config file is the single source of truth for scVI/scANVI settings. The config file is checked during input validation, before scVI/scANVI training starts. The old `--scvi_training_config` parameter is no longer supported.

PanHumanPy optional controls:

- `--panhumanpy_feature_names_col`: optional query `var` column containing gene symbols when `query.var_names` are not gene symbols.

Output naming and merge behavior:

- `--outdir`: output directory. Default: `autoannotate-results-<project_tag>`.
- `--merge_annotations`: merge the per-tool predictions into the query `.h5ad` and write the cluster summary and marker genes. Defaults to `true` for the default workflow and `false` for `-entry ANNOTATE_WITH_MODELS`. See [Skipping the Merge Step](#skipping-the-merge-step).
- `--scanvi_obs_col`: scANVI label column name in output. Default: `scanvi_label`.
- `--scanvi_score_col`: scANVI confidence column name in output. Default: `scanvi_confidence`.
- `--query_cluster_col`: existing query `obs` cluster column used for cluster summaries. Only read by the merge step. Default: `leiden`.
- `--compute_missing_clusters`: compute Leiden clusters in the merge step if `--query_cluster_col` is absent. Default: `false`.
- `--marker_genes_n_top`: number of marker genes to report per cluster. Default: `50`.
- `--marker_genes_method`: differential-expression test for marker ranking, either `wilcoxon` or `t-test`. Default: `wilcoxon`.
- `--marker_dotplot`: write a marker-gene dotplot PNG from the ranked marker table. Default: `true`.
- `--marker_dotplot_n_top`: number of ranked marker genes per cluster to include in the dotplot. Default: `5`.

When missing clusters are explicitly computed, the merge step reuses existing neighbors or `X_pca` if present. Otherwise it runs PCA/neighbors/Leiden on a normalized, log-transformed copy so the final query matrix is not modified.

Parameter defaults live in [nextflow.config](nextflow.config); user-facing parameter validation lives in [nextflow_schema.json](nextflow_schema.json).

## Validation

Validation is split in two. Reference validation runs once before any model step and checks that:

- the reference file exists
- `ref_label_col` exists in `ref_h5ad.obs`
- `batch_key`, if supplied, exists in `ref_h5ad.obs`
- all `scanvi_categorical_covariate_keys` and `scanvi_continuous_covariate_keys`, if supplied, exist in the reference
- reference filters are valid and keep at least one reference cell
- reference cell and gene indices are unique

Query validation runs once per query and checks that:

- the query file exists
- `query_h5ad.X` contains raw non-negative integer counts
- query cell and gene indices are unique
- `query_cluster_col` exists in `query_h5ad.obs`, unless `compute_missing_clusters` is enabled or the merge step is skipped
- `batch_key` and any covariate keys, if supplied, exist in the query `obs`
- `panhumanpy_feature_names_col`, if supplied, exists in `query_h5ad.var`
- the query and reference gene identifiers use the same broad scheme, such as symbols or Ensembl IDs
- the query shares at least one gene with the reference

The whole reference stage is skipped when no reference h5ad is used, so a `--tools panhumanpy` run, a run annotating purely from pretrained models, and every `-entry ANNOTATE_WITH_MODELS` run validate their queries alone. The last four checks above only apply when a reference is used.

The merge step also validates that every annotation CSV has unique `cell_id` values matching `query.obs_names` exactly. Missing or extra cells now fail the merge instead of becoming silent `NaN` values.

`-entry ANNOTATE_WITH_MODELS` validates its sheet before submitting any job; see [Annotating From Existing Models](#annotating-from-existing-models) for what it checks.

Neither stage loads the AnnData object. Each opens the h5ad with `h5py` and reads the cell and gene indices plus the names of the `obs` and `var` columns and of the layers. The only column values read are those of `ref_filter_col`, when reference filtering is requested. The raw-count check reads `X` in chunks, or samples rows of it for a CellTypist-only run. Validation memory therefore tracks the number of cells, not how many layers, embeddings, or metadata columns the file carries. On a 758,815-cell query, query validation peaked under 600 MB, against 33 GB when the file was opened with `anndata.read_h5ad(backed="r")`.

## Outputs

Shared reference and model artifacts are published at the top of the output directory, and each query gets its own subdirectory named after its sample ID:

```
<outdir>/
  reference_manifest.json
  reference_genes.txt
  celltypist_model.pkl
  celltypist_model_metadata.json
  scvi_model/
  scvi_model_metadata.json
  scanvi_model/
  scanvi_model_metadata.json
  donor1/
    donor1.query_manifest.json
    <project_tag>.donor1.celltypist_predictions.csv
    <project_tag>.donor1.scanvi_predictions.csv
    <project_tag>.donor1.scanvi_latent.h5ad
    <project_tag>.donor1.panhumanpy_predictions.csv
    <project_tag>.donor1.annotated.h5ad
    <project_tag>.donor1.cluster_annotation_summary.csv
    <project_tag>.donor1.cluster_marker_genes.tsv
    <project_tag>.donor1.cluster_marker_dotplot.png
  donor2/
    ...
```

The sample ID appears in both the directory name and the filenames so that results stay distinguishable when copied out of the tree. Single-query runs use the same layout.

Each selected tool's own prediction table is published next to the merged object, so a run that skips the merge step still has deliverables. With `--merge_annotations false` the same directory holds the query manifest and the prediction tables alone.

`-entry ANNOTATE_WITH_MODELS` uses one directory per (model, query) pair, named `<sample_id>.<model_id>`, and repeats the pair ID in every filename:

```
<outdir>/
  obj1.gut_ct/
    obj1.gut_ct.query_manifest.json
    <project_tag>.obj1.gut_ct.celltypist_predictions.csv
  obj2.gut_ct/
    ...
  obj3.lung_scanvi/
    obj3.lung_scanvi.query_manifest.json
    <project_tag>.obj3.lung_scanvi.scanvi_predictions.csv
    <project_tag>.obj3.lung_scanvi.scanvi_latent.h5ad
```

No model or reference artifacts are published there, because that entry point trains nothing.

The marker dotplot is written when `--marker_dotplot true` and at least two clusters are available.

The annotated `.h5ad` is the query object with additional `obs` columns and scANVI latent representations in `obsm`.

The lists below describe a run with all three tools selected. A subset run only produces the columns and artifacts belonging to the tools named in `--tools`.

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

Latent representations copied into `query.obsm`:

- `X_scanvi`: latent representation returned by `scvi-tools` `SCANVI.get_latent_representation` during query annotation.
- `X_scvi`: compatibility alias of the same scANVI latent representation for downstream tools that expect the scVI-style key.

The latent arrays are written by `annotate_scanvi.py` to a separate `scanvi_latent.h5ad` artifact and merged into the final annotated `.h5ad`; the `scanvi_predictions.csv` table stays limited to per-cell labels and confidence scores.

Comparison columns created in `query.obs`, only when `--tools` includes both `celltypist` and `scanvi`:

- `celltypist_scanvi_agree`: per-cell boolean; `true` only when both methods have non-missing labels and `celltypist_majority_voting` equals the scANVI label.
- `celltypist_scanvi_cluster_agreement_fraction`: cluster-level fraction copied onto each cell in the cluster; the fraction of cells where `celltypist_scanvi_agree` is `true`.

The cluster summary is grouped by `--query_cluster_col`. Its columns depend on `--tools`: `<query_cluster_col>` and `n_cells` are always present, each selected tool contributes its own pair of columns, and the two `celltypist_scanvi_*` columns need both CellTypist and scANVI. With all three tools selected it includes:

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

If the query contains fewer than two clusters, the marker-gene TSV is still written with headers but no marker rows because there is no between-cluster comparison. The marker dotplot is skipped in that case.

When `--marker_dotplot true`, the merge step writes `<project_tag>.cluster_marker_dotplot.png` with Scanpy `dotplot` using the top `--marker_dotplot_n_top` ranked genes per cluster from the same marker calculation. Set `--marker_dotplot false` to skip the PNG output while keeping the marker-gene TSV unchanged.

A subset run publishes only the artifacts of the tools it ran. `reference_manifest.json` and `reference_genes.txt` are omitted when no reference-based tool is selected.

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

Use the generated files for a demo run by editing a params file:

```bash
cp examples/autoannotate.params.yml demo.params.yml
```

```yaml
ref_h5ad: data/demo/reference.h5ad
query_h5ad: data/demo/query.h5ad
project_tag: demo
ref_label_col: bulk_labels
query_cluster_col: louvain
outdir: results/demo
```

```bash
module load cellgen/nextflow
nextflow run main.nf -params-file demo.params.yml
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

Annotate several queries against one set of trained models:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad '/path/to/queries/*.h5ad' \
  --project_tag cohort \
  --ref_label_col cell_type \
  --batch_key donor \
  --query_cluster_col leiden
```

Train CellTypist plus scVI/scANVI model artifacts only, with no query at all:

```bash
module load cellgen/nextflow

nextflow run main.nf -entry TRAIN_MODELS \
  --ref_h5ad /path/to/reference.h5ad \
  --project_tag model_training \
  --ref_label_col cell_type \
  --batch_key donor \
  --outdir results/model_training
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

Run PanHumanPy only, with no reference dataset:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --tools panhumanpy \
  --query_h5ad /path/to/query.h5ad \
  --project_tag panhumanpy_only \
  --query_cluster_col leiden
```

Run CellTypist and PanHumanPy, skipping the GPU scANVI steps:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --tools celltypist,panhumanpy \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag no_scanvi \
  --ref_label_col cell_type
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

Annotate with CellTypist only and skip the merge step, publishing the predictions alone:

```bash
module load cellgen/nextflow

nextflow run main.nf \
  --tools celltypist \
  --ref_h5ad /path/to/reference.h5ad \
  --query_h5ad /path/to/query.h5ad \
  --project_tag predictions_only \
  --ref_label_col cell_type \
  --merge_annotations false
```

Apply different existing models to different queries in one run:

```bash
module load cellgen/nextflow

nextflow run main.nf -entry ANNOTATE_WITH_MODELS \
  --model_sheet /path/to/model_sheet.csv \
  --project_tag reuse_run \
  --outdir results/reuse_run
```

The same run, with each pair merged back into its own annotated object:

```bash
module load cellgen/nextflow

nextflow run main.nf -entry ANNOTATE_WITH_MODELS \
  --model_sheet /path/to/model_sheet.csv \
  --project_tag reuse_run \
  --outdir results/reuse_run \
  --query_cluster_col leiden \
  --merge_annotations true
```

## Assumptions

- Input `.h5ad` files should already be normalized/preprocessed appropriately for the chosen annotation methods.
- CellTypist and scANVI are trained from the reference alone, never from a reference-query intersection, so a single training run serves any number of queries. Feature selection and HVG selection then narrow each model's feature space to a subset of the reference. At annotation time each query is aligned to the trained model's own feature list and genes it lacks are zero-filled; scANVI query mapping otherwise follows the saved model registry.
- Pretrained scANVI directories must carry the `scanvi_model_metadata.json` written by this pipeline, which records the model's gene space. A reference gene list cannot substitute for it, because the model's features are a post-HVG subset.
- Every query is annotated independently. No query influences another's labels, and no cross-query summary is produced. In `-entry ANNOTATE_WITH_MODELS` the same holds per (model, query) pair: two models over one query produce two independent sets of results.
- PanHumanPy receives the original query feature space, optionally using `--panhumanpy_feature_names_col` for symbols.
- CellTypist, scANVI, and PanHumanPy outputs are kept separate, and each can be enabled or disabled independently with `--tools`.
- PanHumanPy is not compared with CellTypist or scANVI because it uses a different reference.
