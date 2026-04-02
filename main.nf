#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def celltypistMode = params.celltypist_mode?.toString()?.trim()
def scanviMode = params.scanvi_mode?.toString()?.trim()

if( !celltypistMode ) {
    celltypistMode = params.celltypist_model ? 'pretrained' : 'train'
}
if( !scanviMode ) {
    scanviMode = params.scanvi_model_dir ? 'pretrained' : 'train'
}

def refH5ad = file(params.ref_h5ad, checkIfExists: true)
def queryH5ad = file(params.query_h5ad, checkIfExists: true)
def scviReferenceModelDir = params.scvi_reference_model_dir ? file(params.scvi_reference_model_dir, checkIfExists: true) : null

if( !params.project_tag ) {
    error "Missing required parameter: --project_tag"
}
if( !params.ref_label_col ) {
    error "Missing required parameter: --ref_label_col"
}
if( !params.batch_key ) {
    error "Missing required parameter: --batch_key"
}
if( (params.ref_timepoint_col && !params.ref_timepoints) || (!params.ref_timepoint_col && params.ref_timepoints) ) {
    error "Use --ref_timepoint_col and --ref_timepoints together when filtering the reference by timepoint"
}
if( !['standard', 'detailed'].contains(params.celltypist_training_mode) ) {
    error "Invalid --celltypist_training_mode '${params.celltypist_training_mode}'. Expected one of: standard, detailed"
}
if( !['train', 'pretrained'].contains(celltypistMode) ) {
    error "Invalid --celltypist_mode '${params.celltypist_mode}'. Expected one of: train, pretrained"
}
if( !['train', 'pretrained'].contains(scanviMode) ) {
    error "Invalid --scanvi_mode '${params.scanvi_mode}'. Expected one of: train, pretrained"
}
if( celltypistMode == 'pretrained' && !params.celltypist_model ) {
    error "--celltypist_model is required when --celltypist_mode=pretrained"
}
if( scanviMode == 'pretrained' && !params.scanvi_model_dir ) {
    error "--scanvi_model_dir is required when --scanvi_mode=pretrained"
}
if( scanviMode == 'pretrained' && params.scvi_reference_model_dir ) {
    error "Use either --scanvi_model_dir for a pretrained scANVI model or --scvi_reference_model_dir to initialize scANVI training from a pretrained scVI reference model, not both"
}

process VALIDATE_INPUTS {
    label 'normal'

    input:
    path ref_h5ad
    path query_h5ad

    output:
    path 'validation_manifest.json', emit: manifest

    script:
    """
    python3 ${projectDir}/scripts/validate_inputs.py \
        --ref-h5ad ${ref_h5ad} \
        --query-h5ad ${query_h5ad} \
        --ref-label-col '${params.ref_label_col}' \
        --batch-key '${params.batch_key}' \
        ${params.ref_timepoint_col ? "--ref-timepoint-col '${params.ref_timepoint_col}'" : ''} \
        ${params.ref_timepoints ? "--ref-timepoints '${params.ref_timepoints}'" : ''} \
        --output validation_manifest.json
    """
}

process CELLTYPIST_TRAIN_MODEL {
    label 'normal4core'

    input:
    path ref_h5ad
    path validation_manifest

    output:
    path 'celltypist_model.pkl', emit: model
    path 'celltypist_model_metadata.json', emit: metadata

    script:
    """
    python3 ${projectDir}/scripts/train_celltypist.py \
        --ref-h5ad ${ref_h5ad} \
        --validation-manifest ${validation_manifest} \
        --ref-label-col '${params.ref_label_col}' \
        --output-model celltypist_model.pkl \
        --output-metadata celltypist_model_metadata.json \
        --n-jobs ${task.cpus} \
        --training-mode '${params.celltypist_training_mode}' \
        --balance-cells-per-label ${params.celltypist_balance_cells_per_label} \
        --feature-selection-top-genes ${params.celltypist_feature_selection_top_genes} \
        --feature-selection-max-iter ${params.celltypist_feature_selection_max_iter} \
        --max-iter ${params.celltypist_max_iter} \
        ${params.celltypist_feature_selection ? '--feature-selection' : ''}
    """
}

process CELLTYPIST_STAGE_PRETRAINED_MODEL {
    label 'normal'

    input:
    path celltypist_model

    output:
    path 'celltypist_model.pkl', emit: model

    script:
    """
    cp ${celltypist_model} celltypist_model.pkl
    """
}

process CELLTYPIST_ANNOTATE {
    label 'normal4core'

    input:
    path query_h5ad
    path validation_manifest
    path celltypist_model

    output:
    path 'celltypist_predictions.csv', emit: predictions

    script:
    """
    python3 ${projectDir}/scripts/annotate_celltypist.py \
        --query-h5ad ${query_h5ad} \
        --validation-manifest ${validation_manifest} \
        --model-path ${celltypist_model} \
        --output celltypist_predictions.csv
    """
}

process SCANVI_TRAIN_MODEL {
    label 'normalgpu'

    input:
    path ref_h5ad
    path query_h5ad
    path validation_manifest

    output:
    path 'scanvi_model', emit: model_dir
    path 'scanvi_model_metadata.json', emit: metadata

    script:
    """
    python3 ${projectDir}/scripts/train_scanvi.py \
        --ref-h5ad ${ref_h5ad} \
        --query-h5ad ${query_h5ad} \
        --validation-manifest ${validation_manifest} \
        --ref-label-col '${params.ref_label_col}' \
        --batch-key '${params.batch_key}' \
        --unknown-label '${params.unknown_label}' \
        --output-model-dir scanvi_model \
        --output-metadata scanvi_model_metadata.json \
        ${params.scvi_reference_model_dir ? "--scvi-model-dir ${scviReferenceModelDir}" : ''}
    """
}

process SCANVI_STAGE_PRETRAINED_MODEL {
    label 'normal'

    input:
    path scanvi_model_dir

    output:
    path 'scanvi_model', emit: model_dir

    script:
    """
    cp -r ${scanvi_model_dir} scanvi_model
    """
}

process SCANVI_ANNOTATE {
    label 'normalgpu'

    input:
    path query_h5ad
    path validation_manifest
    path scanvi_model_dir

    output:
    path 'scanvi_predictions.csv', emit: predictions

    script:
    """
    python3 ${projectDir}/scripts/annotate_scanvi.py \
        --query-h5ad ${query_h5ad} \
        --validation-manifest ${validation_manifest} \
        --model-dir ${scanvi_model_dir} \
        --label-column '${params.scanvi_obs_col}' \
        --score-column '${params.scanvi_score_col}' \
        --query-max-epochs ${params.scanvi_query_max_epochs} \
        --output scanvi_predictions.csv
    """
}

process MERGE_ANNOTATIONS {
    label 'normal'

    input:
    path query_h5ad
    path celltypist_csv
    path scanvi_csv

    output:
    path "${params.project_tag}.annotated.h5ad", emit: h5ad

    script:
    """
    python3 ${projectDir}/scripts/merge_annotations.py \
        --query-h5ad ${query_h5ad} \
        --celltypist-csv ${celltypist_csv} \
        --scanvi-csv ${scanvi_csv} \
        --scanvi-label-column '${params.scanvi_obs_col}' \
        --scanvi-score-column '${params.scanvi_score_col}' \
        --output ${params.project_tag}.annotated.h5ad
    """
}

workflow {
    Channel
        .value(refH5ad)
        .set { ref_h5ad_ch }

    Channel
        .value(queryH5ad)
        .set { query_h5ad_ch }

    validation_manifest_ch = VALIDATE_INPUTS(ref_h5ad_ch, query_h5ad_ch).manifest

    celltypist_model_ch = celltypistMode == 'pretrained'
        ? CELLTYPIST_STAGE_PRETRAINED_MODEL(Channel.value(file(params.celltypist_model, checkIfExists: true))).model
        : CELLTYPIST_TRAIN_MODEL(ref_h5ad_ch, validation_manifest_ch).model

    scanvi_model_ch = scanviMode == 'pretrained'
        ? SCANVI_STAGE_PRETRAINED_MODEL(Channel.value(file(params.scanvi_model_dir, checkIfExists: true))).model_dir
        : SCANVI_TRAIN_MODEL(ref_h5ad_ch, query_h5ad_ch, validation_manifest_ch).model_dir

    celltypist_predictions_ch = CELLTYPIST_ANNOTATE(query_h5ad_ch, validation_manifest_ch, celltypist_model_ch).predictions
    scanvi_predictions_ch = SCANVI_ANNOTATE(query_h5ad_ch, validation_manifest_ch, scanvi_model_ch).predictions

    MERGE_ANNOTATIONS(query_h5ad_ch, celltypist_predictions_ch, scanvi_predictions_ch)
}
