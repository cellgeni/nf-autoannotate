#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def celltypistMode = params.celltypist_model ? 'pretrained' : 'train'
def scanviMode = params.scanvi_model_dir ? 'pretrained' : 'train'

def refH5ad = file(params.ref_h5ad, checkIfExists: true)
def queryH5ad = file(params.query_h5ad, checkIfExists: true)
def scviReferenceModelDir = params.scvi_reference_model_dir ? file(params.scvi_reference_model_dir, checkIfExists: true) : null
def scviTrainingConfig = params.scvi_training_config ? file(params.scvi_training_config, checkIfExists: true) : null
def useScviTrainingConfig = params.scvi_training_config ? true : false
def outputDir = params.outdir?.toString()?.trim() ?: "autoannotate-results-${params.project_tag}"
def refFilterCol = params.ref_filter_col
def refFilterValues = params.ref_filter_values
def unsupportedFilterParams = ["ref_timepoint_col", "ref_timepoints"].findAll { params.containsKey(it) }

if( unsupportedFilterParams ) {
    error "Unsupported parameter(s): ${unsupportedFilterParams.collect { "--${it}" }.join(", ")}. Use --ref_filter_col and --ref_filter_values instead."
}
if( !params.project_tag ) {
    error "Missing required parameter: --project_tag"
}
if( !params.ref_label_col ) {
    error "Missing required parameter: --ref_label_col"
}
if( (refFilterCol && !refFilterValues) || (!refFilterCol && refFilterValues) ) {
    error "Use --ref_filter_col and --ref_filter_values together when filtering the reference"
}
if( !['standard', 'detailed'].contains(params.celltypist_training_mode) ) {
    error "Invalid --celltypist_training_mode '${params.celltypist_training_mode}'. Expected one of: standard, detailed"
}
if( params.scanvi_hvg_n_top_genes.toInteger() < 0 ) {
    error "Invalid --scanvi_hvg_n_top_genes '${params.scanvi_hvg_n_top_genes}'. Expected a non-negative integer"
}
if( params.scanvi_model_dir && params.scvi_reference_model_dir ) {
    error "Use either --scanvi_model_dir for a pretrained scANVI model or --scvi_reference_model_dir to initialize scANVI training from a pretrained scVI reference model, not both"
}

process VALIDATE_INPUTS {
    label 'validation'

    input:
    path ref_h5ad
    path query_h5ad

    output:
    path 'validation_manifest.json', emit: manifest
    path 'shared_genes.txt', emit: shared_genes

    script:
    """
    python3 ${projectDir}/scripts/validate_inputs.py \
        --ref-h5ad ${ref_h5ad} \
        --query-h5ad ${query_h5ad} \
        --ref-label-col '${params.ref_label_col}' \
        ${(!useScviTrainingConfig && params.batch_key) ? "--batch-key '${params.batch_key}'" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_categorical_covariate_keys) ? "--categorical-covariate-keys '${params.scanvi_categorical_covariate_keys}'" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_continuous_covariate_keys) ? "--continuous-covariate-keys '${params.scanvi_continuous_covariate_keys}'" : ''} \
        ${params.panhumanpy_feature_names_col ? "--panhumanpy-feature-names-col '${params.panhumanpy_feature_names_col}'" : ''} \
        ${refFilterCol ? "--ref-filter-col '${refFilterCol}'" : ''} \
        ${refFilterValues ? "--ref-filter-values '${refFilterValues}'" : ''} \
        --output validation_manifest.json \
        --shared-genes-output shared_genes.txt
    """
}

process CELLTYPIST_TRAIN_MODEL {
    label 'cpu4'

    input:
    path ref_h5ad
    path validation_manifest
    path shared_genes

    output:
    path 'celltypist_model.pkl', emit: model
    path 'celltypist_model_metadata.json', emit: metadata

    script:
    """
    python3 ${projectDir}/scripts/train_celltypist.py \
        --ref-h5ad ${ref_h5ad} \
        --validation-manifest ${validation_manifest} \
        --shared-genes ${shared_genes} \
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

process CELLTYPIST_ANNOTATE {
    label 'normal'

    input:
    path query_h5ad
    path shared_genes
    path celltypist_model

    output:
    path 'celltypist_predictions.csv', emit: predictions

    script:
    """
    python3 ${projectDir}/scripts/annotate_celltypist.py \
        --query-h5ad ${query_h5ad} \
        --shared-genes ${shared_genes} \
        --model-path ${celltypist_model} \
        --output celltypist_predictions.csv
    """
}

process SCANVI_TRAIN_MODEL {
    label 'gpu'

    input:
    path ref_h5ad
    path query_h5ad
    path validation_manifest
    path shared_genes

    output:
    path 'scanvi_model', emit: model_dir
    path 'scanvi_model_metadata.json', emit: metadata

    script:
    """
    python3 ${projectDir}/scripts/train_scanvi.py \
        --ref-h5ad ${ref_h5ad} \
        --query-h5ad ${query_h5ad} \
        --validation-manifest ${validation_manifest} \
        --shared-genes ${shared_genes} \
        --ref-label-col '${params.ref_label_col}' \
        ${(!useScviTrainingConfig && params.batch_key) ? "--batch-key '${params.batch_key}'" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_categorical_covariate_keys) ? "--categorical-covariate-keys '${params.scanvi_categorical_covariate_keys}'" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_continuous_covariate_keys) ? "--continuous-covariate-keys '${params.scanvi_continuous_covariate_keys}'" : ''} \
        ${(!useScviTrainingConfig && params.scvi_layer) ? "--layer '${params.scvi_layer}'" : ''} \
        --unknown-label '${params.unknown_label}' \
        --hvg-n-top-genes ${params.scanvi_hvg_n_top_genes} \
        ${useScviTrainingConfig ? "--training-config ${scviTrainingConfig}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_n_hidden != null) ? "--scvi-n-hidden ${params.scvi_n_hidden}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_n_latent != null) ? "--scvi-n-latent ${params.scvi_n_latent}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_n_layers != null) ? "--scvi-n-layers ${params.scvi_n_layers}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_dropout_rate != null) ? "--scvi-dropout-rate ${params.scvi_dropout_rate}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_dispersion) ? "--scvi-dispersion '${params.scvi_dispersion}'" : ''} \
        ${(!useScviTrainingConfig && params.scvi_gene_likelihood) ? "--scvi-gene-likelihood '${params.scvi_gene_likelihood}'" : ''} \
        ${(!useScviTrainingConfig && params.scvi_latent_distribution) ? "--scvi-latent-distribution '${params.scvi_latent_distribution}'" : ''} \
        ${(!useScviTrainingConfig && params.scvi_train_max_epochs != null) ? "--scvi-train-max-epochs ${params.scvi_train_max_epochs}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_train_batch_size != null) ? "--scvi-train-batch-size ${params.scvi_train_batch_size}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_train_train_size != null) ? "--scvi-train-train-size ${params.scvi_train_train_size}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_train_validation_size != null) ? "--scvi-train-validation-size ${params.scvi_train_validation_size}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_train_early_stopping != null) ? "--scvi-train-early-stopping ${params.scvi_train_early_stopping}" : ''} \
        ${(!useScviTrainingConfig && params.scvi_train_early_stopping_patience != null) ? "--scvi-train-early-stopping-patience ${params.scvi_train_early_stopping_patience}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_train_max_epochs != null) ? "--scanvi-train-max-epochs ${params.scanvi_train_max_epochs}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_train_batch_size != null) ? "--scanvi-train-batch-size ${params.scanvi_train_batch_size}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_train_train_size != null) ? "--scanvi-train-train-size ${params.scanvi_train_train_size}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_train_validation_size != null) ? "--scanvi-train-validation-size ${params.scanvi_train_validation_size}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_train_early_stopping != null) ? "--scanvi-train-early-stopping ${params.scanvi_train_early_stopping}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_train_early_stopping_patience != null) ? "--scanvi-train-early-stopping-patience ${params.scanvi_train_early_stopping_patience}" : ''} \
        --output-model-dir scanvi_model \
        --output-metadata scanvi_model_metadata.json \
        ${params.scvi_reference_model_dir ? "--scvi-model-dir ${scviReferenceModelDir}" : ''}
    """
}

process SCANVI_ANNOTATE {
    label 'gpu'

    input:
    path query_h5ad
    path shared_genes
    path scanvi_model_dir

    output:
    path 'scanvi_predictions.csv', emit: predictions

    script:
    """
    python3 ${projectDir}/scripts/annotate_scanvi.py \
        --query-h5ad ${query_h5ad} \
        --shared-genes ${shared_genes} \
        --model-dir ${scanvi_model_dir} \
        --label-column '${params.scanvi_obs_col}' \
        --score-column '${params.scanvi_score_col}' \
        ${useScviTrainingConfig ? "--training-config ${scviTrainingConfig}" : "--query-max-epochs ${params.scanvi_query_max_epochs}"} \
        ${(!useScviTrainingConfig && params.scanvi_query_batch_size != null) ? "--query-batch-size ${params.scanvi_query_batch_size}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_query_train_size != null) ? "--query-train-size ${params.scanvi_query_train_size}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_query_validation_size != null) ? "--query-validation-size ${params.scanvi_query_validation_size}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_query_early_stopping != null) ? "--query-early-stopping ${params.scanvi_query_early_stopping}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_query_early_stopping_patience != null) ? "--query-early-stopping-patience ${params.scanvi_query_early_stopping_patience}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_query_check_val_every_n_epoch != null) ? "--query-check-val-every-n-epoch ${params.scanvi_query_check_val_every_n_epoch}" : ''} \
        ${(!useScviTrainingConfig && params.scanvi_query_plan_weight_decay != null) ? "--query-plan-weight-decay ${params.scanvi_query_plan_weight_decay}" : ''} \
        --output scanvi_predictions.csv
    """
}

process PANHUMANPY_ANNOTATE {
    label 'normal'

    input:
    path query_h5ad

    output:
    path 'panhumanpy_predictions.csv', emit: predictions

    script:
    """
    python3 ${projectDir}/scripts/annotate_panhumanpy.py \
        --query-h5ad ${query_h5ad} \
        ${params.panhumanpy_feature_names_col ? "--feature-names-col '${params.panhumanpy_feature_names_col}'" : ''} \
        --output panhumanpy_predictions.csv
    """
}

process MERGE_ANNOTATIONS {
    label 'normal'
    publishDir outputDir, mode: 'copy', overwrite: true

    input:
    path query_h5ad
    path celltypist_csv
    path scanvi_csv
    path panhumanpy_csv

    output:
    path "${params.project_tag}.annotated.h5ad", emit: h5ad
    path "${params.project_tag}.cluster_annotation_summary.csv", emit: cluster_summary

    script:
    """
    python3 ${projectDir}/scripts/merge_annotations.py \
        --query-h5ad ${query_h5ad} \
        --celltypist-csv ${celltypist_csv} \
        --scanvi-csv ${scanvi_csv} \
        --panhumanpy-csv ${panhumanpy_csv} \
        --cluster-column '${params.query_cluster_col}' \
        --scanvi-label-column '${params.scanvi_obs_col}' \
        --scanvi-score-column '${params.scanvi_score_col}' \
        ${params.compute_missing_clusters ? '--compute-missing-clusters' : ''} \
        --output ${params.project_tag}.annotated.h5ad \
        --cluster-summary-output ${params.project_tag}.cluster_annotation_summary.csv
    """
}

workflow {
    Channel
        .value(refH5ad)
        .set { ref_h5ad_ch }

    Channel
        .value(queryH5ad)
        .set { query_h5ad_ch }

    validation_outputs = VALIDATE_INPUTS(ref_h5ad_ch, query_h5ad_ch)
    validation_manifest_ch = validation_outputs.manifest
    shared_genes_ch = validation_outputs.shared_genes

    celltypist_model_ch = celltypistMode == 'pretrained'
        ? Channel.value(file(params.celltypist_model, checkIfExists: true))
        : CELLTYPIST_TRAIN_MODEL(ref_h5ad_ch, validation_manifest_ch, shared_genes_ch).model

    scanvi_model_ch = scanviMode == 'pretrained'
        ? Channel.value(file(params.scanvi_model_dir, checkIfExists: true))
        : SCANVI_TRAIN_MODEL(ref_h5ad_ch, query_h5ad_ch, validation_manifest_ch, shared_genes_ch).model_dir

    celltypist_predictions_ch = CELLTYPIST_ANNOTATE(query_h5ad_ch, shared_genes_ch, celltypist_model_ch).predictions
    scanvi_predictions_ch = SCANVI_ANNOTATE(query_h5ad_ch, shared_genes_ch, scanvi_model_ch).predictions
    panhumanpy_predictions_ch = PANHUMANPY_ANNOTATE(query_h5ad_ch).predictions

    MERGE_ANNOTATIONS(query_h5ad_ch, celltypist_predictions_ch, scanvi_predictions_ch, panhumanpy_predictions_ch)
}
