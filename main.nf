#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { validateParameters } from 'plugin/nf-schema'

validateParameters()

def knownTools = ['celltypist', 'scanvi', 'panhumanpy']
def selectedTools = (params.tools instanceof List ? params.tools : params.tools.toString().split(','))
    .collect { it.toString().trim().toLowerCase() }
    .findAll { it } as Set
def unknownTools = selectedTools - (knownTools as Set)
if( unknownTools ) {
    error "Unknown --tools entries: ${unknownTools.sort().join(', ')}. Expected one or more of: ${knownTools.join(', ')}"
}
if( !selectedTools ) {
    error "Select at least one annotation tool with --tools (${knownTools.join(', ')})"
}

def runCelltypist = selectedTools.contains('celltypist')
def runScanvi = selectedTools.contains('scanvi')
def runPanhumanpy = selectedTools.contains('panhumanpy')
def toolsArg = selectedTools.sort().join(',')

def celltypistMode = params.celltypist_model ? 'pretrained' : 'train'
def scanviMode = params.scanvi_model_dir ? 'pretrained' : 'train'

// A pretrained model carries its own feature space, so the reference is only strictly
// required when something is actually trained from it. It stays useful when supplied
// anyway, since it is what the per-query gene coverage diagnostic compares against.
def trainsCelltypist = runCelltypist && celltypistMode == 'train'
def trainsScanvi = runScanvi && scanviMode == 'train'
def needsReference = trainsCelltypist || trainsScanvi
def useReference = needsReference || (params.ref_h5ad as boolean)

if( !params.project_tag ) {
    error "Missing required parameter: --project_tag"
}

// Resolved from the parameter alone rather than from the tool selection, because
// ANNOTATE_WITH_MODELS never touches a reference and must not be held to the training
// rules. The entry points that do train check for a reference in their own body.
def refH5ad = params.ref_h5ad ? file(params.ref_h5ad, checkIfExists: true) : null

// query_h5ad accepts a single path, a glob, or a list of either. Resolve it eagerly so
// that an empty match or a sample-ID collision fails before any task is submitted. It is
// left unresolved for TRAIN_MODELS, which never reads a query, and for
// ANNOTATE_WITH_MODELS, which takes its queries from --model_sheet.
def sampleIdPattern = ~/^[A-Za-z0-9][A-Za-z0-9._-]*$/
def queryFiles = []
if( params.query_h5ad ) {
    def queryPatterns = params.query_h5ad instanceof Collection ? params.query_h5ad : [params.query_h5ad]
    queryFiles = queryPatterns.collectMany { pattern ->
        def resolved = file(pattern.toString(), checkIfExists: true)
        resolved instanceof Collection ? resolved.toList() : [resolved]
    }
    // Sorted so channel order, staging order, and -resume hashes do not depend on the
    // order the filesystem happened to return glob matches in.
    queryFiles = queryFiles.sort { it.toString() }
    if( !queryFiles ) {
        error "No query h5ad files matched --query_h5ad"
    }
    // Sample IDs are interpolated into output paths and shell arguments, so keep them to
    // characters that cannot escape a directory or break out of a quoted argument.
    def invalidSampleIds = queryFiles.collect { it.baseName }.findAll { !sampleIdPattern.matcher(it).matches() }
    if( invalidSampleIds ) {
        error "Query filenames produce unusable sample IDs: ${invalidSampleIds.sort().unique().join(', ')}. A sample ID is the filename without its final extension and must match ${sampleIdPattern}."
    }
    def duplicateSampleIds = queryFiles.groupBy { it.baseName }.findAll { id, files -> files.size() > 1 }
    if( duplicateSampleIds ) {
        error "Query files produce duplicate sample IDs: ${duplicateSampleIds.keySet().sort().join(', ')}. Sample IDs come from the filename without its extension, so rename the conflicting files."
    }
}

// Every per-query process is keyed by a meta map rather than a bare sample ID, because the
// unit of work is one query in the default workflow but one (model, query) pair in
// ANNOTATE_WITH_MODELS. `id` names the output directory and every output filename, and
// `tools` is what that single unit runs, which is narrower than --tools only in the
// model-sheet workflow.
def querySamples = queryFiles.collect { queryFile ->
    def sampleId = queryFile.baseName
    [[id: sampleId, sample_id: sampleId, tools: toolsArg], queryFile]
}
def scviReferenceModelDir = params.scvi_reference_model_dir ? file(params.scvi_reference_model_dir, checkIfExists: true) : null
def nfScviMetricsConfig = params.nf_scvi_metrics_config ? file(params.nf_scvi_metrics_config, checkIfExists: true) : null
def useNfScviMetricsConfig = params.nf_scvi_metrics_config ? true : false
def outputDir = params.outdir?.toString()?.trim() ?: "autoannotate-results-${params.project_tag}"
def refFilterCol = params.ref_filter_col
def refFilterValues = params.ref_filter_values

// Fail fast on options that only make sense for a tool that was not selected, rather
// than silently ignoring them and leaving the user to wonder why the tool never ran.
if( !runCelltypist && params.celltypist_model ) {
    error "--celltypist_model was provided but 'celltypist' is not in --tools"
}
if( !runScanvi && params.scanvi_model_dir ) {
    error "--scanvi_model_dir was provided but 'scanvi' is not in --tools"
}
if( !runScanvi && params.scvi_reference_model_dir ) {
    error "--scvi_reference_model_dir was provided but 'scanvi' is not in --tools"
}
if( !runScanvi && params.nf_scvi_metrics_config ) {
    error "--nf_scvi_metrics_config was provided but 'scanvi' is not in --tools"
}
if( !runPanhumanpy && params.panhumanpy_feature_names_col ) {
    error "--panhumanpy_feature_names_col was provided but 'panhumanpy' is not in --tools"
}

if( (refFilterCol && !refFilterValues) || (!refFilterCol && refFilterValues) ) {
    error "Use --ref_filter_col and --ref_filter_values together when filtering the reference"
}
if( runCelltypist && !['standard', 'detailed'].contains(params.celltypist_training_mode) ) {
    error "Invalid --celltypist_training_mode '${params.celltypist_training_mode}'. Expected one of: standard, detailed"
}
if( runScanvi && params.scanvi_hvg_n_top_genes.toInteger() < 0 ) {
    error "Invalid --scanvi_hvg_n_top_genes '${params.scanvi_hvg_n_top_genes}'. Expected a non-negative integer"
}
if( params.scanvi_model_dir && params.scvi_reference_model_dir ) {
    error "Use either --scanvi_model_dir for a pretrained scANVI model or --scvi_reference_model_dir to initialize scANVI training from a pretrained scVI reference model, not both"
}
if( params.marker_genes_n_top.toInteger() < 1 ) {
    error "Invalid --marker_genes_n_top ${params.marker_genes_n_top}. Expected a positive integer"
}
if( !["wilcoxon", "t-test"].contains(params.marker_genes_method) ) {
    error "Invalid --marker_genes_method ${params.marker_genes_method}. Expected one of: wilcoxon, t-test"
}
if( params.marker_dotplot_n_top.toInteger() < 1 ) {
    error "Invalid --marker_dotplot_n_top ${params.marker_dotplot_n_top}. Expected a positive integer"
}

// --- Model sheet -----------------------------------------------------------------------
// ANNOTATE_WITH_MODELS annotates queries from models that already exist, with one
// (model, query) pair per sheet row. The sheet is the only place those pairs are named, so
// it is parsed and fully validated here: a bad path, an unusable ID, or a pair that would
// overwrite another pair's outputs fails before any task is submitted.
def modelSheetColumns = ['tool', 'model', 'model_id', 'query_h5ad']
def modelSheetUnits = []
if( params.model_sheet ) {
    def sheetFile = file(params.model_sheet, checkIfExists: true)
    def sheetLines = sheetFile.readLines()
        .collect { it.trim() }
        .findAll { it && !it.startsWith('#') }
    if( sheetLines.size() < 2 ) {
        error "--model_sheet '${sheetFile}' needs a header row and at least one data row"
    }

    // A tab anywhere in the header means a TSV; anything else is read as a CSV. Fields are
    // split on the delimiter with no quoting rules, so a path cannot contain it.
    def delimiter = sheetLines[0].contains('\t') ? '\t' : ','
    def header = sheetLines[0].split(delimiter, -1).collect { it.trim().toLowerCase().replaceAll(/^"|"$/, '') }
    def unknownColumns = (header.findAll { it } as Set) - (modelSheetColumns as Set)
    if( unknownColumns ) {
        error "Unknown --model_sheet columns: ${unknownColumns.sort().join(', ')}. Expected any of: ${modelSheetColumns.join(', ')}"
    }
    ['tool', 'model', 'query_h5ad'].each { requiredColumn ->
        if( !header.contains(requiredColumn) ) {
            error "--model_sheet is missing the required '${requiredColumn}' column. Found: ${header.findAll { it }.join(', ')}"
        }
    }

    def seenPairs = [:]
    def sampleQueryPaths = [:]
    sheetLines.drop(1).eachWithIndex { line, index ->
        // Numbered over the data rows rather than the file's own lines, since blank and
        // commented lines were already dropped.
        def rowLabel = "row ${index + 1}"
        def values = line.split(delimiter, -1).collect { it.trim().replaceAll(/^"|"$/, '') }
        if( values.size() != header.size() ) {
            error "--model_sheet ${rowLabel} has ${values.size()} fields but the header declares ${header.size()}"
        }
        def row = [header, values].transpose().collectEntries { column, value -> [column, value] }

        def tool = row.tool?.toLowerCase()
        if( !knownTools.contains(tool) ) {
            error "--model_sheet ${rowLabel} names an unknown tool '${row.tool}'. Expected one of: ${knownTools.join(', ')}"
        }
        if( tool == 'panhumanpy' && row.model ) {
            error "--model_sheet ${rowLabel} gives a model for panhumanpy, which carries its own model. Leave the model column empty."
        }
        if( tool != 'panhumanpy' && !row.model ) {
            error "--model_sheet ${rowLabel} has no model path, which ${tool} requires"
        }

        def modelPath = row.model ? file(row.model, checkIfExists: true) : null
        if( tool == 'celltypist' && modelPath.isDirectory() ) {
            error "--model_sheet ${rowLabel} points celltypist at a directory. Give it the celltypist_model.pkl file."
        }
        if( tool == 'scanvi' && !modelPath.isDirectory() ) {
            error "--model_sheet ${rowLabel} points scanvi at a file. Give it the scanvi_model directory saved by this pipeline."
        }

        // Every model this pipeline trains is called celltypist_model.pkl or scanvi_model,
        // so the default ID collides as soon as two training runs are combined. The
        // uniqueness check below turns that collision into an instruction to set model_id.
        def modelId = row.model_id ?: (modelPath ? modelPath.baseName : tool)
        if( !sampleIdPattern.matcher(modelId).matches() ) {
            error "--model_sheet ${rowLabel} produces an unusable model ID '${modelId}'. Model IDs name output paths and must match ${sampleIdPattern}."
        }

        def resolvedQueries = file(row.query_h5ad, checkIfExists: true)
        resolvedQueries = resolvedQueries instanceof Collection ? resolvedQueries.toList() : [resolvedQueries]
        if( !resolvedQueries ) {
            error "--model_sheet ${rowLabel} matched no query h5ad files with '${row.query_h5ad}'"
        }

        resolvedQueries.sort { it.toString() }.each { queryFile ->
            def sampleId = queryFile.baseName
            if( !sampleIdPattern.matcher(sampleId).matches() ) {
                error "--model_sheet ${rowLabel} produces an unusable sample ID '${sampleId}' from '${queryFile.name}'. A sample ID is the filename without its final extension and must match ${sampleIdPattern}."
            }
            def previousPath = sampleQueryPaths.put(sampleId, queryFile.toString())
            if( previousPath && previousPath != queryFile.toString() ) {
                error "--model_sheet uses sample ID '${sampleId}' for two different files: ${previousPath} and ${queryFile}. Sample IDs come from the filename without its extension, so rename one of them."
            }
            def unitId = "${sampleId}.${modelId}".toString()
            def previousRow = seenPairs.put(unitId, rowLabel)
            if( previousRow ) {
                error "--model_sheet pairs sample '${sampleId}' with model ID '${modelId}' twice (${previousRow} and ${rowLabel}). Each pair gets its own output directory, so give the models distinct model_id values."
            }
            modelSheetUnits << [
                [id: unitId, sample_id: sampleId, tools: tool],
                queryFile,
                modelPath ?: [],
            ]
        }
    }

    if( !modelSheetUnits ) {
        error "--model_sheet '${sheetFile}' produced no (model, query) pairs"
    }
    // Sorted so channel order and -resume hashes do not depend on sheet row order.
    modelSheetUnits = modelSheetUnits.sort { it[0].id }
}

process VALIDATE_REFERENCE {
    label 'validation'
    publishDir outputDir, mode: 'copy', overwrite: true

    input:
    path ref_h5ad

    output:
    path 'reference_manifest.json', emit: manifest
    path 'reference_genes.txt', emit: genes

    script:
    """
    python3 ${projectDir}/scripts/validate_reference.py \
        --ref-h5ad ${ref_h5ad} \
        --ref-label-col '${params.ref_label_col}' \
        ${(runScanvi && !useNfScviMetricsConfig && params.batch_key) ? "--batch-key '${params.batch_key}'" : ''} \
        ${(runScanvi && !useNfScviMetricsConfig && params.scanvi_categorical_covariate_keys) ? "--categorical-covariate-keys '${params.scanvi_categorical_covariate_keys}'" : ''} \
        ${(runScanvi && !useNfScviMetricsConfig && params.scanvi_continuous_covariate_keys) ? "--continuous-covariate-keys '${params.scanvi_continuous_covariate_keys}'" : ''} \
        ${useNfScviMetricsConfig ? "--nf-scvi-metrics-config ${nfScviMetricsConfig}" : ''} \
        ${refFilterCol ? "--ref-filter-col '${refFilterCol}'" : ''} \
        ${refFilterValues ? "--ref-filter-values '${refFilterValues}'" : ''} \
        --output reference_manifest.json \
        --genes-output reference_genes.txt
    """
}

process VALIDATE_QUERY {
    label 'validation'
    tag "${meta.id}"
    publishDir "${outputDir}/${meta.id}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(query_h5ad)
    // Empty lists when no reference h5ad takes part in the run, which is what the falsy
    // checks in the command detect: an unstaged path input arrives as an empty list.
    path reference_genes
    path reference_manifest
    // Kept out of the meta map on purpose. Query validation is the only step whose
    // requirements depend on the merge step, so only its own hash should change when the
    // merge step is turned on or off, not the hash of every expensive annotation task.
    val merge_annotations

    output:
    tuple val(meta), path("${meta.id}.query_manifest.json"), emit: manifest

    script:
    def unitTools = meta.tools.tokenize(',')
    def scanviSetup = unitTools.contains('scanvi') && !useNfScviMetricsConfig
    """
    python3 ${projectDir}/scripts/validate_query.py \
        --query-h5ad ${query_h5ad} \
        --sample-id '${meta.sample_id}' \
        --tools '${meta.tools}' \
        ${reference_genes ? "--reference-genes ${reference_genes}" : ''} \
        ${reference_manifest ? "--reference-manifest ${reference_manifest}" : ''} \
        ${(scanviSetup && params.batch_key) ? "--batch-key '${params.batch_key}'" : ''} \
        ${(scanviSetup && params.scanvi_categorical_covariate_keys) ? "--categorical-covariate-keys '${params.scanvi_categorical_covariate_keys}'" : ''} \
        ${(scanviSetup && params.scanvi_continuous_covariate_keys) ? "--continuous-covariate-keys '${params.scanvi_continuous_covariate_keys}'" : ''} \
        ${(unitTools.contains('scanvi') && useNfScviMetricsConfig) ? "--nf-scvi-metrics-config ${nfScviMetricsConfig}" : ''} \
        ${(unitTools.contains('panhumanpy') && params.panhumanpy_feature_names_col) ? "--panhumanpy-feature-names-col '${params.panhumanpy_feature_names_col}'" : ''} \
        --query-cluster-col '${params.query_cluster_col}' \
        ${params.compute_missing_clusters ? '--compute-missing-clusters' : ''} \
        ${merge_annotations ? '' : '--skip-merge'} \
        --output ${meta.id}.query_manifest.json
    """
}

process CELLTYPIST_TRAIN_MODEL {
    label 'cpu4'
    publishDir outputDir, mode: 'copy', overwrite: true

    input:
    path ref_h5ad
    path reference_manifest

    output:
    path 'celltypist_model.pkl', emit: model
    path 'celltypist_model_metadata.json', emit: metadata

    script:
    """
    python3 ${projectDir}/scripts/train_celltypist.py \
        --ref-h5ad ${ref_h5ad} \
        --reference-manifest ${reference_manifest} \
        --ref-label-col '${params.ref_label_col}' \
        --output-model celltypist_model.pkl \
        --output-metadata celltypist_model_metadata.json \
        --n-jobs ${task.cpus} \
        --training-mode '${params.celltypist_training_mode}' \
        --balance-cells-per-label ${params.celltypist_balance_cells_per_label} \
        --feature-selection-top-genes ${params.celltypist_feature_selection_top_genes} \
        --feature-selection-max-iter ${params.celltypist_feature_selection_max_iter} \
        --max-iter ${params.celltypist_max_iter} \
        ${params.celltypist_feature_selection ? '--feature-selection' : ''} \
        ${params.celltypist_use_gpu ? '--use-gpu' : ''}
    """
}

process CELLTYPIST_ANNOTATE {
    label 'normal'
    tag "${meta.id}"
    // Published so that a run which skips the merge step still has deliverables outside the
    // work directory. Prediction tables are small next to the annotated object.
    publishDir "${outputDir}/${meta.id}", mode: 'copy', overwrite: true,
        saveAs: { filename -> "${params.project_tag}.${meta.id}.${filename}" }

    input:
    tuple val(meta), path(query_h5ad), path(celltypist_model)

    output:
    tuple val(meta), path('celltypist_predictions.csv'), emit: predictions

    script:
    """
    python3 ${projectDir}/scripts/annotate_celltypist.py \
        --query-h5ad ${query_h5ad} \
        --model-path ${celltypist_model} \
        --output celltypist_predictions.csv
    """
}

process SCANVI_TRAIN_MODEL {
    label 'gpu'
    publishDir outputDir, mode: 'copy', overwrite: true

    input:
    path ref_h5ad
    path reference_manifest

    output:
    path 'scvi_model', emit: scvi_model_dir
    path 'scvi_model_metadata.json', emit: scvi_metadata
    path 'scanvi_model', emit: model_dir
    path 'scanvi_model_metadata.json', emit: metadata

    script:
    """
    python3 ${projectDir}/scripts/train_scanvi.py \
        --ref-h5ad ${ref_h5ad} \
        --reference-manifest ${reference_manifest} \
        --ref-label-col '${params.ref_label_col}' \
        ${(!useNfScviMetricsConfig && params.batch_key) ? "--batch-key '${params.batch_key}'" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_categorical_covariate_keys) ? "--categorical-covariate-keys '${params.scanvi_categorical_covariate_keys}'" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_continuous_covariate_keys) ? "--continuous-covariate-keys '${params.scanvi_continuous_covariate_keys}'" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_layer) ? "--layer '${params.scvi_layer}'" : ''} \
        --unknown-label '${params.unknown_label}' \
        --hvg-n-top-genes ${params.scanvi_hvg_n_top_genes} \
        ${useNfScviMetricsConfig ? "--training-config ${nfScviMetricsConfig}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_n_hidden != null) ? "--scvi-n-hidden ${params.scvi_n_hidden}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_n_latent != null) ? "--scvi-n-latent ${params.scvi_n_latent}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_n_layers != null) ? "--scvi-n-layers ${params.scvi_n_layers}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_dropout_rate != null) ? "--scvi-dropout-rate ${params.scvi_dropout_rate}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_dispersion) ? "--scvi-dispersion '${params.scvi_dispersion}'" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_gene_likelihood) ? "--scvi-gene-likelihood '${params.scvi_gene_likelihood}'" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_latent_distribution) ? "--scvi-latent-distribution '${params.scvi_latent_distribution}'" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_train_max_epochs != null) ? "--scvi-train-max-epochs ${params.scvi_train_max_epochs}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_train_batch_size != null) ? "--scvi-train-batch-size ${params.scvi_train_batch_size}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_train_train_size != null) ? "--scvi-train-train-size ${params.scvi_train_train_size}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_train_validation_size != null) ? "--scvi-train-validation-size ${params.scvi_train_validation_size}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_train_early_stopping != null) ? "--scvi-train-early-stopping ${params.scvi_train_early_stopping}" : ''} \
        ${(!useNfScviMetricsConfig && params.scvi_train_early_stopping_patience != null) ? "--scvi-train-early-stopping-patience ${params.scvi_train_early_stopping_patience}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_train_max_epochs != null) ? "--scanvi-train-max-epochs ${params.scanvi_train_max_epochs}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_train_batch_size != null) ? "--scanvi-train-batch-size ${params.scanvi_train_batch_size}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_train_train_size != null) ? "--scanvi-train-train-size ${params.scanvi_train_train_size}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_train_validation_size != null) ? "--scanvi-train-validation-size ${params.scanvi_train_validation_size}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_train_early_stopping != null) ? "--scanvi-train-early-stopping ${params.scanvi_train_early_stopping}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_train_early_stopping_patience != null) ? "--scanvi-train-early-stopping-patience ${params.scanvi_train_early_stopping_patience}" : ''} \
        --output-scvi-model-dir scvi_model \
        --output-scvi-metadata scvi_model_metadata.json \
        --output-model-dir scanvi_model \
        --output-metadata scanvi_model_metadata.json \
        ${params.scvi_reference_model_dir ? "--scvi-model-dir ${scviReferenceModelDir}" : ''}
    """
}

process SCANVI_ANNOTATE {
    label 'gpu'
    tag "${meta.id}"
    publishDir "${outputDir}/${meta.id}", mode: 'copy', overwrite: true,
        saveAs: { filename -> "${params.project_tag}.${meta.id}.${filename}" }

    input:
    tuple val(meta), path(query_h5ad), path(scanvi_model_dir)

    output:
    tuple val(meta), path('scanvi_predictions.csv'), emit: predictions
    tuple val(meta), path('scanvi_latent.h5ad'), emit: latent

    script:
    """
    python3 ${projectDir}/scripts/annotate_scanvi.py \
        --query-h5ad ${query_h5ad} \
        --model-dir ${scanvi_model_dir} \
        --label-column '${params.scanvi_obs_col}' \
        --score-column '${params.scanvi_score_col}' \
        ${useNfScviMetricsConfig ? "--training-config ${nfScviMetricsConfig}" : "--query-max-epochs ${params.scanvi_query_max_epochs}"} \
        ${(!useNfScviMetricsConfig && params.scanvi_query_batch_size != null) ? "--query-batch-size ${params.scanvi_query_batch_size}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_query_train_size != null) ? "--query-train-size ${params.scanvi_query_train_size}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_query_validation_size != null) ? "--query-validation-size ${params.scanvi_query_validation_size}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_query_early_stopping != null) ? "--query-early-stopping ${params.scanvi_query_early_stopping}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_query_early_stopping_patience != null) ? "--query-early-stopping-patience ${params.scanvi_query_early_stopping_patience}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_query_check_val_every_n_epoch != null) ? "--query-check-val-every-n-epoch ${params.scanvi_query_check_val_every_n_epoch}" : ''} \
        ${(!useNfScviMetricsConfig && params.scanvi_query_plan_weight_decay != null) ? "--query-plan-weight-decay ${params.scanvi_query_plan_weight_decay}" : ''} \
        --output scanvi_predictions.csv \
        --latent-output scanvi_latent.h5ad
    """
}

process PANHUMANPY_ANNOTATE {
    label 'normal'
    tag "${meta.id}"
    publishDir "${outputDir}/${meta.id}", mode: 'copy', overwrite: true,
        saveAs: { filename -> "${params.project_tag}.${meta.id}.${filename}" }

    input:
    tuple val(meta), path(query_h5ad)

    output:
    tuple val(meta), path('panhumanpy_predictions.csv'), emit: predictions

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
    tag "${meta.id}"
    publishDir "${outputDir}/${meta.id}", mode: 'copy', overwrite: true

    input:
    // Per-tool prediction files, staged under the fixed names declared by the
    // annotation processes. Only the tools this unit ran contribute entries.
    tuple val(meta), path(query_h5ad), path(annotation_files)

    output:
    tuple val(meta), path("${params.project_tag}.${meta.id}.annotated.h5ad"), emit: h5ad
    path "${params.project_tag}.${meta.id}.cluster_annotation_summary.csv", emit: cluster_summary
    path "${params.project_tag}.${meta.id}.cluster_marker_genes.tsv", emit: marker_genes
    path "${params.project_tag}.${meta.id}.cluster_marker_dotplot.png", emit: marker_dotplot, optional: true

    script:
    def unitTools = meta.tools.tokenize(',')
    def prefix = "${params.project_tag}.${meta.id}"
    """
    python3 ${projectDir}/scripts/merge_annotations.py \
        --query-h5ad ${query_h5ad} \
        ${unitTools.contains('celltypist') ? '--celltypist-csv celltypist_predictions.csv' : ''} \
        ${unitTools.contains('scanvi') ? '--scanvi-csv scanvi_predictions.csv --scanvi-latent-h5ad scanvi_latent.h5ad' : ''} \
        ${unitTools.contains('panhumanpy') ? '--panhumanpy-csv panhumanpy_predictions.csv' : ''} \
        --cluster-column '${params.query_cluster_col}' \
        --scanvi-label-column '${params.scanvi_obs_col}' \
        --scanvi-score-column '${params.scanvi_score_col}' \
        ${params.compute_missing_clusters ? '--compute-missing-clusters' : ''} \
        --output ${prefix}.annotated.h5ad \
        --cluster-summary-output ${prefix}.cluster_annotation_summary.csv \
        --marker-genes-output ${prefix}.cluster_marker_genes.tsv \
        --marker-genes-n-top ${params.marker_genes_n_top} \
        --marker-genes-method ${params.marker_genes_method} \
        ${params.marker_dotplot ? "--marker-dotplot-output ${prefix}.cluster_marker_dotplot.png --marker-dotplot-n-top ${params.marker_dotplot_n_top}" : ""}
    """
}

workflow {
    if( params.model_sheet ) {
        error "--model_sheet is only read by the ANNOTATE_WITH_MODELS entry point. Run it with -entry ANNOTATE_WITH_MODELS, or drop --model_sheet to train and annotate in one run."
    }
    if( !querySamples ) {
        error "Missing required parameter: --query_h5ad"
    }
    if( needsReference && !params.ref_h5ad ) {
        error "Missing required parameter: --ref_h5ad. It is required when a reference-based tool is trained (celltypist or scanvi without a pretrained model)"
    }
    if( useReference && !params.ref_label_col ) {
        error "Missing required parameter: --ref_label_col. It is required whenever a reference h5ad is used"
    }
    if( !needsReference && (refFilterCol || refFilterValues) ) {
        error "Reference filtering was requested but no reference-based tool is being trained; filtering only affects model training"
    }

    // Merge is the default here: the annotated h5ad, cluster summary, and marker genes are
    // what this workflow exists to produce.
    def mergeAnnotations = params.merge_annotations == null ? true : params.merge_annotations.toString().toBoolean()
    log.info "Merge step: ${mergeAnnotations ? 'enabled' : 'skipped, publishing per-tool predictions only'}"

    // An empty list stages no reference file, which is what a query-only run
    // (for example --tools panhumanpy, or annotation from pretrained models) needs.
    Channel
        .value(useReference ? refH5ad : [])
        .set { ref_h5ad_ch }

    Channel
        .fromList(querySamples)
        .set { query_ch }

    // Every reference-side process takes only value channels as input, so its outputs are
    // value channels too and can be read once per query. Keep it that way: if one of these
    // ever became a queue channel it would pair off against the first query and silently
    // leave the rest unannotated.
    if( useReference ) {
        reference_outputs = VALIDATE_REFERENCE(ref_h5ad_ch)
        reference_manifest_ch = reference_outputs.manifest
        reference_genes_ch = reference_outputs.genes
    }
    else {
        reference_manifest_ch = Channel.value([])
        reference_genes_ch = Channel.value([])
    }

    // Joining the per-sample manifest back onto the query makes annotation wait for that
    // sample's validation, while still letting samples progress independently.
    validated_query_ch = query_ch
        .join(VALIDATE_QUERY(
            query_ch,
            reference_genes_ch,
            reference_manifest_ch,
            Channel.value(mergeAnnotations),
        ).manifest)
        .map { meta, query_h5ad, query_manifest -> [meta, query_h5ad] }

    annotation_channels = []

    if( runCelltypist ) {
        celltypist_model_ch = celltypistMode == 'pretrained'
            ? Channel.value(file(params.celltypist_model, checkIfExists: true))
            : CELLTYPIST_TRAIN_MODEL(ref_h5ad_ch, reference_manifest_ch).model

        // Combining against a value channel hands every query the same model without
        // consuming it, which is what lets one training run serve all of the queries.
        annotation_channels << CELLTYPIST_ANNOTATE(validated_query_ch.combine(celltypist_model_ch)).predictions
    }

    if( runScanvi ) {
        scanvi_model_ch = scanviMode == 'pretrained'
            ? Channel.value(file(params.scanvi_model_dir, checkIfExists: true))
            : SCANVI_TRAIN_MODEL(ref_h5ad_ch, reference_manifest_ch).model_dir

        scanvi_annotation_outputs = SCANVI_ANNOTATE(validated_query_ch.combine(scanvi_model_ch))
        annotation_channels << scanvi_annotation_outputs.predictions
        annotation_channels << scanvi_annotation_outputs.latent
    }

    if( runPanhumanpy ) {
        annotation_channels << PANHUMANPY_ANNOTATE(validated_query_ch).predictions
    }

    if( mergeAnnotations ) {
        // Each selected tool contributes a known number of files per sample, so grouping by
        // that size releases a sample as soon as its own annotations land instead of waiting
        // for every other sample. Sorted for a stable staging order, which keeps -resume warm.
        def expectedAnnotationFiles = (runCelltypist ? 1 : 0) + (runScanvi ? 2 : 0) + (runPanhumanpy ? 1 : 0)
        annotation_files_ch = annotation_channels
            .tail()
            .inject(annotation_channels.head()) { merged, channel -> merged.mix(channel) }
            .groupTuple(size: expectedAnnotationFiles, sort: true)

        // failOnMismatch turns a sample whose annotation group never completed into a hard
        // error. Without it groupTuple simply never emits that key and the sample would
        // vanish from the results while the run still reported success.
        MERGE_ANNOTATIONS(validated_query_ch.join(annotation_files_ch, failOnMismatch: true, failOnDuplicate: true))
    }
}

workflow TRAIN_MODELS {
    if( params.model_sheet ) {
        error "--model_sheet is only read by the ANNOTATE_WITH_MODELS entry point, which annotates from models that already exist"
    }
    if( !runCelltypist && !runScanvi ) {
        error "The TRAIN_MODELS entry point needs --tools to include celltypist or scanvi"
    }
    if( !trainsCelltypist && !trainsScanvi ) {
        error "Nothing to train: every selected reference-based tool was given a pretrained model"
    }
    if( !params.ref_h5ad ) {
        error "Missing required parameter: --ref_h5ad. It is required when a reference-based tool is trained (celltypist or scanvi without a pretrained model)"
    }
    if( !params.ref_label_col ) {
        error "Missing required parameter: --ref_label_col. It is required whenever a reference h5ad is used"
    }

    Channel
        .value(refH5ad)
        .set { ref_h5ad_ch }

    // Training is independent of the query objects, so this entry point never reads them.
    reference_manifest_ch = VALIDATE_REFERENCE(ref_h5ad_ch).manifest

    if( trainsCelltypist ) {
        CELLTYPIST_TRAIN_MODEL(ref_h5ad_ch, reference_manifest_ch)
    }
    if( trainsScanvi ) {
        SCANVI_TRAIN_MODEL(ref_h5ad_ch, reference_manifest_ch)
    }
}

workflow ANNOTATE_WITH_MODELS {
    if( !params.model_sheet ) {
        error "The ANNOTATE_WITH_MODELS entry point needs --model_sheet: a CSV or TSV naming one tool, model, and query per row"
    }
    // Every tool, model, and query comes from the sheet here, so reject the parameters that
    // describe the same things globally instead of quietly ignoring them.
    def sheetOnlyParams = [
        '--query_h5ad': params.query_h5ad,
        '--ref_h5ad': params.ref_h5ad,
        '--ref_label_col': params.ref_label_col,
        '--ref_filter_col': params.ref_filter_col,
        '--ref_filter_values': params.ref_filter_values,
        '--celltypist_model': params.celltypist_model,
        '--scanvi_model_dir': params.scanvi_model_dir,
        '--scvi_reference_model_dir': params.scvi_reference_model_dir,
    ].findAll { name, value -> value as boolean }.keySet().sort()
    if( sheetOnlyParams ) {
        error "ANNOTATE_WITH_MODELS takes its tools, models, and queries from --model_sheet, so ${sheetOnlyParams.join(', ')} would be ignored. Remove ${sheetOnlyParams.size() > 1 ? 'them' : 'it'} and describe the run in the sheet instead."
    }

    // Merge is opt-in here. Each pair runs a single tool, so the usual reason to merge -
    // collecting several tools into one object - does not apply, and the merge step is by
    // far the most expensive part of a run.
    def mergeAnnotations = params.merge_annotations == null ? false : params.merge_annotations.toString().toBoolean()
    log.info "Annotating ${modelSheetUnits.size()} (model, query) pair(s) from ${params.model_sheet}"
    log.info "Merge step: ${mergeAnnotations ? 'enabled' : 'skipped, publishing per-tool predictions only'}"

    // One item per (model, query) pair. The model is an empty list for panhumanpy, which
    // brings its own model and takes no model input.
    Channel
        .fromList(modelSheetUnits)
        .set { unit_ch }

    // No reference takes part: every model already carries the gene space it was trained on,
    // and each annotation step reports its own overlap against that space.
    validated_unit_ch = unit_ch
        .join(VALIDATE_QUERY(
            unit_ch.map { meta, query_h5ad, model -> [meta, query_h5ad] },
            Channel.value([]),
            Channel.value([]),
            Channel.value(mergeAnnotations),
        ).manifest)
        .map { meta, query_h5ad, model, query_manifest -> [meta, query_h5ad, model] }

    unit_by_tool = validated_unit_ch.branch { meta, query_h5ad, model ->
        celltypist: meta.tools == 'celltypist'
        scanvi: meta.tools == 'scanvi'
        panhumanpy: meta.tools == 'panhumanpy'
    }

    celltypist_outputs = CELLTYPIST_ANNOTATE(unit_by_tool.celltypist)
    scanvi_outputs = SCANVI_ANNOTATE(unit_by_tool.scanvi)
    panhumanpy_outputs = PANHUMANPY_ANNOTATE(
        unit_by_tool.panhumanpy.map { meta, query_h5ad, model -> [meta, query_h5ad] }
    )

    if( mergeAnnotations ) {
        // One tool per pair, so there is nothing to group across processes: a pair's
        // annotation files are complete as soon as its own annotation step finishes.
        annotation_files_ch = celltypist_outputs.predictions
            .map { meta, predictions -> [meta, [predictions]] }
            .mix(
                scanvi_outputs.predictions
                    .join(scanvi_outputs.latent)
                    .map { meta, predictions, latent -> [meta, [predictions, latent]] }
            )
            .mix(panhumanpy_outputs.predictions.map { meta, predictions -> [meta, [predictions]] })

        MERGE_ANNOTATIONS(
            validated_unit_ch
                .map { meta, query_h5ad, model -> [meta, query_h5ad] }
                .join(annotation_files_ch, failOnMismatch: true, failOnDuplicate: true)
        )
    }
}
