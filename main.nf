#!/usr/bin/env nextflow

process CELLTYPIST_CREATEREFERENCE {   
    input:
    path(ref_h5ad)
    
    output:
    path(ctp_model)
    
    script:
    """
    celltypist ....
    """
}

process CELLTYPIST_ANNOTATE {
    input:
    tuple path(ref_h5ad), path(query_h5ad), path(ctp_model)

    output:
    path(ctp_csv)

    script:
    """
    python  ...
    """
}

process SCANVI_ANNOTATE {
    input:
    tuple path(ref_h5ad), path(query_h5ad)

    output:
    path(scv_csv)

    script:
    """
    python ...
    """
}

process MERGE_ANNOTATIONS {
    input:
    tuple path(query_h5ad), path(ctp_csv), path(scv_csv)

    output:
    path(final_h5ad)

    script:
    """
    python ...
    """
}

workflow {

    CELLTYPIST_CREATEREFERENCE(params.ref_h5ad)
    CELLTYPIST_ANNOTATE(params.ref_h5ad, params.query_h5ad, CELLTYPIST_CREATEREFERENCE.out)

    SCANVI_ANNOTATE(params.ref_h5ad, params.query_h5ad)
    MERGE_ANNOTATIONS(params.query_h5ad, CELLTYPIST_ANNOTATE.out, SCANVI_ANNOTATE.out)

}
