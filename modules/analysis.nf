#!/usr/bin/env nextflow

process SPLIT_MALIGNANT {
    input:
        path qc_data
    output:
        tuple(path("${params.outdir}/nonmalignant/seurat/*"), path("${params.outdir}/malignant/seurat/*")), emit: (nonmal_ch, mal_ch)
    script:
        """
        Rscript scripts/03-split_malignant_nonmalignant.R --input ${qc_data} --outdir ${params.outdir}
        """
}

process NONMALIGNANT_DIMREDUCE {
    input:
        path nonmal_ch
    output:
        path "${params.outdir}/nonmalignant/dimreduce/*", emit: nonmal_res
    script:
        """
        Rscript scripts/04-nonmalignant_dimensionality_reduction.R --input ${nonmal_ch} --outdir ${params.outdir}/nonmalignant/dimreduce
        """
}

process NONMALIGNANT_CLUSTER {
    input:
        path nonmal_res
    output:
        path "${params.outdir}/nonmalignant/cluster/*"
    script:
        """
        Rscript scripts/05-nonmalignant_clustering.R --input ${nonmal_res} --outdir ${params.outdir}/nonmalignant/cluster
        """
}

process MALIGNANT_DIMREDUCE {
    input:
        path mal_ch
    output:
        path "${params.outdir}/malignant/dimreduce/*", emit: mal_res
    script:
        """
        Rscript scripts/06-malignant_dimensionality_reduction.R --input ${mal_ch} --outdir ${params.outdir}/malignant/dimreduce
        """
}

process TRAJECTORY {
    input:
        path mal_res
    output:
        path "${params.outdir}/malignant/trajectory/*"
    script:
        """
        Rscript scripts/07-malignant_trajectory_analysis.R --input ${mal_res} --outdir ${params.outdir}/malignant/trajectory
        """
}
