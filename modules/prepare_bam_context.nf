#!/usr/bin/env nextflow

process PREPARE_BAM_CONTEXT {
    label 'process_medium'

    input:
    val bam_bai_collection

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), path("*.control.bam"), path("*.control.bai"), emit: bam_context

    script:
    def sample_id = meta.sample_id
    
    """
    samtools \\
        flagstat \\
        --threads ${task.cpus} \\
        $bam \\
        > ${sample_id}.flagstat
    """
}
