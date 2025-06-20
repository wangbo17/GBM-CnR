#!/usr/bin/env nextflow

process SAMTOOLS_STATS {
    label 'process_medium'
    
    publishDir "results/samtools/stats", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)

    output:
    path("*.stats"), emit: 'stats'

    script:
    def sample_id = meta.sample_id
    
    """
    samtools \\
        stats \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        ${bam} \\
        > ${sample_id}.stats
    """
}
