#!/usr/bin/env nextflow

process PICARD_MARKDUPLICATES {
    label 'process_medium'
    
    publishDir "results/picard/markduplicates", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("*.filtered.dedup.bam"), path("*.filtered.dedup.bam.bai"), emit: bam_bai
    path("*.MarkDuplicates.metrics.txt"), emit: metrics

    script:
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072

    def sample_id = meta.sample_id
    
    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${sample_id}.filtered.dedup.bam \\
        --REFERENCE_SEQUENCE ${fasta} \\
        --METRICS_FILE ${sample_id}.MarkDuplicates.metrics.txt \\
        --ASSUME_SORT_ORDER coordinate \\
        --REMOVE_DUPLICATES true

    samtools index ${sample_id}.filtered.dedup.bam
    """
}