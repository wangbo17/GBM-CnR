#!/usr/bin/env nextflow

process BAM_FILTER_INDEX {
    label 'process_medium'

    container = 'oras://community.wave.seqera.io/library/samtools:1.22--105e5e643c53f059'
    publishDir "results/samtools/bam_filter_index", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.filtered.dedup.bam"), path("*.filtered.dedup.bam.bai"), emit: bam_bai

    script:
    def group = meta.group
    def sample_id = meta.sample_id

    """
    samtools view -f 3 -F 268 -q 20 -@ ${task.cpus} -h -b ${bam} > ${sample_id}.filtered.dedup.bam
    # -f 3     : include only reads that are properly paired
    # -F 268   : exclude reads that are unmapped, have an unmapped mate, or are secondary alignments
    # -q 20    : require mapping quality ≥ 20 to exclude ambiguous alignments
    # -h       : retain the BAM header in the output

    samtools index ${sample_id}.filtered.dedup.bam
    """
}
