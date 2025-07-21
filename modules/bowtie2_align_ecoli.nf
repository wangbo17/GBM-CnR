#!/usr/bin/env nextflow

process BOWTIE2_ALIGN_ECOLI {
    label 'process_medium'

    container = 'oras://community.wave.seqera.io/library/bowtie2_samtools:6df3a3213a70e258'
    publishDir "results/bowtie2/align/ecoli", mode: 'copy'

    input:
    tuple val(meta), path(read1), path(read2)
    path bt2_index

    output:
    tuple val(meta), path("*.unmapped.1.gz"), path("*.unmapped.2.gz"), emit: 'unmapped_reads'
    tuple val(meta), path("*.sam"), emit: 'sam'
    path "*.log", emit: 'log'
    path "*.seqDepth.txt", emit: 'seqdepth'

    script:
    def sample_id = meta.sample_id

    """
    prefix_name=\$(basename \$(ls ${bt2_index}/*.1.bt2) | sed -E 's/\\.1\\.bt2\$//')
    prefix_path=${bt2_index}/\$prefix_name

    bowtie2 \\
        --end-to-end --very-sensitive \\
        --no-mixed --no-discordant \\
        --un-conc-gz ${sample_id}.unmapped.gz \\
        --phred33 -I 10 -X 700 \\
        -x \$prefix_path \\
        -1 ${read1} -2 ${read2} \\
        --threads ${task.cpus} \\
        -S ${sample_id}.bowtie2.spikein.sam \\
        2> ${sample_id}.bowtie2.spikein.log
    
    seqDepthDouble=\$(samtools view -F 0x04 ${sample_id}.bowtie2.spikein.sam | wc -l)
    seqDepth=\$((seqDepthDouble / 2))
    echo \$seqDepth > ${sample_id}.seqDepth.txt
    """
}