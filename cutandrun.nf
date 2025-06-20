#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { BOWTIE2_BUILD } from './modules/bowtie2_build.nf'
include { BOWTIE2_ALIGN } from './modules/bowtie2_align.nf'
include { BAM_FILTER } from './modules/bam_filter.nf'

include { PICARD_MARKDUPLICATES } from './modules/picard_markduplicates.nf'
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf'
include { SAMTOOLS_IDXSTATS } from './modules/samtools_idxstats.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools_flagstat.nf'

// include { MACS2_CALLPEAK } from './modules/macs2_callpeak.nf'

include { MULTIQC } from './modules/multiqc.nf'

/*
 * Pipeline parameters
 */

// Primary input
params.input_csv = "data/samplesheet.csv"
params.blacklist = "data/hg38-blacklist.bed"
params.fasta = "test-datasets/cutandrun/reference/genomes/hg38-chr20.fa.gz"
params.gtf = "test-datasets/cutandrun/reference/genomes/hg38-chr20-genes.gtf.gz"
params.bowtie2 = "test-datasets/cutandrun/reference/genomes/hg38-chr20-bowtie2.tar.gz"
params.gene_bed = "test-datasets/cutandrun/reference/genomes/hg38-chr20-genes.bed.gz"
params.report_id = "cutandrun"

log.info """\
==================================================================================================
\033[36m        ____ ____  __  __    ____ _   _ _____ ___   ____  _   _ _   _   _____           _ 
       / ___| __ )|  \\/  |  / ___| | | |_   _( _ ) |  _ \\| | | | \\ | | |_   _|__   ___ | |
      | |  _|  _ \\| |\\/| | | |   | | | | | | / _ \\/\\ |_) | | | |  \\| |   | |/ _ \\ / _ \\| |
      | |_| | |_) | |  | | | |___| |_| | | || (_>  <  _ <| |_| | |\\  |   | | (_) | (_) | |
       \\____|____/|_|  |_|  \\____|\\___/  |_| \\___/\\/_| \\_\\\\___/|_| \\_|   |_|\\___/ \\___/|_|\033[0m
                                                                                     
==================================================================================================

  \033[1;37m▸ Key Parameters:\033[0m

    ▹ \033[33mExecuted By\033[0m       : ${System.getProperty('user.name')}
    ▹ \033[33mRun Date\033[0m          : ${workflow.start.format("yyyy-MM-dd HH:mm 'UTC'")}

                                                                Author: Bo Wang | Version: Beta
==================================================================================================
"""
.stripIndent(true)


process SAVE_MATCHED_SAMPLES {
    publishDir "results/matched_samples", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(control_bam), path(control_bai)

    output:
    path "${meta.sample_id}_matched_samples.tsv"

    script:
    """
    echo -e "sample_id\tdonor_id\tbam\tbai\tcontrol_bam\tcontrol_bai" > ${meta.sample_id}_matched_samples.tsv
    echo -e "${meta.sample_id}\t${meta.donor_id}\t${bam}\t${bai}\t${control_bam}\t${control_bai}" >> ${meta.sample_id}_matched_samples.tsv
    """
}


workflow {

    // INPUT CHANNEL CREATION
    read_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row ->
            def meta = [
                group    : row.group,
                sample_id: row.sample_id,
                donor_id : row.donor_id,
                control  : row.control
            ]
            tuple(meta, file(row.fastq_1), file(row.fastq_2))
        }

    // CORE PREPROCESSING AND ALIGNMENT
    FASTQC(read_ch)

    TRIM_GALORE(read_ch)

    BOWTIE2_BUILD(file(params.fasta))

    BOWTIE2_ALIGN(TRIM_GALORE.out.trimmed_reads, BOWTIE2_BUILD.out.bt2_index)

    BAM_FILTER(BOWTIE2_ALIGN.out.bam)

    // DUPLICATE REMOVAL AND QC STATISTICS

    PICARD_MARKDUPLICATES(BAM_FILTER.out.bam, file(params.fasta))
    
    SAMTOOLS_STATS(PICARD_MARKDUPLICATES.out.bam_bai, file(params.fasta))

    SAMTOOLS_IDXSTATS(PICARD_MARKDUPLICATES.out.bam_bai)

    SAMTOOLS_FLAGSTAT(PICARD_MARKDUPLICATES.out.bam_bai)

    // DOWNSTREAM ANALYSIS

    def exp_ch = PICARD_MARKDUPLICATES.out.bam_bai
        .filter { meta, bam, bai -> meta.control.toString() == '0' }
        .map { meta, bam, bai -> 
            tuple(meta.donor_id, meta, bam, bai)
        }

    def ctrl_ch = PICARD_MARKDUPLICATES.out.bam_bai
        .filter { meta, bam, bai -> meta.control.toString() == '1' }
        .map { meta, bam, bai -> 
            tuple(meta.donor_id, bam, bai)
        }

    def bam_paired_ch = exp_ch.combine(ctrl_ch, by: 0)
        .map { donor_id, exp_meta, exp_bam, exp_bai, ctrl_bam, ctrl_bai ->
            tuple(exp_meta, exp_bam, exp_bai, ctrl_bam, ctrl_bai)
        }

    bam_paired_ch.view { meta, exp_bam, exp_bai, ctrl_bam, ctrl_bai ->
        """
        [Matched Pair]
        Group     : ${meta.group}
        Sample ID : ${meta.sample_id}
        Donor ID  : ${meta.donor_id}
        Exp BAM   : ${exp_bam}
        Ctrl BAM  : ${ctrl_bam}
        Exp BAI   : ${exp_bai}
        Ctrl BAI  : ${ctrl_bai}
        """
    }

    // SUMMARY REPORT GENERATION

    MULTIQC(
        FASTQC.out.zip.mix(
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            BOWTIE2_ALIGN.out.log,
            PICARD_MARKDUPLICATES.out.metrics,
            SAMTOOLS_STATS.out.stats,
            SAMTOOLS_IDXSTATS.out.idxstats,
            SAMTOOLS_FLAGSTAT.out.flagstat
        ).collect(),
        params.report_id
    )



}