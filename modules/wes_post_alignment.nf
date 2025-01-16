// Post_alignment processing for WES
process wes_post_alignment {
    publishDir "${params.outDir}/post_alignment/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(bam_file), path(bai_file) // Filtrated BAM file with its BAI
    path ref_files                                // Reference genome files
    path ref_genome                               // Referenge genome fasta file
    path bed_files                                // Targed sequence BED input
    path targets_bed                              // BED file that contains interesting regions
    path ref_dict                                 // Reference genome dict file required by picard

    output:
    tuple val(ID), path("${ID}_mdup.bam"), emit: mdup_bam_bai // Markduplicates BAM file
    path "${ID}_mdup_metrics.txt"                             // Markduplicates metrics file
    path "${ID}_hs_metrics.txt"                               // Picard CollectHsMetric metrics file

    script:
    """
    # Convert the BED file into a interval_list for picard
    java -jar /usr/picard/picard.jar BedToIntervalList \
    I=${targets_bed} \
    O="${ID}_targets.interval_list" \
    SD=${ref_dict}

    # Run Picard MarkDuplicates
    java -jar /usr/picard/picard.jar MarkDuplicates \
        I=${bam_file} \
        O="${ID}_mdup.bam" \
        M="${ID}_mdup_metrics.txt"

    # Run Picard CollectHsMetrics
    java -jar /usr/picard/picard.jar CollectHsMetrics \
        I="${ID}_mdup.bam" \
        O="${ID}_hs_metrics.txt" \
        R=${ref_genome} \
        BAIT_INTERVALS="${ID}_targets.interval_list" \
        TARGET_INTERVALS="${ID}_targets.interval_list"

    rm ${bam_file}
    """
}
