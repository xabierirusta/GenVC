// Post_alignment processing for WGS
process wgs_post_alignment {
    publishDir "${params.outDir}/post_alignment/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(sorted_bam)                         // Sorted BAM files input
    path ref_files                                          // Reference genome related files                             
    path ref_genome                                         // Reference genome fasta

    output:
    tuple val(ID), path("${ID}_mdup.bam"), emit: mdup_bam_bai // Markduplicates BAM file
    path "${ID}_mdup_metrics.txt"                             // Picard Markduplicates metrics text file

    script:
    """
    # Run Picard MarkDuplicates
    java -jar /usr/picard/picard.jar MarkDuplicates \
        I="${ID}_RG.bam" \
        O="${ID}_mdup.bam" \
        M="${ID}_mdup_metrics.txt"

    rm ${sorted_bam}
    """
}
