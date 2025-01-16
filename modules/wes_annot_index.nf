// Indexing of annotated vcf files
process wes_annot_index {
    publishDir "${params.outDir}/post_alignment/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(mdup_bam)                   // Markduplicates processed BAM file

    output:
    tuple val(ID), path("${ID}_mdup.bam"), path("${ID}_mdup.bam.bai"), emit: mdup_bam_bai   // Output channel contianing the markduplicates BAM and BAI files
    path "${ID}_stats.txt"                                                                  // Output samtools stats file
    script:
    """
    # Run samtools stats
    samtools stats ${mdup_bam} > "${ID}_stats.txt"
    samtools index ${mdup_bam}
    """
}
