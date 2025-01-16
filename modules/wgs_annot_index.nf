// Indexing of annotated vcf files
process wgs_annot_index {
    publishDir "${params.outDir}/post_alignment/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(mdup_bam)            // Input markduplicates BAM file

    output:
    tuple val(ID), path("${ID}_mdup.bam"), path("${ID}_mdup.bam.bai"), emit: mdup_bam_bai  // Output channel containing BAM and BAI files
    path "${ID}_stats.txt"                                                                 // Samtools stats output txt file

    script:
    """
    # Run samtools stats
    samtools stats ${mdup_bam} > "${ID}_stats.txt"
    samtools index ${mdup_bam}
    """
}
