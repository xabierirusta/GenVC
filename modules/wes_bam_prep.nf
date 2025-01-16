// Bam_processing for WES 
process wes_bam_prep {
    publishDir "${params.outDir}/alignment/${ID}", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(ID), path(sam_file)                             // Aligned SAM file input with its ID

    output:
    tuple val(ID), path("${ID}.sorted.bam"), emit: sorted_bam // Output channel with sorted BAM files

    script:
    """
    # .sam to .bam
    samtools view -S -b ${sam_file} > "${ID}.bam"

    # Sort .bam file
    samtools sort -o "${ID}.sorted.bam" "${ID}.bam"

    rm ${sam_file}
    rm "${ID}.bam"
    """
}
