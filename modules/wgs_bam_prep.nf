// Bam_processing for WGS
process wgs_bam_prep {
    publishDir "${params.outDir}/alignment/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(sam_file)     // Aligned SAM file input with its ID
    path ref_dict                     // Reference genome dict file

    output:
    tuple val(ID), path("${ID}_RG_filt.bam"), emit: bam_bai // Output sorted BAM files with BAI

    script:
    """
    # .sam to .bam
    samtools view -S -b ${sam_file} > "${ID}.bam"

    # Sort .bam file
    samtools sort -o "${ID}.sorted.bam" "${ID}.bam"

    # Extract chromosome names from the dictionary, add 'chr' prefix
    chromosomes=\$(grep "^@SQ" ${ref_dict} | cut -f2 | sed 's/SN://g' | sed 's/^chr//')

    # Modify the header of the filtered BAM file (replace SN:1 with SN:chr1, SN:2 with SN:chr2, etc.)
    samtools view -H "${ID}.sorted.bam" > "${ID}_header.sam"
    for i in \$chromosomes; do
        sed -i "s/SN:\$i/SN:chr\$i/g" "${ID}_header.sam"
    done

    # Reapply the modified header to the reheadered BAM file
    samtools reheader "${ID}_header.sam" "${ID}.sorted.bam" > "${ID}_RG_filt.bam"

    rm "${ID}_header.sam"
    rm "${ID}.sorted.bam"
    rm "${ID}.bam"
    rm ${sam_file}
    """
}
