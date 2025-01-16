// Exonic alignment filtration for WES
process wes_alignment_filtration {
    publishDir "${params.outDir}/alignment/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(sorted_bam_file) // Sorted BAM file with its ID
    path bed_files                       // BED files required in the process
    path targets_bed                     // Target regions BED file
    path ref_dict                        // Cleaned reference genome dict file

    output:
    tuple val(ID), path("${ID}_RG_targets.bam"), path("${ID}_RG_targets.bam.bai"), emit: filtrated_bam_bai // BAM & BAI files filtrated with target exonic sequences
   
    script:
    """
    # Extract chromosome names from the dictionary, add 'chr' prefix
    chromosomes=\$(grep "^@SQ" ${ref_dict} | cut -f2 | sed 's/SN://g' | sed 's/^chr//')

    # Modify the header of the filtered BAM file (replace SN:1 with SN:chr1, SN:2 with SN:chr2, etc.)
    samtools view -H ${sorted_bam_file} > "${ID}_header.sam"
    for i in \$chromosomes; do
        sed -i "s/SN:\$i/SN:chr\$i/g" "${ID}_header.sam"
    done

    # Reapply the modified header to the reheadered BAM file
    samtools reheader "${ID}_header.sam" ${sorted_bam_file} > "${ID}_RG_filt.bam"

    # Filter the BAM file using the target regions BED file
    samtools view -L ${targets_bed} -b "${ID}_RG_filt.bam" > "${ID}_RG_targets.bam"
    samtools index "${ID}_RG_targets.bam"

    rm ${sorted_bam_file}
    rm "${ID}_header.sam"
    rm "${ID}_RG_filt.bam"
    """
}
