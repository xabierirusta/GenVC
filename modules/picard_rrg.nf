// Bam_processing for WES 
process picard_rrg {
    publishDir "${params.outDir}/alignment/${ID}", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(ID), path(bam_file)                      // Aligned SAM file input with its ID

    output:
    tuple val(ID), path("${ID}_RG.bam"), emit: rrg_bam // Output channel with sorted BAM files

    script:
    """
    java -jar /usr/picard/picard.jar AddOrReplaceReadGroups \
        I=${bam_file} \
        O="${ID}_RG.bam" \
        RGID=${ID} \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        RGSM=${ID}

    rm ${bam_file}
    """
}