// Filtering mutect2 variants
process filter_mutect2 {
    publishDir "${params.outDir}/variant_call/mutect2/${ID}", mode: 'copy', pattern: "*.vcf"

    input:
    tuple val(ID), path(input)                            // Tumor BAM and BAI files input
    path gnomad_ann                                       // gnomAD database vcf file
    path ref_files                                        // Cleaned Reference genome files 
    path ref_cleaned_genome                               // Cleaned reference genome fasta
    path ref_dict                                         // Cleaned reference genome dict for GATK processes
    path vcf_files                                        // VCF files required by GATK
    path contigs                                          // .list file to filter the interesing contigs from the VCF files
    path mutect2stats                                     // Stats file obtained from GATK Mutect2

    output:
    tuple val(ID), path("${ID}_filtered_mutect2.vcf"), emit: filt_mutect2_vcf // Filtered mutect2 variants

    script:
    """
    gatk IndexFeatureFile -I ${input[2]}

    gatk --java-options "-Xmx16g" GetPileupSummaries \
        -I ${input[0]} \
        -V ${gnomad_ann} \
        -L ${contigs} \
        -O "${ID}_tumor-pileups.table"

    gatk --java-options "-Xmx16g" CalculateContamination \
        -I "${ID}_tumor-pileups.table" \
        -O "${ID}_contamination.table"

    gatk --java-options "-Xmx16g" FilterMutectCalls \
        -V ${input[2]} \
        --contamination-table "${ID}_contamination.table" \
        --reference ${ref_cleaned_genome} \
        --stats ${mutect2stats} \
        -O "${ID}_filtered_mutect2.vcf"
    """
}
