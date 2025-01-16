// Filtering haplotypecaller variants
process filter_haplotypecaller {
    publishDir "${params.outDir}/variant_call/haplotypecaller/${ID}", mode: 'copy', pattern: "*.vcf"
    
    input:
    tuple val(ID), path(input_vcf)                        // Haplotypecaller called variants vcf
    path ref_files                                        // Cleaned Reference genome files 
    path ref_genome                                       // Cleaned reference genome fasta
    path ref_dict                                         // Cleaned reference genome dict for GATK processes
    path vcf_files                                        // Complementary VCF files needed by GATK

    output:
    tuple val(ID), path("${ID}_filtered_htcaller.vcf"), emit: filt_htc_vcf // Filtered haplotypecaller variants

    script:
    """
    gatk VariantFiltration \
        -R ${ref_genome} \
        -V ${input_vcf} \
        -O "${ID}_filtered_htcaller.vcf" \
        --filter-name "QD_filter" --filter "QD < 2.0" \
        --filter-name "FS_filter" --filter "FS > 60.0" \
        --filter-name "MQ_filter" --filter "MQ < 40.0" \
        --filter-name "DP_filter" --filter "DP < 10" \
        --filter-expression "SOR > 4.0" --filter-name "HighSOR"
    """
}
