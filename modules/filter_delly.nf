// Filtration of SV variants obtained from Delly tools
process filter_delly {
    publishDir "${params.outDir}/variant_call/delly/${ID}", mode: 'copy', pattern: "*.vcf"

    input:
    tuple val(ID), path(input_vcf)                        // Input delly vcf file
    path ref_files                                        // Cleaned Reference genome files 
    path ref_genome                                       // Cleaned reference genome fasta
    path ref_dict                                         // Cleaned reference genome dict for GATK processes

    output:
    tuple val(ID), path("${ID}_filtered_delly.vcf"), emit: filt_delly_vcf // Output vcf file with filtered SV variants from Delly

    script:
    """
    gatk VariantFiltration \
        -R ${ref_genome} \
        -V ${input_vcf} \
        -O "${ID}_filtered_delly.vcf" \
        --filter-name "LowQual" \
        --filter-expression 'QUAL < 50' \
        --filter-name "ShortSV" \
        --filter-expression 'SVLEN < 500' \
        --filter-name "LowPESE" \
        --filter-expression 'PE < 10 || SE < 10' \
        --filter-name "LargeCIPOS" \
        --filter-expression 'CIPOS > 1000' \
        --filter-name "LargeCISTD" \
        --filter-expression 'CISTD > 1000'
    """
}
