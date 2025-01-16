// Germline SNP and indel variant calling using GATK haplotypecaller
process vc_haplotypecaller {
    publishDir "${params.outDir}/variant_call/haplotypecaller/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(bam_file), path(bai_file)         // Input BAM file (sorted and indexed)
    path ref_files                                        // Cleaned reference genome related files 
    path dbsnp_vcf                                        // dbSNP VCF file for variant annotation (optional)
    path ref_genome                                       // Cleaned referenge genome fasta
    path ref_dict                                         // Cleaned referece genome dict file required by GATK
    path vcf_files                                        // VCF_files that are required by GATK

    output:
    tuple val(ID), path("${ID}_htcaller.vcf"), emit: htcaller_vcf // HaplotypeCaller VCF file
    path "${ID}_recal_data.table"                                 // Recal data table of all samples

    script:
    """
    # Generate recalibration table for Base Quality Score Recalibration (BQSR)
    gatk BaseRecalibrator \
        -I ${bam_file} \
        -R ${ref_genome}\
        --known-sites ${dbsnp_vcf} \
        -O "${ID}_recal_data.table"

    # Apply the model to adjust the base quality scores
    gatk ApplyBQSR \
        -I ${bam_file} \
        -R ${ref_genome}\
        --bqsr-recal-file "${ID}_recal_data.table" \
        -O "${ID}.bqsr.bam"

    # Run GATK HaplotypeCaller to call SNPs and Indels
    gatk HaplotypeCaller \
        -R ${ref_genome} \
        -I "${ID}.bqsr.bam" \
        -O "${ID}_htcaller.vcf" 

    rm "${ID}.bqsr.bam"
    """
}
