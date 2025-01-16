// Somatic SNP and indel calling using GATK Mutect2
process vc_mutect2 {
    publishDir "${params.outDir}/variant_call/mutect2/${ID}", mode: 'copy', patter: "*.vcf"

    input:
    tuple val(ID), path(bam_file), path(bai_file)         // Input BAM file (sorted and indexed)
    tuple val(id), path(normal_bam), path(normal_bai)     // Input Normal BAM and BAI files 
    path ref_files                                        // Cleaned reference genome realted files
    path ref_genome                                       // Cleaned reference genome fasta file
    path dbsnp_vcf                                        // dbSNP VCF file for variant annotation (optional)
    path gnomad_ann                                       // gnomAD database vcf file
    path pon_hg38                                         // Panel of normals for hg38 genome
    path ref_dict                                         // Cleaned reference genome dict file required by GATK
    path vcf_files                                        // VCF files that may be required by GATK

    output:
    tuple val(ID), path("${ID}_mutect2.vcf.gz"), emit: mutect2_vcf // Mutect2 vcf file
    path "${ID}_mutect2.vcf.gz.stats", emit: mutect2_stats         // Mutect2 stats file
    
    script:
    """
    # Generate recalibration table for Base Quality Score Recalibration (BQSR) for the samples and normal files
    gatk BaseRecalibrator \
    -I ${bam_file} \
    -R ${ref_genome} \
    --known-sites ${dbsnp_vcf}  \
    -O "${ID}_recal_data.table"

    gatk BaseRecalibrator \
    -I ${normal_bam} \
    -R ${ref_genome} \
    --known-sites ${dbsnp_vcf}  \
    -O "${id}_recal_data.table"

    # Apply the model to adjust the base quality scores for the samples and normal files
    gatk ApplyBQSR \
    -I ${bam_file} \
    -R ${ref_genome} \
    --bqsr-recal-file "${ID}_recal_data.table" \
    -O "${ID}.bqsr.bam"

    gatk ApplyBQSR \
    -I ${normal_bam} \
    -R ${ref_genome} \
    --bqsr-recal-file "${id}_recal_data.table" \
    -O "${id}.bqsr.bam"

    # Run GATK Mutect2 
    gatk Mutect2 \
     -R ${ref_genome} \
     -I "${ID}.bqsr.bam" \
     -I "${id}.bqsr.bam" \
     -normal ${id} \
     --germline-resource ${gnomad_ann} \
     --panel-of-normals ${pon_hg38} \
     -O "${ID}_mutect2.vcf.gz"

     rm "${ID}.bqsr.bam"
     rm "${id}.bqsr.bam"
    """
}
