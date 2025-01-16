// Process to index all the annotated vcf files 
process final_index {
    publishDir "${params.outDir}/annotated/${ID}", mode: 'copy'
    fair true                                                // Keep the order of the files

    input:
    tuple val(ID), path(ann_vcf)                             // Input annotated VCF files

    output:
    tuple val(ID), path("${ID}_ann_gatk.vcf.gz"), path("${ID}_ann_gatk.vcf.gz.tbi"), emit: gatk_ann       //Output annotated GATK variants + index
    tuple val(ID), path("${ID}_ann_cnvkit.vcf.gz"), path("${ID}_ann_cnvkit.vcf.gz.tbi"), emit: cnvkit_ann   //Output annotated CNVkit variants + index
    tuple val(ID), path("${ID}_ann_delly.vcf.gz"), path("${ID}_ann_delly.vcf.gz.tbi"), emit: delly_ann     //Output annotated Delly variants + index

    script:
    """
    # Process GATK VCF file
    bcftools view -Oz -f PASS -o "${ID}_ann_gatk.vcf.gz" ${ann_vcf[2]}
    bcftools index -t "${ID}_ann_gatk.vcf.gz" -o "${ID}_ann_gatk.vcf.gz.tbi"

    # Process CNVkit VCF file
    bcftools view -Oz -o "${ID}_ann_cnvkit.vcf.gz" ${ann_vcf[0]}
    bcftools index -t "${ID}_ann_cnvkit.vcf.gz" -o "${ID}_ann_cnvkit.vcf.gz.tbi"

    # Process Delly VCF file
    bcftools view -Oz -f PASS -o "${ID}_ann_delly.vcf.gz" ${ann_vcf[1]}
    bcftools index -t "${ID}_ann_delly.vcf.gz" -o "${ID}_ann_delly.vcf.gz.tbi"
    """
    
}
