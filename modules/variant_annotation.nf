// Process to annotate variants from Delly, CNVkit and GATK haplotypecaller and mutect2
process variant_annotation {
    publishDir "${params.outDir}/annotated/${ID}", mode: 'copy', pattern: "*.vcf"
    fair true  

    input:
    tuple val(ID), path(vcf_files)                                  // Tuple with sample ID and all VCFs
    path cache_dir
    
    output:
    tuple val(ID), path("${ID}_ann_gatk.vcf"), emit: gatk_ann       // Output annotated GATK variants
    tuple val(ID), path("${ID}_ann_cnvkit.vcf"), emit: cnvkit_ann   // Output annotated CNVkit variants
    tuple val(ID), path("${ID}_ann_delly.vcf"), emit: delly_ann     // Output annotated Delly variants

    script:
    """
    # Annotate GATK variants
    /opt/vep/src/ensembl-vep/vep -i ${vcf_files[0]} -o "${ID}_ann_cnvkit.vcf" --vcf --cache --dir_cache ${cache_dir} --merged --everything --assembly GRCh38                                 
    /opt/vep/src/ensembl-vep/vep -i ${vcf_files[1]} -o "${ID}_ann_delly.vcf" --vcf --cache --dir_cache ${cache_dir}  --merged --everything --assembly GRCh38
    /opt/vep/src/ensembl-vep/vep -i ${vcf_files[2]} -o "${ID}_ann_gatk.vcf" --vcf --cache --dir_cache ${cache_dir} --merged --everything --assembly GRCh38

    """
}
