// CNV variant calling for WGS
process vc_cnvkit_wgs {
    publishDir "${params.outDir}/variant_call/cnvkit/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(bam_file), path(bai_file)     // BAM and BAI files input for CNV variant callling
    path ref_files                                    // Reference genome files directory
    path refFlat                                      // Refflat annotations input
    path bed_files                                    // Target region Bed file input
    path cleaned_ref_genome                           // Cleaned reference genome fasta file
    path targets_bed                                  // Interesting regions BED file
    tuple val(id), path(normal_bam), path(normal_bai) // Normal BAM and BAI for CNVkit

    output:
    tuple val(ID), path("${ID}_mdup.cns"), path("${ID}_mdup.cnr"), emit: cnvkit_output  // Output with the cns and cnr for plotting
    tuple val(ID), path("${ID}_cnvkit.vcf"), emit: cnvkit_vcf                 // VCF file with the variants obtained using CNVkit

    script:
    """
    /opt/conda/bin/cnvkit.py batch ${bam_file} -n ${normal_bam} \
        -m wgs -f ${cleaned_ref_genome} -t ${targets_bed} --annotate ${refFlat}

    # Export the variants information file as VCF
    /opt/conda/bin/cnvkit.py export vcf "${ID}_mdup.call.cns" -i ${ID} -o "${ID}_cnvkit.vcf"
    """
}
