// CNV variant calling for WES without having a reference BAM file
process vc_cnvkit_wes_noref {
    publishDir "${params.outDir}/variant_call/cnvkit/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(bam_file), path(bai_file) // BAM and BAI files input for CNV variant callling
    path ref_files                                // Reference genome files directory
    path bed_files                                // Target and reference BED files directory
    path refFlat                                  // Refflat annotations input
    path access_bed                               // hg38 access bed needed by CNVkit
    path targets_bed                              // Interesting regions BED file
    path cleaned_ref_genome                       // Cleaned reference genome fasta file

    output:
    tuple val(ID), path("${ID}_call.cns"), path("${ID}_cleaned.cnr"), emit: cnvkit_output // Output with the cns and cnr for plotting
    tuple val(ID), path("${ID}_cnvkit.vcf"), emit: cnvkit_vcf                             // VCF file with the variants obtained using CNVkit
 
    script:
    """
    # Create target and antitarget BED files
    /opt/conda/bin/cnvkit.py target ${targets_bed} --annotate ${refFlat} --split -o "${ID}_targets.bed"
    /opt/conda/bin/cnvkit.py antitarget "${ID}_targets.bed" -g ${access_bed} -o "${ID}_antitargets.bed"

    # Calculate the coverage for the samples
    /opt/conda/bin/cnvkit.py coverage ${bam_file} "${ID}_targets.bed" -o "${ID}.targetcoverage.cnn"
    /opt/conda/bin/cnvkit.py coverage ${bam_file} "${ID}_antitargets.bed" -o "${ID}.antitargetcoverage.cnn"

    # Create a flat reference without normal control samples
    /opt/conda/bin/cnvkit.py reference "${ID}.targetcoverage.cnn" "${ID}.antitargetcoverage.cnn" -f ${cleaned_ref_genome} -o "${ID}_Reference.cnn"

    # Obtain copy number ratios table
    /opt/conda/bin/cnvkit.py fix "${ID}.targetcoverage.cnn" "${ID}.antitargetcoverage.cnn" "${ID}_Reference.cnn" -o "${ID}.cnr"
    # Remove reads with 0 coverage depth and extreme log2
    awk 'NR == 1 || (\$5 > 0 && \$6 > -10 && \$6 < 10)' "${ID}.cnr" > "${ID}_cleaned.cnr"

    # Infer discrete copy number segments from the given coverage table
    /opt/conda/bin/cnvkit.py segment "${ID}_cleaned.cnr" -o "${ID}.cns"
    # Given segmented log2 ratio estimates (.cns), derive each segment’s absolute integer copy number
    /opt/conda/bin/cnvkit.py call "${ID}.cns" -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o "${ID}_call.cns"

    # Export the variants information file as VCF
    /opt/conda/bin/cnvkit.py export vcf "${ID}_call.cns" -i ${ID} -o "${ID}_cnvkit.vcf"
    """
}
