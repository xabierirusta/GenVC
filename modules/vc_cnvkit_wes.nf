// CNV variant calling for WES with a reference BAM file
process vc_cnvkit_wes {
    publishDir "${params.outDir}/variant_call/cnvkit/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(bam_file), path(bai_file)     // BAM and BAI files input for CNV variant callling
    path ref_files                                    // Reference genome files directory
    path bed_files                                    // Target and reference BED files directory
    path refFlat                                      // Refflat annotations input
    path access_bed                                   // hg38 access bed needed by CNVkit
    path targets_bed                                  // Interesting regions BED file
    path cleaned_ref_genome                           // Cleaned reference genome fasta file
    tuple val(id), path(normal_bam), path(normal_bai) // Normal BAM and BAI files for CNVKit

    output:
    tuple val(ID), path("${ID}_call.cns"), path("${ID}.cnr"), emit: cnvkit_output // Output with the cns and cnr for plotting
    tuple val(ID), path("${ID}_cnvkit.vcf"), emit: cnvkit_vcf                     // VCF file with the variants obtained using CNVkit

    script:
    """
    # Create target and antitarget BED files for samples and the normal control
    /opt/conda/bin/cnvkit.py target ${targets_bed} --annotate ${refFlat} --split -o "${ID}.targets.bed"
    /opt/conda/bin/cnvkit.py antitarget "${ID}.targets.bed" -g ${access_bed} -o "${ID}.antitargets.bed"
    /opt/conda/bin/cnvkit.py target ${targets_bed} --annotate ${refFlat} --split -o "${id}.targets.bed"
    /opt/conda/bin/cnvkit.py antitarget "${id}.targets.bed" -g ${access_bed} -o "${id}.antitargets.bed"

    # Calculate coverage for tumor and normal files
    /opt/conda/bin/cnvkit.py coverage ${bam_file} "${ID}.targets.bed" -o "${ID}.targetcoverage.cnn"
    /opt/conda/bin/cnvkit.py coverage ${bam_file} "${ID}.antitargets.bed" -o "${ID}.antitargetcoverage.cnn"
    /opt/conda/bin/cnvkit.py coverage ${normal_bam} "${id}.targets.bed" -o "${id}_normal.targetcoverage.cnn"
    /opt/conda/bin/cnvkit.py coverage ${normal_bam} "${id}.antitargets.bed" -o "${id}_normal.antitargetcoverage.cnn"

    # Create references for normal files
    /opt/conda/bin/cnvkit.py reference "${id}_normal.targetcoverage.cnn" "${id}_normal.antitargetcoverage.cnn" -f ${cleaned_ref_genome} -o "normal_reference.cnn"

    # Fix the obtained files for tumor and normal files
    /opt/conda/bin/cnvkit.py fix "${ID}.targetcoverage.cnn" "${ID}.antitargetcoverage.cnn" "normal_reference.cnn" -o "${ID}.cnr"

    # Remove reads with 0 coverage depth and extreme log2
    awk 'NR == 1 || (\$5 > 0 && \$6 > -10 && \$6 < 10)' "${ID}.cnr" > "${ID}_cleaned.cnr"

    # Infer discrete copy number segments from the given coverage table
    /opt/conda/bin/cnvkit.py segment "${ID}_cleaned.cnr" -o "${ID}.cns"
    # Given segmented log2 ratio estimates (.cns), derive each segment’s absolute integer copy number
    /opt/conda/bin/cnvkit.py call "${ID}.cns" -o "${ID}_call.cns"

    # Export the variants information file as VCF
    /opt/conda/bin/cnvkit.py export vcf "${ID}_call.cns" -i ${ID} -o "${ID}_cnvkit.vcf"
    """
}
