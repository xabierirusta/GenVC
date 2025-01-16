// Structural variants calling using Delly
process vc_delly {
    publishDir "${params.outDir}/variant_call/delly/${ID}", mode: 'copy', pattern: "*.vcf"

    input:
    tuple val(ID), path(bam_file), path(bai_file)         // Input BAM file (sorted and indexed)
    path ref_files                                        // Cleaned Reference genome related files 
    path ref_genome                                       // Cleaned Reference genome fasta

    output:
    tuple val(ID), path("${ID}_delly.vcf"), emit: delly_vcf // Output vcf file with SV variants from Delly

    script:
    """
    # Variant calling of SV using Delly
    /opt/delly/bin/delly call -g ${ref_genome} ${bam_file} > "${ID}_delly.vcf"
    """
}
