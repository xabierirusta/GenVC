// Alignment for paired-ended data using bwa-mem 
process alignment {
  publishDir "${params.outDir}/alignment/${ID}", mode: 'copy', pattern: "*.sam"

    input:
    tuple val(ID), path(trimmed_files_1), path(trimmed_files_2) // Trimmed files input
    path ref_files                                              // Reference genome index files input
    path ref_genome                                             // Reference genome fasta input

    output:
    tuple val(ID), path("${ID}.sam"), emit: aligned_sam         // Create an output channel for the aligned .SAM files

    script:
    """
    bwa mem -M -t ${params.t} ${ref_genome} ${trimmed_files_1} ${trimmed_files_2} > "${ID}.sam"
    """
}
