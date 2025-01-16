// Run Trimmomatic for Illumina samples LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 
process wes_trimming_s {
    publishDir "${params.outDir}/trimmed/${ID}", mode: 'copy', pattern: "*trimmed.fastq.gz"
    
    input:
    tuple val(ID), path(fastq)                                              // Input for raw samples

    output:
      tuple val(ID), path("${ID}_trimmed.fastq.gz"), emit: trimmed_fastq_s  // Output channel to store the trimmed_files
    
    script:
    """
    # Single end Trimmomatic application
    if [ "${params.seq}" == "s" ]; then
        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar SE \
            -threads ${params.t} \
            "${fastq[0]}" \
            "${ID}_trimmed.fastq.gz" \
            ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
    fi
    """
}
