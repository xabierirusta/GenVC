// Run Trimmomatic for Illumina samples LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
process wgs_trimming_s {
  publishDir "${params.outDir}/trimmed/${ID}", mode: 'copy', pattern: "*.fastq.gz"

  input:
  tuple val(ID), path(fastq)  // Input for raw samples and sequencing type

  output:
    tuple val(ID), path("${ID}_trimmed.fastq.gz"), emit: trimmed_fastq_s  // Output channel for trimmed files
    
    script:
    """   
    # Single-ended trimmomatic process
    if [ "${params.seq}" == "s" ]; then
        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar SE \
            -threads ${params.t} \
            "${fastq[0]}" \
            "${ID}_trimmed.fastq.gz" \
            ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
    fi
    """
}