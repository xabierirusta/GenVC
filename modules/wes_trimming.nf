// Run Trimmomatic for Illumina samples LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 
process wes_trimming {
    publishDir "${params.outDir}/trimmed/${ID}", mode: 'copy', pattern: "*paired.fastq.gz"

    input:
    tuple val(ID), path(fastq)  // Input for raw samples and sequencing type

    output:
    tuple val(ID), path("${ID}_1_trimmed_paired.fastq.gz"), path("${ID}_2_trimmed_paired.fastq.gz"), emit: trimmed_fastq  // Output channel to store the trimmed_files
    
    script:
    """
    # Paired-ended Trimmomatic application
    if [ "${params.seq}" == "p" ]; then
        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE \
            -threads ${params.t} \
            "${fastq[0]}" "${fastq[1]}" \
            "${ID}_1_trimmed_paired.fastq.gz" \
            "${ID}_1_trimmed_unpaired.fastq.gz" \
            "${ID}_2_trimmed_paired.fastq.gz" \
            "${ID}_2_trimmed_unpaired.fastq.gz" \
            ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
    fi
    """
}
