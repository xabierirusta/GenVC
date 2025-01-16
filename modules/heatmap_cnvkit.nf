// Process to create plots and heatmaps for the data obtained from CNVkit
process heatmap_cnvkit {
    publishDir "${params.outDir}/plots/cnvkit", mode: 'copy'

    input:
    path cns_files                                  // cns files of all samples obtained from CNVkit

    output:
    path "cnvkit_heatmap.pdf", emit: cnvkit_heatmap // Output heatmap

    script:
    """
    /opt/conda/bin/cnvkit.py heatmap ${cns_files.join(' ')} -o "cnvkit_heatmap.pdf"
    """
}
