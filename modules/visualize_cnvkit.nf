// Process to create plots from the data obtained from CNVkit
process visualize_cnvkit {
    publishDir "${params.outDir}/plots/cnvkit/${ID}", mode: 'copy'

    input:
    tuple val(ID), path(cns_file), path(cnr_file)  // Input of files obtained from CNVkit analysis

    output:
    tuple val(ID), path("${ID}_scatter.pdf"), path("${ID}_diagram.pdf"), emit: cnvkit_plots // Output plots

    script:
    """
    # Plots 
    /opt/conda/bin/cnvkit.py scatter -s ${cns_file} ${cnr_file} -o "${ID}_scatter.pdf"
    /opt/conda/bin/cnvkit.py diagram -s ${cns_file} ${cnr_file} -o "${ID}_diagram.pdf"
    """
}
