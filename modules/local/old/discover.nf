process DISCOVER_DATA {
    tag "Discovering data content"
    label 'process_low'
    shell "/bin/bash"
    container 'foodnet.sif'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path mmwrFile

    output:
    path "data_inventory.json", emit: inventory

    script:
    """
    # Run the discovery script
    Rscript ${workflow.projectDir}/bin/discover_data.R ${mmwrFile} data_inventory.json
    """
}
