process TRENDY {
    tag "$pathogen"
    label 'process_large'
    shell "/bin/bash"
    container 'foodnet.sif'

    publishDir "${params.outdir}/${projID}/spline_results", mode: 'copy'

    input:
    val pathogen
    path mmwrFile
    path censusFileB
    path censusFileP
    val travel
    val cidt
    val projID
    val whichScript
    val preprocessed
    path metadataFile  // Changed to path to properly stage the file
    val states

    output:
    path "${pathogen}_brm.Rds", emit: rds, optional: true
    path "${pathogen}_IRCatch.csv", emit: csv, optional: true
    path "${pathogen}*.png", emit: png, optional: true
    path "${pathogen}*_EstIRRCatch_*.csv", emit: irr, optional: true
    path "${pathogen}_summary.txt", emit: summary, optional: true
    path "${pathogen}_error.txt", optional: true, emit: errors

    errorStrategy { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries 3

    script:
    // Simplified metadata handling - use the staged file directly
    def discoveryDataParam = metadataFile ? "--discovery_data ${metadataFile}" : ""
    
    // Add cleanFile parameter for CSV data
    def cleanFileParam = preprocessed.toString() == 'true' ? "--cleanFile ${mmwrFile}" : ""

    """
    # Copy functions.R to the current directory
    cp ${workflow.projectDir}/bin/functions.R .

    # Ensure the file was copied successfully
    if [ ! -f functions.R ]; then
        echo "Error: Failed to copy functions.R"
        exit 1
    fi

    # Debug: Show available files
    echo "Working directory contents:"
    ls -la
    echo "MMWR file path: ${mmwrFile}"
    echo "Metadata file: ${metadataFile ?: 'none'}"

    Rscript ${whichScript} \\
      --mmwrFile ${mmwrFile} \\
      --censusFileB ${censusFileB} \\
      --censusFileP ${censusFileP} \\
      --travel ${travel} \\
      --cidt ${cidt} \\
      --projID ${projID} \\
      --outDir . \\
      --pathogen ${pathogen} \\
      --preprocessed ${preprocessed} \\
      ${cleanFileParam} \\
      ${discoveryDataParam} \\
      --cores ${task.cpus} \\
      --chains ${params.chains} \\
      --iterations ${params.iterations} \\
      --adapt_delta ${params.adapt_delta} \\
      --max_treedepth ${params.max_treedepth} \\
      --seed ${params.seed} \\
      --states ${states} \\
      --debug FALSE
    """
}
