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
    path cleanFile

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
    // Properly handle the cleanFile parameter
    def cleanFileParam = ""
    if (preprocessed) {
        if (cleanFile) {
            cleanFileParam = "--cleanFile ${cleanFile}"
        } else {
            log.warn "Preprocessing enabled but no clean file provided for pathogen: ${pathogen}"
        }
    }

    """
    # Copy functions.R to the current directory
    cp ${workflow.projectDir}/bin/functions.R .

    # Ensure the file was copied successfully
    if [ ! -f functions.R ]; then
        echo "Error: Failed to copy functions.R"
        exit 1
    fi

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
      --cores ${task.cpus} \\
      --chains ${params.chains} \\
      --iterations ${params.iterations} \\
      --adapt_delta ${params.adapt_delta} \\
      --max_treedepth ${params.max_treedepth} \\
      --seed ${params.seed} \\
      --debug FALSE
    """
}
