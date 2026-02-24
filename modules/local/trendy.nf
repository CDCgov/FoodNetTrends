process TRENDY {
    tag "$pathogen"
    label 'process_large'
    shell "/bin/bash"
    container 'foodnet.sif'

    publishDir "${params.outdir}/${projID}/spline_results", mode: 'copy'

    input:
    tuple val(pathogenGrouping), val(pathogen), val(subgroup), val(dataMetrics)
    path mmwrFile
    path censusFileB
    path censusFileP
    val travel
    val cidt
    val states
    val projID
    val whichScript
    val preprocessed
    path cleanFile
    path catchmentConfig

    output:
    path "${pathogenGrouping.replaceAll('[^a-zA-Z0-9_-]', '_').replaceAll('_+', '_').replaceAll('_$', '')}_brm.Rds", emit: rds, optional: true
    path "${pathogenGrouping.replaceAll('[^a-zA-Z0-9_-]', '_').replaceAll('_+', '_').replaceAll('_$', '')}_IRCatch.csv", emit: csv, optional: true
    path "${pathogenGrouping.replaceAll('[^a-zA-Z0-9_-]', '_').replaceAll('_+', '_').replaceAll('_$', '')}*.png", emit: png, optional: true
    path "${pathogenGrouping.replaceAll('[^a-zA-Z0-9_-]', '_').replaceAll('_+', '_').replaceAll('_$', '')}*_EstIRRCatch_*.csv", emit: irr, optional: true
    path "${pathogenGrouping.replaceAll('[^a-zA-Z0-9_-]', '_').replaceAll('_+', '_').replaceAll('_$', '')}_summary.txt", emit: summary, optional: true
    path "${pathogenGrouping.replaceAll('[^a-zA-Z0-9_-]', '_').replaceAll('_+', '_').replaceAll('_$', '')}_error.txt", optional: true, emit: errors

    errorStrategy { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries 3

    script:
    // Set data metrics in task.ext for resource allocation
    task.ext.dataMetrics = dataMetrics
    
    // Log resource allocation for this pathogen
    log.info "Pathogen: ${pathogen}, Rows: ${dataMetrics?.rows ?: 'unknown'}, " +
             "Complexity: ${dataMetrics?.complexity ?: 'unknown'}, " +
             "Allocated CPUs: ${task.cpus}, Memory: ${task.memory}"
    
    // Properly handle the cleanFile parameter
    def cleanFileParam = ""
    if (preprocessed) {
        if (cleanFile && cleanFile.name != 'NO_FILE') {
            cleanFileParam = "--cleanFile ${cleanFile}"
        } else {
            error "Preprocessing enabled but no clean file provided for pathogen: ${pathogen}"
        }
    }
    
    // Handle catchment config parameter
    def catchmentConfigArg = catchmentConfig.name != 'NO_FILE' ? "--catchment-config ${catchmentConfig}" : ""
    
    // Handle states parameter (empty string means all states)
    def statesArg = states ? "--states ${states}" : ""

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
      ${statesArg} \\
      --projID ${projID} \\
      --outDir . \\
      --pathogen ${pathogen} \\
      --subgroup '${subgroup}' \\
      --preprocessed ${preprocessed} \\
      ${cleanFileParam} \\
      --cores ${task.cpus} \\
      --chains ${params.chains} \\
      --iterations ${params.iterations} \\
      --adapt_delta ${params.adapt_delta} \\
      --max_treedepth ${params.max_treedepth} \\
      --seed ${params.seed} \\
      ${catchmentConfigArg} \\
      --debug FALSE
    """
}
