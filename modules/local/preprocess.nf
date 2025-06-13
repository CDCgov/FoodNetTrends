process PREPROCESS {
    tag "Preprocessing MMWR data"
    label 'process_medium'
    shell "/bin/bash"
    container 'foodnet.sif'

    publishDir "${params.outdir}/${projID}/preprocessed", mode: 'copy'

    input:
    path mmwrFile
    val projID
    path serotypeConfig

    output:
    path "clean_mmwr.csv", emit: cleanFile
    path "clean_mmwr_preprocessing_report.csv", emit: preprocessingReport

    script:
    // Use absolute path to the script or a relative path from the current directory
    def scriptPath = "${workflow.projectDir}/bin/preprocess.R"
    def serotypeConfigArg = serotypeConfig.name != 'NO_FILE' ? "--serotype-config ${serotypeConfig}" : ""
    
    // Validate matching sensitivity
    def validSensitivities = ['STRICT', 'MEDIUM', 'RELAXED']
    if (!validSensitivities.contains(params.matching_sensitivity)) {
        error "Invalid matching_sensitivity: ${params.matching_sensitivity}. Must be one of: ${validSensitivities.join(', ')}"
    }

    """
    Rscript ${scriptPath} \\
      --mmwrFile ${mmwrFile} \\
      --outputFile clean_mmwr.csv \\
      --matching-sensitivity ${params.matching_sensitivity} \\
      ${serotypeConfigArg}
    """
}
