#!/usr/bin/env nextflow

// Import module
include { PREPROCESS } from '../modules/local/preprocess'

workflow PREPROCESS_WORKFLOW {
    // Define input channel
    mmwrFile = file(params.mmwrFile)
    
    // Set output base name (derived from file or parameter)
    outputBase = params.outputBase ?: file(params.mmwrFile).getBaseName()

    // Set metadata generation flag (default to true for this workflow)
    generateMetadata = params.generateMetadata ?: true

    // Check if file exists
    if (!mmwrFile.exists()) {
        error "MMWR file not found: ${params.mmwrFile}"
    }

    // Log preprocessing start
    log.info """
    ==============================================
    FoodNet Trends Preprocessing
    ==============================================
    MMWR File         : ${params.mmwrFile}
    Output Base       : ${outputBase}
    Generate Metadata : ${generateMetadata}
    Output Dir        : ${params.outdir}/preprocessed
    ==============================================
    """

    // Run the preprocessing
    PREPROCESS(
        mmwrFile,
        outputBase,
        generateMetadata
    )
}
