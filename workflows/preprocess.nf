#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FoodNet Trends Preprocessing Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow handles data preprocessing for the FoodNet Trends pipeline.
    
    Main steps:
    1. Validate input parameters
    2. Process MMWR surveillance data
    3. Generate clean CSV and optional metadata
    
    Last updated: 2025-05-18
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import module
include { PREPROCESS } from '../modules/local/preprocess'

workflow PREPROCESS_WORKFLOW {
    // Log workflow version at startup
    log.info "Running FoodNet Trends Preprocessing Workflow v1.0"
    
    // Validate required parameters
    if (!params.mmwrFile) {
        error "Missing required parameter: --mmwrFile must be specified"
    }
    if (!params.censusFileB) {
        error "Missing required parameter: --censusFileB must be specified"
    }
    if (!params.censusFileP) {
        error "Missing required parameter: --censusFileP must be specified"
    }
    if (!params.outdir) {
        error "Missing required parameter: --outdir must be specified"
    }
    
    // Define input channel
    mmwrFile = file(params.mmwrFile, checkIfExists: true)
    censusFileB = file(params.censusFileB, checkIfExists: true)
    censusFileP = file(params.censusFileP, checkIfExists: true)
    
    // Set output base name (derived from file or parameter)
    outputBase = params.outputBase ?: file(params.mmwrFile).getBaseName()

    // Set metadata generation flag (default to true for this workflow)
    generateMetadata = params.generateMetadata ?: true

    // Check if file exists
    if (!mmwrFile.exists()) {
        error "MMWR file not found: ${params.mmwrFile}"
    }
    if (!censusFileB.exists()) {
        error "Census bacterial file not found: ${params.censusFileB}"
    }
    if (!censusFileP.exists()) {
        error "Census parasitic file not found: ${params.censusFileP}"
    }

    // Log preprocessing start
    log.info """
    ==============================================
    FoodNet Trends Preprocessing
    ==============================================
    MMWR File         : ${params.mmwrFile} 
                        (${mmwrFile.size() >= 1024*1024 ? 
                            String.format('%.2f MB', mmwrFile.size()/(1024*1024)) : 
                            String.format('%.2f KB', mmwrFile.size()/1024)})
    Census File (B)   : ${params.censusFileB}
    Census File (P)   : ${params.censusFileP}
    Output Base       : ${outputBase}
    Generate Metadata : ${generateMetadata}
    Output Dir        : ${params.outdir}/preprocessed
    Nextflow Version  : ${nextflow.version}
    Starting time     : ${new Date()}
    ==============================================
    """

    // Run the preprocessing
    PREPROCESS(
        mmwrFile,
        censusFileB,
        censusFileP,
        outputBase,
        generateMetadata
    )
    
    // Handle workflow completion
    workflow.onComplete {
        log.info """
        ==============================================
        FoodNet Trends Preprocessing: ${workflow.success ? 'COMPLETED' : 'FAILED'}
        ==============================================
        Completed at     : ${new Date()}
        Duration         : ${workflow.duration}
        Success          : ${workflow.success}
        Work directory   : ${workflow.workDir}
        Exit status      : ${workflow.exitStatus}
        ==============================================
        """
    }
    
    // Return output channels for potential downstream use
    emit:
    cleanedData = PREPROCESS.out.cleanedData
    metadata = PREPROCESS.out.metadata
    logs = PREPROCESS.out.logs
}
