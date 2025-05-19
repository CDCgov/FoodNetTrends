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
    // Census files are now optional with warnings
    if (!params.censusFileB) {
        log.warn "No census bacterial file specified. Using placeholder."
    }
    if (!params.censusFileP) {
        log.warn "No census parasitic file specified. Using placeholder."
    }
    if (!params.outdir) {
        error "Missing required parameter: --outdir must be specified"
    }
    
    // Define input channel
    def mmwrFilePath = params.mmwrFile
    mmwrFile = file(mmwrFilePath, checkIfExists: false)
    
    // Handle empty census file parameters
    def censusFileB = ""
    def censusFileP = ""
    
    // Check if census file parameters are not empty strings
    if (params.censusFileB && params.censusFileB != "") {
        censusFileB = file(params.censusFileB, checkIfExists: false)
        if (!censusFileB.exists()) {
            log.warn "WARNING: Census bacterial file does not exist: ${params.censusFileB}"
            log.warn "Will proceed with placeholder census bacterial file"
        }
    } else {
        log.warn "WARNING: Census bacterial file parameter is empty"
        log.warn "Will proceed with placeholder census bacterial file"
    }
    
    if (params.censusFileP && params.censusFileP != "") {
        censusFileP = file(params.censusFileP, checkIfExists: false)
        if (!censusFileP.exists()) {
            log.warn "WARNING: Census parasitic file does not exist: ${params.censusFileP}"
            log.warn "Will proceed with placeholder census parasitic file"
        }
    } else {
        log.warn "WARNING: Census parasitic file parameter is empty"
        log.warn "Will proceed with placeholder census parasitic file"
    }
    
    // Set output base name (derived from file or parameter)
    def outputBase = ""
    try {
        outputBase = params.outputBase ?: new File(params.mmwrFile).getName().replaceFirst("[.][^.]+\$", "")
    } catch (Exception e) {
        log.warn "Could not determine base name from file: ${e.message}"
        outputBase = "foodnet_data_" + new Date().format('yyyyMMdd_HHmmss')
    }
    
    // Set metadata generation flag (default to true for this workflow)
    def generateMetadata = params.generateMetadata ?: true
    
    // Check if MMWR file exists
    if (!mmwrFile.exists()) {
        error "MMWR file not found: ${params.mmwrFile}"
    }
    
    // Log preprocessing start
    log.info """
    ==============================================
    FoodNet Trends Preprocessing
    ==============================================
    MMWR File         : ${params.mmwrFile}
    Census File (B)   : ${params.censusFileB} ${censusFileB && censusFileB.exists() ? "✓" : "✗"}
    Census File (P)   : ${params.censusFileP} ${censusFileP && censusFileP.exists() ? "✓" : "✗"}
    Output Base       : ${outputBase}
    Generate Metadata : ${generateMetadata}
    Output Dir        : ${params.outdir}/preprocessed
    Nextflow Version  : ${nextflow.version}
    Starting time     : ${new Date()}
    ==============================================
    """
    
    // Run the preprocessing
    try {
        PREPROCESS(
            mmwrFile,
            censusFileB,
            censusFileP,
            outputBase,
            generateMetadata
        )
    } catch (Exception e) {
        log.error "Error in preprocessing: ${e.message}"
        throw e
    }
    
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
