#!/usr/bin/env nextflow
/*
 * ==================================================================
 * FoodNetTrends v1.0 - Preprocessing Workflow
 * ==================================================================
 *
 * Purpose:
 *   Standalone workflow for data preprocessing operations.
 *   Validates and standardizes raw MMWR surveillance data for
 *   downstream analysis without running the full modeling pipeline.
 *
 * Workflow Steps:
 *   1. Input parameter validation
 *   2. Raw data processing and standardization
 *   3. Metadata generation and output verification
 *
 * Use Cases:
 *   - Data preparation for multiple analysis runs
 *   - Quality control and data validation
 *   - Preprocessing for external analysis tools
 *
 * Last updated: 2025-05-22
 * ==================================================================
 */

// Import module
include { PREPROCESS } from '../modules/local/preprocess'

workflow PREPROCESS_WORKFLOW {
    // Log workflow version at startup
    log.info "Running FoodNetTrends Preprocessing Workflow v1.0"
    
    // Validate required parameters
    if (!params.mmwrFile) {
        error "Missing required parameter: --mmwrFile must be specified"
    }
    // Census files are required
    if (!params.censusFileB) {
        error "Census bacterial file is required. Please specify with --censusFileB parameter."
    }
    if (!params.censusFileP) {
        error "Census parasitic file is required. Please specify with --censusFileP parameter."
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
            error "ERROR: Census bacterial file does not exist: ${params.censusFileB}"
        }
    } else {
        error "ERROR: Census bacterial file parameter is required but empty"
    }
    
    if (params.censusFileP && params.censusFileP != "") {
        censusFileP = file(params.censusFileP, checkIfExists: false)
        if (!censusFileP.exists()) {
            error "ERROR: Census parasitic file does not exist: ${params.censusFileP}"
        }
    } else {
        error "ERROR: Census parasitic file parameter is required but empty"
    }
    
    // Set output base name (derived from file or parameter)
    def outputBase = ""
    try {
        if (params.outputBase) {
            outputBase = params.outputBase
        } else {
            // Extract the filename using string operations instead of File object
            def mmwrPath = params.mmwrFile.toString()
            def lastSlash = mmwrPath.lastIndexOf('/')
            if (lastSlash == -1) {
                lastSlash = mmwrPath.lastIndexOf('\\')
            }
            
            def fileName = lastSlash > -1 ? mmwrPath.substring(lastSlash + 1) : mmwrPath
            def lastDot = fileName.lastIndexOf('.')
            outputBase = lastDot > -1 ? fileName.substring(0, lastDot) : fileName
        }
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
    FoodNetTrends Preprocessing
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
    
    // Handle workflow completion - using null-safe syntax to avoid NPE
    workflow.onComplete = {
        def w = workflow
        def success = w?.success ?: false
        def status = success ? 'COMPLETED' : 'FAILED'
        def now = new Date()

        log.info """
        ==============================================
        FoodNetTrends Preprocessing: ${status}
        ==============================================
        Completed at     : ${now}
        Duration         : ${w?.duration ?: 'unknown'}
        Success          : ${success}
        Work directory   : ${w?.workDir ?: 'unknown'}
        Exit status      : ${w?.exitStatus ?: 'unknown'}
        Output path      : ${params.outdir}
        ==============================================
        """
    }
    
    // Return output channels for potential downstream use
    emit:
    cleanedData = PREPROCESS.out.cleanedData
    metadata = PREPROCESS.out.metadata
    logs = PREPROCESS.out.logs
}
