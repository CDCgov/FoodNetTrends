#!/usr/bin/env nextflow

// Import modules
include { DISCOVER_DATA } from '../modules/local/discover'

workflow DISCOVER {
    // Define input channel
    mmwrFile = file(params.mmwrFile)

    // Check if file exists
    if (!mmwrFile.exists()) {
        error "MMWR file not found: ${params.mmwrFile}"
    }

    // Log discovery start
    log.info """
    ==============================================
    FoodNet Trends Data Discovery
    ==============================================
    MMWR File     : ${params.mmwrFile}
    Output Dir    : ${params.outdir}
    ==============================================
    """

    // Run the discovery process
    DISCOVER_DATA(
        mmwrFile
    )
}
