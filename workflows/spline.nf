#!/usr/bin/env nextflow

// Import modules
include { TRENDY } from '../modules/local/trendy'
include { PREPROCESS } from '../modules/local/preprocess'

workflow SPLINE {
    // Define input channels
    if (params.pathogen) {
        // If specific pathogen is requested, use it
        pathogens = Channel.of(params.pathogen)
    } else {
        // Default to CAMPYLOBACTER and CYCLOSPORA for testing
        pathogens = Channel.of('CAMPYLOBACTER', 'CYCLOSPORA')
    }

    // Input files
    mmwrFile = file(params.mmwrFile)
    censusFileB = file(params.censusFileB)
    censusFileP = file(params.censusFileP)

    // Check if files exist
    if (!mmwrFile.exists()) {
        error "MMWR file not found: ${params.mmwrFile}"
    }
    if (!censusFileB.exists()) {
        error "Census bacterial file not found: ${params.censusFileB}"
    }
    if (!censusFileP.exists()) {
        error "Census parasitic file not found: ${params.censusFileP}"
    }

    // Log pipeline start
    log.info """
    ==============================================
    FoodNet Trends Pipeline
    ==============================================
    Project ID    : ${params.projID}
    MMWR File     : ${params.mmwrFile}
    Census Files  : ${params.censusFileB}, ${params.censusFileP}
    Travel        : ${params.travel}
    CIDT          : ${params.cidt}
    Cores         : ${params.cpus ?: 'default'}
    Chains        : ${params.chains}
    Iterations    : ${params.iterations}
    Adapt Delta   : ${params.adapt_delta}
    Max Treedepth : ${params.max_treedepth}
    Seed          : ${params.seed}
    Output Dir    : ${params.outdir}/${params.projID}
    ==============================================
    """

    // Conditional preprocessing
    if (params.preprocessed) {
        log.info "Using preprocessed data from: ${params.cleanFile}"
        cleanFile = file(params.cleanFile)
        if (!cleanFile.exists()) {
            error "Preprocessed file not found: ${params.cleanFile}"
        }

        // Run TRENDY with preprocessed data
        TRENDY(
            pathogens,
            mmwrFile,
            censusFileB,
            censusFileP,
            params.travel,
            params.cidt,
            params.projID,
            params.trendyScript,
            params.preprocessed,
            cleanFile
        )
    } else {
        log.info "Preprocessing raw data files"

        // Run preprocessing step
        PREPROCESS(
            mmwrFile,
            params.projID
        )

        // Create a proper channel from the preprocessed file
        processedFile = PREPROCESS.out.cleanFile.first()

        // Run TRENDY with processed data
        TRENDY(
            pathogens,
            mmwrFile,
            censusFileB,
            censusFileP,
            params.travel,
            params.cidt,
            params.projID,
            params.trendyScript,
            params.preprocessed,
            processedFile
        )
    }

    // Log completion
    log.info "Pipeline completed successfully"
}
