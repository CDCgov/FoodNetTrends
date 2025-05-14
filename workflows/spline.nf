#!/usr/bin/env nextflow

// Import modules
include { TRENDY } from '../modules/local/trendy'
include { PREPROCESS } from '../modules/local/preprocess'

workflow SPLINE {
    // Define input channels
    if (params.pathogen) {
        // Convert comma-separated string to a channel of pathogens
        def pathogenList = params.pathogen.tokenize(',')
        pathogens = Channel.fromList(pathogenList)
    } else {
        // Default to CAMPYLOBACTER and CYCLOSPORA for testing
        pathogens = Channel.of('CAMPYLOBACTER', 'CYCLOSPORA')
    }

    // Check that required parameters are provided
    if (!params.mmwrFile) {
        error "Missing required parameter: --mmwrFile must be specified"
    }
    if (!params.censusFileB) {
        error "Missing required parameter: --censusFileB must be specified"
    }
    if (!params.censusFileP) {
        error "Missing required parameter: --censusFileP must be specified"
    }

    // Input files
    mmwrFile = file(params.mmwrFile)
    censusFileB = file(params.censusFileB)
    censusFileP = file(params.censusFileP)
    metadataFile = params.metadata ? file(params.metadata) : null
    
    // Flag for preprocessed data
    isPreprocessed = params.preprocessed ?: false

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
    if (metadataFile != null && !metadataFile.exists()) {
        error "Metadata file not found: ${params.metadata}"
    }

    // Log pipeline start
    log.info """
    ==============================================
    FoodNet Trends Pipeline
    ==============================================
    Project ID    : ${params.projID ?: new Date().format('yyyyMMdd_HHmmss')}
    MMWR File     : ${params.mmwrFile}
    Preprocessed  : ${isPreprocessed}
    Census Files  : ${params.censusFileB}, ${params.censusFileP}
    Travel        : ${params.travel}
    CIDT          : ${params.cidt}
    Pathogens     : ${params.pathogen ?: 'default (CAMPYLOBACTER,CYCLOSPORA)'}
    States        : ${params.states ?: 'all'}
    Cores         : ${params.cpus ?: 'default'}
    Chains        : ${params.chains}
    Iterations    : ${params.iterations}
    Adapt Delta   : ${params.adapt_delta}
    Max Treedepth : ${params.max_treedepth}
    Seed          : ${params.seed}
    Output Dir    : ${params.outdir}/${params.projID ?: new Date().format('yyyyMMdd_HHmmss')}
    ==============================================
    """

    // Set default projID if not specified
    def projID = params.projID ?: new Date().format('yyyyMMdd_HHmmss')

    // Always run preprocessing first if not using preprocessed data
    if (!isPreprocessed) {
        log.info "Preprocessing raw data files"
        PREPROCESS(
            mmwrFile,
            projID,
            true  // Generate metadata
        )
        
        // Use the preprocessed file for downstream analysis
        processedFile = PREPROCESS.out.cleanedData.first()
        metadataFromProcess = PREPROCESS.out.metadata.first()
    } else {
        // If using preprocessed data, skip preprocessing step
        log.info "Using preprocessed data: ${mmwrFile}"
        processedFile = mmwrFile
    }

    // Default script path if not provided
    def trendyScript = params.trendyScript ?: "${workflow.projectDir}/bin/trendy.R"

    // Run TRENDY with input data - stage metadataFile as a real input file
    TRENDY(
        pathogens,
        processedFile,
        censusFileB,
        censusFileP,
        params.travel,
        params.cidt,
        projID,
        trendyScript,
        isPreprocessed,
        metadataFile, // Pass the actual file, not just the name
        params.states
    )

    // Log completion
    log.info "Pipeline completed successfully"
}
