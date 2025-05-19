#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FoodNet Trends Spline Analysis Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow implements Bayesian hierarchical spline models to analyze trends
    in FoodNet surveillance data.
    
    Main steps:
    1. Validate input parameters and files
    2. Preprocess data if using raw input
    3. Run Bayesian modeling for each pathogen in parallel
    4. Generate incidence rate estimates and visualizations
    5. Create an interactive HTML dashboard for result exploration
    
    Last updated: 2025-05-18
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import modules
include { TRENDY } from '../modules/local/trendy'
include { PREPROCESS } from '../modules/local/preprocess'
include { GENERATE_DASHBOARD } from '../modules/local/generate_dashboard'

workflow SPLINE {
    // Log workflow version at startup
    log.info "Running FoodNet Trends Spline Analysis Workflow v1.0"
    
    // Define input channels
    if (params.pathogen) {
        // Convert comma-separated string to a channel of pathogens
        def pathogenList = params.pathogen.tokenize(',')
        pathogens = Channel.fromList(pathogenList).map { pathogen -> [pathogen, mmwrFile] }
        
        // Log the pathogens being analyzed
        log.info "Analyzing ${pathogenList.size()} pathogens: ${pathogenList.join(', ')}"
    } else {
        // Default to CAMPYLOBACTER and CYCLOSPORA for testing
        pathogens = Channel.of(
            ['CAMPYLOBACTER', mmwrFile], 
            ['CYCLOSPORA', mmwrFile]
        )
        log.info "No pathogens specified, using defaults: CAMPYLOBACTER, CYCLOSPORA"
    }

    // Validate required parameters
    def missingParams = []
    if (!params.mmwrFile) missingParams << "--mmwrFile"
    if (!params.censusFileB) missingParams << "--censusFileB"
    if (!params.censusFileP) missingParams << "--censusFileP"
    
    if (missingParams.size() > 0) {
        error "Missing required parameter(s): ${missingParams.join(', ')}"
    }

    // Input files
    mmwrFile = file(params.mmwrFile, checkIfExists: true)
    censusFileB = file(params.censusFileB, checkIfExists: true)
    censusFileP = file(params.censusFileP, checkIfExists: true)
    
    // Flag for preprocessed data
    isPreprocessed = params.preprocessed ?: false

    // Check file sizes to catch obvious issues
    if (mmwrFile.size() == 0) {
        error "MMWR file is empty: ${params.mmwrFile}"
    }
    if (censusFileB.size() == 0) {
        error "Census bacterial file is empty: ${params.censusFileB}"
    }
    if (censusFileP.size() == 0) {
        error "Census parasitic file is empty: ${params.censusFileP}"
    }
    
    // Get metadata file with alternate path fallback
    if (params.metadata) {
        metadataFile = findMetadataFile(params.metadata)
    } else {
        metadataFile = null
    }
    
    // Dashboard template file
    dashboardTemplate = file("${workflow.projectDir}/assets/dashboard_template.html", checkIfExists: true)
    dashboardScript = file("${workflow.projectDir}/bin/generate_dashboard.R", checkIfExists: true)

    // Set default projID if not specified
    def projID = params.projID ?: new Date().format('yyyyMMdd_HHmmss')

    // Log pipeline start
    log.info """
    ==============================================
    FoodNet Trends Spline Analysis
    ==============================================
    Project ID    : ${projID}
    MMWR File     : ${params.mmwrFile} (${formatSize(mmwrFile.size())})
    Preprocessed  : ${isPreprocessed}
    Census Files  : 
      Bacterial   : ${params.censusFileB} (${formatSize(censusFileB.size())})
      Parasitic   : ${params.censusFileP} (${formatSize(censusFileP.size())})
    Travel        : ${params.travel}
    CIDT          : ${params.cidt}
    Pathogens     : ${params.pathogen ?: 'default (CAMPYLOBACTER,CYCLOSPORA)'}
    States        : ${params.states ?: 'all'}
    Dashboard     : Enabled
    Cores         : ${params.cpus ?: 'default'}
    Chains        : ${params.chains}
    Iterations    : ${params.iterations}
    Adapt Delta   : ${params.adapt_delta}
    Max Treedepth : ${params.max_treedepth}
    Seed          : ${params.seed}
    Output Dir    : ${params.outdir}/${projID}
    Nextflow Ver  : ${nextflow.version}
    Starting time : ${new Date()}
    ==============================================
    """

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
        
        // Replace mmwrFile for downstream processes with the preprocessed version
        pathogens = pathogens.map { pathogen, file -> [pathogen, processedFile] }
        
        // Log successful preprocessing
        log.info "Raw data preprocessing complete, proceeding to analysis"
    }

    // Scripts directory path
    scripts_path = "${workflow.projectDir}/bin"

    // Run TRENDY with input data
    (model, summary, ir_outputs, plots, irr_outputs, dashboard_trendy, logs) = TRENDY(
        pathogens,
        censusFileB,
        censusFileP,
        projID,
        scripts_path,
        params.travel,
        params.cidt,
        dashboardTemplate
    )

    // Run dashboard generation after all modeling is complete
    dashboard = GENERATE_DASHBOARD(
        ir_outputs.collect(),
        params.outdir + '/' + projID,
        projID,
        dashboardTemplate,
        dashboardScript
    )

    // Handle workflow completion
    workflow.onComplete {
        def status = workflow.success ? 'COMPLETED' : 'FAILED'
        def dashboardInfo = ""
        if (dashboard) {
            dashboardInfo = "\nDashboard       : Available in output directory"
        }
        log.info """
        ==============================================
        FoodNet Trends Analysis: ${status}
        ==============================================
        Completed at     : ${new Date()}
        Duration         : ${workflow.duration}
        Success          : ${workflow.success}
        Work directory   : ${workflow.workDir}
        Exit status      : ${workflow.exitStatus}
        Output directory : ${params.outdir}/${projID}${dashboardInfo}
        ==============================================
        """
    }
}

// Helper function to format file sizes
def formatSize(size) {
    if (size < 1024) return "${size} B"
    else if (size < 1024*1024) return String.format("%.2f KB", size/1024)
    else if (size < 1024*1024*1024) return String.format("%.2f MB", size/(1024*1024))
    else return String.format("%.2f GB", size/(1024*1024*1024))
}

// Helper function to find metadata file in standard or alternate locations
def findMetadataFile(String path) {
    def mainFile = file(path)
    
    if (mainFile.exists()) {
        return mainFile
    }
    
    // Try to find it in the metadata subdirectory
    def parentDir = file(mainFile.getParent())
    def metadataDir = file("${parentDir}/metadata")
    def altFile = file("${metadataDir}/${mainFile.getName()}")
    
    if (altFile.exists()) {
        log.info "Found metadata file in alternate location: ${altFile}"
        return altFile
    }
    
    // If we get here, the file doesn't exist in either location
    error "Metadata file not found: ${path} (also checked in ${metadataDir})"
}
