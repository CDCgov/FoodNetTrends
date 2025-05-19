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
    
    Last updated: 2025-05-19
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import modules
include { TRENDY } from '../modules/local/trendy'
include { PREPROCESS } from '../modules/local/preprocess'
include { GENERATE_DASHBOARD } from '../modules/local/generate_dashboard'

workflow SPLINE {
    // Log workflow version at startup
    log.info "Running FoodNet Trends Spline Analysis Workflow v1.0"
    
    // Validate required parameters
    if (!params.mmwrFile) {
        error "Missing required parameter: --mmwrFile must be specified"
    }
    
    // Read input files - handle carefully
    mmwrFile = file(params.mmwrFile, checkIfExists: true)
    
    // Use simple string handling for census files
    censusFileB = params.censusFileB ? file(params.censusFileB, checkIfExists: false) : file("NO_FILE")
    censusFileP = params.censusFileP ? file(params.censusFileP, checkIfExists: false) : file("NO_FILE")
    
    // Dashboard templates
    dashboardTemplate = file("${workflow.projectDir}/assets/dashboard_template.html", checkIfExists: false)
    dashboardScript = file("${workflow.projectDir}/bin/generate_dashboard.R", checkIfExists: false)
    
    // Define pathogen list from parameter
    def pathogenList = params.pathogen ? params.pathogen.tokenize(',') : ['CAMPYLOBACTER', 'CYCLOSPORA']
    log.info "Analyzing ${pathogenList.size()} pathogens: ${pathogenList.join(', ')}"
    
    // Create pathogen channel
    pathogens = Channel.fromList(pathogenList)
        .map { pathogen -> tuple(pathogen, mmwrFile) }
    
    // Set default projID if not specified
    projID = params.projID ?: new Date().format('yyyyMMdd_HHmmss')
    
    // Flag for preprocessed data
    isPreprocessed = params.preprocessed ?: false
    
    // Scripts directory path
    scripts_path = file("${workflow.projectDir}/bin", checkIfExists: true)
    
    // Log pipeline start
    log.info """
    ==============================================
    FoodNet Trends Spline Analysis
    ==============================================
    Project ID    : ${projID}
    MMWR File     : ${params.mmwrFile}
    Preprocessed  : ${isPreprocessed}
    Census Files  : 
      Bacterial   : ${params.censusFileB ?: 'Not provided'}
      Parasitic   : ${params.censusFileP ?: 'Not provided'}
    Travel        : ${params.travel}
    CIDT          : ${params.cidt}
    Pathogens     : ${pathogenList.join(', ')}
    States        : ${params.states ?: 'all'}
    Output Dir    : ${params.outdir}/${projID}
    Nextflow Ver  : ${nextflow.version}
    Starting time : ${new Date()}
    ==============================================
    """
    
    // Preprocessing step if needed
    if (!isPreprocessed) {
        log.info "Preprocessing raw data files"
        
        PREPROCESS(
            mmwrFile,
            censusFileB,
            censusFileP,
            projID,
            true  // Generate metadata
        )
        
        // Use the preprocessed file for downstream analysis
        processedFile = PREPROCESS.out.cleanedData.first()
        metadataFromProcess = PREPROCESS.out.metadata.first()
        
        // Create new channel for downstream processes
        pathogens = Channel.fromList(pathogenList)
            .map { pathogen -> tuple(pathogen, processedFile) }
        
        log.info "Raw data preprocessing complete, proceeding to analysis"
    }
    
    // Run TRENDY with input data
    TRENDY(
        pathogens,
        censusFileB,
        censusFileP,
        projID,
        scripts_path,
        params.travel,
        params.cidt,
        dashboardTemplate
    )
    
    // Extract outputs for downstream use
    model = TRENDY.out.model
    summary = TRENDY.out.summary
    ir_outputs = TRENDY.out.ir_outputs
    plots = TRENDY.out.plots
    irr_outputs = TRENDY.out.irr_outputs
    dashboard_trendy = TRENDY.out.dashboard
    logs = TRENDY.out.logs
    
    // Run dashboard generation after all modeling is complete
    if (params.enable_dashboard) {
        GENERATE_DASHBOARD(
            ir_outputs.collect(),
            ".",
            projID,
            dashboardTemplate,
            dashboardScript
        )
        dashboard = GENERATE_DASHBOARD.out.dashboard
    }
    
    // Handle workflow completion
    workflow.onComplete {
        log.info """
        ==============================================
        FoodNet Trends Analysis: ${workflow.success ? 'COMPLETED' : 'FAILED'}
        ==============================================
        Completed at     : ${new Date()}
        Duration         : ${workflow.duration}
        Success          : ${workflow.success}
        Work directory   : ${workflow.workDir}
        Exit status      : ${workflow.exitStatus}
        Output directory : ${params.outdir}/${projID}
        ==============================================
        """
    }
}

// Helper function to find metadata file in standard or alternate locations
def findMetadataFile(String path) {
    // First try the exact path as provided
    def mainFilePath = path
    def mainFile = file(mainFilePath, checkIfExists: false)
    
    if (mainFile.exists()) {
        log.info "Using metadata file: ${mainFilePath}"
        return mainFile
    }
    
    // Try to find it in the metadata subdirectory by constructing the path manually
    // Use string operations instead of File objects to avoid getFileSystem errors
    def lastSlash = mainFilePath.lastIndexOf('/')
    if (lastSlash == -1) {
        lastSlash = mainFilePath.lastIndexOf('\\')
    }
    
    def mainFileName = lastSlash > -1 ? mainFilePath.substring(lastSlash + 1) : mainFilePath
    def parentDir = lastSlash > -1 ? mainFilePath.substring(0, lastSlash) : "."
    def altFilePath = "${parentDir}/metadata/${mainFileName}"
    def altFile = file(altFilePath, checkIfExists: false)
    
    if (altFile.exists()) {
        log.info "Found metadata file in alternate location: ${altFilePath}"
        return altFile
    }
    
    // If we get here, the file doesn't exist in either location
    log.warn "Metadata file not found: ${path} (also checked in ${altFilePath})"
    // Return an empty string or null to indicate not found
    return ""
}
