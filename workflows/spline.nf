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
        log.info "Analyzing ${pathogenList.size()} pathogens: ${pathogenList.join(', ')}"
        
        // First define mmwrFile to avoid getFileSystem errors
        def mmwrFilePath = params.mmwrFile
        mmwrFile = file(mmwrFilePath, checkIfExists: false)
        
        // Then create pathogens channel 
        pathogens = Channel.fromList(pathogenList.collect { pathogen -> [pathogen, mmwrFile] })
    } else {
        // Default to CAMPYLOBACTER and CYCLOSPORA for testing
        def mmwrFilePath = params.mmwrFile
        mmwrFile = file(mmwrFilePath, checkIfExists: false)
        
        pathogens = Channel.fromList([
            ['CAMPYLOBACTER', mmwrFile], 
            ['CYCLOSPORA', mmwrFile]
        ])
        log.info "No pathogens specified, using defaults: CAMPYLOBACTER, CYCLOSPORA"
    }

    // Validate required parameters
    def missingParams = []
    if (!params.mmwrFile) missingParams << "--mmwrFile"
    // Census files are now optional with warnings

    if (missingParams.size() > 0) {
        error "Missing required parameter(s): ${missingParams.join(', ')}"
    }

    // Input files - handle carefully to avoid getFileSystem errors
    if (!mmwrFile.exists()) {
        error "MMWR file does not exist: ${params.mmwrFile}"
    }

    // Handle empty census file parameters
    def censusFileB
    def censusFileP
    
    // Create placeholder files for empty census parameters
    if (!params.censusFileB || params.censusFileB == "") {
        log.warn "WARNING: Census bacterial file parameter is empty"
        log.warn "Will proceed with empty census bacterial file"
        
        // Create empty placeholder file in workDir
        def placeholderB = "${workflow.workDir}/empty_census_bacterial.csv"
        def placeholderBFile = new File(placeholderB)
        if (!placeholderBFile.exists()) {
            placeholderBFile.text = "state,population\n"
        }
        censusFileB = file(placeholderB)
    } else {
        censusFileB = file(params.censusFileB, checkIfExists: false)
        if (!censusFileB.exists()) {
            log.warn "WARNING: Census bacterial file does not exist: ${params.censusFileB}"
            log.warn "Will proceed with empty census bacterial file"
            
            // Create empty placeholder file in workDir
            def placeholderB = "${workflow.workDir}/empty_census_bacterial.csv"
            def placeholderBFile = new File(placeholderB)
            if (!placeholderBFile.exists()) {
                placeholderBFile.text = "state,population\n"
            }
            censusFileB = file(placeholderB)
        }
    }

    if (!params.censusFileP || params.censusFileP == "") {
        log.warn "WARNING: Census parasitic file parameter is empty"
        log.warn "Will proceed with empty census parasitic file"
        
        // Create empty placeholder file in workDir
        def placeholderP = "${workflow.workDir}/empty_census_parasitic.csv"
        def placeholderPFile = new File(placeholderP)
        if (!placeholderPFile.exists()) {
            placeholderPFile.text = "state,population\n"
        }
        censusFileP = file(placeholderP)
    } else {
        censusFileP = file(params.censusFileP, checkIfExists: false)
        if (!censusFileP.exists()) {
            log.warn "WARNING: Census parasitic file does not exist: ${params.censusFileP}"
            log.warn "Will proceed with empty census parasitic file"
            
            // Create empty placeholder file in workDir
            def placeholderP = "${workflow.workDir}/empty_census_parasitic.csv"
            def placeholderPFile = new File(placeholderP)
            if (!placeholderPFile.exists()) {
                placeholderPFile.text = "state,population\n" 
            }
            censusFileP = file(placeholderP)
        }
    }
    
    // Flag for preprocessed data
    isPreprocessed = params.preprocessed ?: false

    // Check file size to catch obvious issues
    try {
        if (mmwrFile.size() == 0) {
            error "MMWR file is empty: ${params.mmwrFile}"
        }
    } catch (Exception e) {
        log.warn "Could not check MMWR file size: ${e.message}"
    }

    // Only check census files if they exist
    if (censusFileB && censusFileB.exists() && censusFileB.size() == 0) {
        log.warn "Census bacterial file is empty: ${params.censusFileB}"
    }

    if (censusFileP && censusFileP.exists() && censusFileP.size() == 0) {
        log.warn "Census parasitic file is empty: ${params.censusFileP}"
    }
    
    // Get metadata file with alternate path fallback
    def metadataFile = null
    if (params.metadata && params.metadata != "") {
        metadataFile = findMetadataFile(params.metadata)
    } else {
        log.warn "No metadata file provided. Some features may be limited."
    }
    
    // Dashboard template file
    def dashboardTemplateFile = "${workflow.projectDir}/assets/dashboard_template.html"
    def dashboardScriptFile = "${workflow.projectDir}/bin/generate_dashboard.R"

    dashboardTemplate = file(dashboardTemplateFile, checkIfExists: false)
    dashboardScript = file(dashboardScriptFile, checkIfExists: false)

    if (!dashboardTemplate.exists()) {
        log.warn "Dashboard template file not found: ${dashboardTemplateFile}"
    }

    if (!dashboardScript.exists()) {
        log.warn "Dashboard script not found: ${dashboardScriptFile}"
    }

    // Set default projID if not specified
    def projID = params.projID ?: new Date().format('yyyyMMdd_HHmmss')

    // Log pipeline start
    log.info """
    ==============================================
    FoodNet Trends Spline Analysis
    ==============================================
    Project ID    : ${projID}
    MMWR File     : ${params.mmwrFile}
    Preprocessed  : ${isPreprocessed}
    Census Files  : 
      Bacterial   : ${params.censusFileB} ${censusFileB && censusFileB.exists() ? "✓" : "✗"}
      Parasitic   : ${params.censusFileP} ${censusFileP && censusFileP.exists() ? "✓" : "✗"}
    Travel        : ${params.travel}
    CIDT          : ${params.cidt}
    Pathogens     : ${params.pathogen ?: 'default (CAMPYLOBACTER,CYCLOSPORA)'}
    States        : ${params.states ?: 'all'}
    Dashboard     : ${dashboardTemplate.exists() ? "Enabled" : "Template not found"}
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
        try {
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
            
            // Replace mmwrFile for downstream processes with the preprocessed version
            pathogens = pathogens.map { pathogen, file -> [pathogen, processedFile] }
            
            // Log successful preprocessing
            log.info "Raw data preprocessing complete, proceeding to analysis"
        } catch (Exception e) {
            log.error "Error in preprocessing: ${e.message}"
            throw e
        }
    }

    // Scripts directory path
    def scripts_path = "${workflow.projectDir}/bin"

    // Run TRENDY with input data
    try {
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
    } catch (Exception e) {
        log.error "Error running TRENDY process: ${e.message}"
        throw e
    }

    // Run dashboard generation after all modeling is complete
    try {
        dashboard = GENERATE_DASHBOARD(
            ir_outputs.collect(),
            ".",
            projID,
            dashboardTemplate,
            dashboardScript
        )
    } catch (Exception e) {
        log.error "Error in dashboard generation: ${e.message}"
        log.warn "Continuing without dashboard"
    }

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
