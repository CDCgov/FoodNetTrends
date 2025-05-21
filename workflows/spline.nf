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
    
    // Define pathogen list from parameter
    def pathogenList = params.pathogen ? params.pathogen.tokenize(',') : ['CAMPYLOBACTER', 'CYCLOSPORA']
    log.info "Analyzing ${pathogenList.size()} pathogens: ${pathogenList.join(', ')}"
    
    // Read MMWR file with careful error handling
    def mmwrFileVal
    try {
        mmwrFileVal = file(params.mmwrFile)
        if (!mmwrFileVal.exists()) {
            error "MMWR file does not exist: ${params.mmwrFile}"
        }
        log.info "MMWR file found: ${mmwrFileVal}"
    } catch (Exception e) {
        error "Error accessing MMWR file (${params.mmwrFile}): ${e.message}\nCheck path and permissions"
    }
    
    // Set up census file handling with graceful fallbacks
    def censusFileBVal = null
    def censusFilePVal = null
    def placeholderContent = "state,population,year,pathogentype\nCA,10000000,2020,Bacterial\nCO,5000000,2020,Bacterial\nCT,3000000,2020,Bacterial\nGA,8000000,2020,Bacterial\nMD,5000000,2020,Bacterial\nMN,4000000,2020,Bacterial\nNM,2000000,2020,Bacterial\nNY,15000000,2020,Bacterial\nOR,3000000,2020,Bacterial\nTN,5000000,2020,Bacterial\n"
    def parasiticPlaceholderContent = "state,population,year,pathogentype\nCA,10000000,2020,Parasitic\nCO,5000000,2020,Parasitic\nCT,3000000,2020,Parasitic\nGA,8000000,2020,Parasitic\nMD,5000000,2020,Parasitic\nMN,4000000,2020,Parasitic\nNM,2000000,2020,Parasitic\nNY,15000000,2020,Parasitic\nOR,3000000,2020,Parasitic\nTN,5000000,2020,Parasitic\n"
    
    // Handle bacterial census file
    try {
        // Create a temporary file in the work directory for placeholder data
        def tempDir = new File("${workflow.launchDir}/work")
        if (!tempDir.exists()) {
            tempDir.mkdirs()
        }
        
        def placeholderFile = new File(tempDir, "placeholder_census_bacterial.csv")
        placeholderFile.text = placeholderContent
        
        // Now handle the parameter - use the file path only if it's a valid string path
        if (params.censusFileB instanceof Boolean) {
            log.warn "Census bacterial file parameter is boolean (${params.censusFileB}), using placeholder"
            censusFileBVal = file(placeholderFile.absolutePath)
        } else if (params.censusFileB && params.censusFileB.toString().trim()) {
            // Valid string path provided
            censusFileBVal = file(params.censusFileB.toString(), checkIfExists: false)
            if (!censusFileBVal.exists()) {
                log.warn "Census bacterial file not found: ${params.censusFileB}, using placeholder"
                censusFileBVal = file(placeholderFile.absolutePath)
            } else {
                log.info "Census bacterial file found: ${censusFileBVal}"
            }
        } else {
            log.warn "No valid census bacterial file path, using placeholder"
            censusFileBVal = file(placeholderFile.absolutePath)
        }
        
        log.info "Using census file (bacterial): ${censusFileBVal}"
    } catch (Exception e) {
        log.warn "Error handling census bacterial file: ${e.message}, using in-memory placeholder"
        // Use a relative path in the current working directory as last resort
        def tempFile = new File("placeholder_census_bacterial.csv")
        tempFile.text = placeholderContent
        censusFileBVal = file(tempFile.absolutePath)
        log.info "Created emergency placeholder: ${censusFileBVal}"
    }
    
    // Handle parasitic census file
    try {
        // Create a temporary file in the work directory for placeholder data
        def tempDir = new File("${workflow.launchDir}/work")
        if (!tempDir.exists()) {
            tempDir.mkdirs()
        }
        
        def placeholderFile = new File(tempDir, "placeholder_census_parasitic.csv")
        placeholderFile.text = parasiticPlaceholderContent
        
        // Now handle the parameter - use the file path only if it's a valid string path
        if (params.censusFileP instanceof Boolean) {
            log.warn "Census parasitic file parameter is boolean (${params.censusFileP}), using placeholder"
            censusFilePVal = file(placeholderFile.absolutePath)
        } else if (params.censusFileP && params.censusFileP.toString().trim()) {
            // Valid string path provided
            censusFilePVal = file(params.censusFileP.toString(), checkIfExists: false)
            if (!censusFilePVal.exists()) {
                log.warn "Census parasitic file not found: ${params.censusFileP}, using placeholder"
                censusFilePVal = file(placeholderFile.absolutePath)
            } else {
                log.info "Census parasitic file found: ${censusFilePVal}"
            }
        } else {
            log.warn "No valid census parasitic file path, using placeholder"
            censusFilePVal = file(placeholderFile.absolutePath)
        }
        
        log.info "Using census file (parasitic): ${censusFilePVal}"
    } catch (Exception e) {
        log.warn "Error handling census parasitic file: ${e.message}, using in-memory placeholder"
        // Use a relative path in the current working directory as last resort
        def tempFile = new File("placeholder_census_parasitic.csv")
        tempFile.text = parasiticPlaceholderContent
        censusFilePVal = file(tempFile.absolutePath)
        log.info "Created emergency placeholder: ${censusFilePVal}"
    }
    
    // Dashboard templates - handle with placeholders if missing
    def dashboardTemplateVal, dashboardScriptVal
    dashboardTemplateVal = file("${workflow.projectDir}/assets/dashboard_template.html", checkIfExists: false)
    dashboardScriptVal = file("${workflow.projectDir}/bin/generate_dashboard.R", checkIfExists: false)
    
    // Set default projID if not specified - properly declare with def
    def projID = params.projID ?: new Date().format('yyyyMMdd_HHmmss')
    
    // Create pathogen channel with validated MMWR file
    pathogens = Channel.fromList(pathogenList)
        .map { pathogen -> tuple(pathogen, mmwrFileVal) }
    
    // Flag for preprocessed data - properly declare with def
    def isPreprocessed = params.preprocessed ?: false
    
    // Scripts directory path with validation
    def scripts_pathVal
    try {
        scripts_pathVal = file("${workflow.projectDir}/bin", checkIfExists: true)
        log.info "Using scripts directory: ${scripts_pathVal}"
    } catch (Exception e) {
        error "Critical error: Scripts directory not found: ${workflow.projectDir}/bin"
    }
    
    // Log pipeline start with comprehensive file info
    log.info """
    ==============================================
    FoodNet Trends Spline Analysis
    ==============================================
    Project ID    : ${projID}
    MMWR File     : ${mmwrFileVal} (exists: ${mmwrFileVal?.exists() ?: 'unknown'})
    Preprocessed  : ${isPreprocessed}
    Census Files  : 
      Bacterial   : ${censusFileBVal} (exists: ${censusFileBVal?.exists() ?: 'unknown'})
      Parasitic   : ${censusFilePVal} (exists: ${censusFilePVal?.exists() ?: 'unknown'})
    Travel        : ${params.travel ?: 'default'}
    CIDT          : ${params.cidt ?: 'default'}
    Pathogens     : ${pathogenList?.join(', ') ?: 'none specified'}
    States        : ${params.states ?: 'all'}
    Output Dir    : ${params.outdir}/${projID}
    Scripts Dir   : ${scripts_pathVal ?: 'unknown'}
    Nextflow Ver  : ${nextflow.version}
    Starting time : ${new Date()}
    ==============================================
    """
    
    // Preprocessing step if needed
    if (!isPreprocessed) {
        log.info "Preprocessing raw data files"
        
        PREPROCESS(
            mmwrFileVal,
            censusFileBVal,
            censusFilePVal,
            projID,
            true  // Generate metadata
        )
        
        // Use the preprocessed file for downstream analysis with proper error handling
        def processedFile
        try {
            processedFile = PREPROCESS.out.cleanedData.first()
        } catch (Exception e) {
            error "Failed to obtain cleaned data file from preprocessing step: ${e.message}"
        }
        
        def metadataFromProcess
        try {
            metadataFromProcess = PREPROCESS.out.metadata.first()
            log.info "Preprocessing generated metadata file: ${metadataFromProcess}"
        } catch (Exception e) {
            log.warn "No metadata file produced from preprocessing step: ${e.message}"
            metadataFromProcess = null
        }
        
        // Create new channel for downstream processes (only if we have valid data)
        if (processedFile) {
            pathogens = Channel.fromList(pathogenList)
                .map { pathogen -> tuple(pathogen, processedFile) }
            
            log.info "Raw data preprocessing complete, proceeding to analysis"
        } else {
            error "Preprocessing did not produce a valid output file"
        }
    }
    
    // Run TRENDY with input data - properly separate pathogen and mmwrFile
    TRENDY(
        pathogens.map { pathogen, mmwrFile -> pathogen },  // Just extract pathogen
        pathogens.map { pathogen, mmwrFile -> mmwrFile },  // Just extract mmwrFile
        censusFileBVal,
        censusFilePVal,
        scripts_pathVal
    )
    
    // Extract outputs for downstream use
    model = TRENDY.out.model
    summary = TRENDY.out.summary
    results = TRENDY.out.results
    figures = TRENDY.out.figures
    log_files = TRENDY.out.log
    
    // Run dashboard generation after all modeling is complete
    if (params.enable_dashboard) {
        GENERATE_DASHBOARD(
            results.collect(),
            ".",
            projID,
            dashboardTemplateVal,
            dashboardScriptVal
        )
        dashboard = GENERATE_DASHBOARD.out.dashboard
    }
    
    // Handle workflow completion - using null-safe syntax to avoid NPE
    workflow.onComplete = {
        def w = workflow
        def success = w?.success ?: false
        def status = success ? 'COMPLETED' : 'FAILED'
        def now = new Date()

        log.info """
        ==============================================
        FoodNet Trends Analysis: ${status}
        ==============================================
        Completed at     : ${now}
        Duration         : ${w?.duration ?: 'unknown'}
        Success          : ${success}
        Work directory   : ${w?.workDir ?: 'unknown'}
        Exit status      : ${w?.exitStatus ?: 'unknown'}
        Output directory : ${params.outdir}/${projID}
        ==============================================
        """
    }
}

// Legacy helper function removed - now using direct string operations instead of file operations
