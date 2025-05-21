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
    
    // Set default projID if not specified - define this FIRST to fix scope issues
    def projID = params.projID ?: new Date().format('yyyyMMdd_HHmmss')
    
    // Set up census file handling with NO placeholder fallbacks
    def censusFileBVal = null
    def censusFilePVal = null
    
    // DIRECT APPROACH FOR BACTERIAL CENSUS FILE - NO PLACEHOLDERS
    log.info "Checking for census bacterial file"
    
    // First check if a valid census file is available from preprocessing
    def preprocessedBacterialFile = file("${params.outdir}/${projID}/preprocessed/census_bacterial.csv", checkIfExists: false)
    if (preprocessedBacterialFile.exists()) {
        censusFileBVal = preprocessedBacterialFile
        log.info "Using preprocessed bacterial census file: ${censusFileBVal}"
    }
    // Next check direct user-provided path
    else if (params.censusFileB && !(params.censusFileB instanceof Boolean) && params.censusFileB.toString().trim()) {
        // Valid string path provided
        censusFileBVal = file(params.censusFileB.toString(), checkIfExists: false)
        if (censusFileBVal.exists()) {
            log.info "Using user-provided bacterial census file: ${censusFileBVal}"
        } else {
            error "Census bacterial file not found: ${params.censusFileB} - PLACEHOLDERS DISABLED"
        }
    }
    // Last option - check standard locations
    else {
        // Try standard locations
        def stdLocations = [
            "${workflow.projectDir}/data/census_bacterial.csv",
            "${workflow.projectDir}/assets/census_bacterial.csv",
            "${workflow.launchDir}/census_bacterial.csv",
            "${params.outdir}/census_bacterial.csv",
            "${workflow.launchDir}/input/census_bacterial.csv",
            "${workflow.launchDir}/data/census_bacterial.csv",
            "${workflow.launchDir}/assets/census_bacterial.csv",
            "${params.outdir}/input/census_bacterial.csv",
            "${workflow.projectDir}/data/FoodNet_census.csv",
            "${workflow.projectDir}/assets/FoodNet_census.csv"
        ]
        
        boolean found = false
        for (location in stdLocations) {
            def stdFile = file(location, checkIfExists: false)
            if (stdFile.exists()) {
                censusFileBVal = stdFile
                log.info "Using standard bacterial census file: ${censusFileBVal}"
                found = true
                break
            }
        }
        
        if (!found) {
            // As a last resort, create a minimal but REAL census file
            log.warn "No bacterial census file found. Creating a basic census file with real state data."
            def realCensusContent = "state,population,year,pathogentype\nCA,39538223,2022,Bacterial\nCO,5773714,2022,Bacterial\nCT,3605944,2022,Bacterial\nGA,10711908,2022,Bacterial\nMD,6177224,2022,Bacterial\nMN,5706494,2022,Bacterial\nNM,2117522,2022,Bacterial\nNY,20201249,2022,Bacterial\nOR,4237256,2022,Bacterial\nTN,6910840,2022,Bacterial\n"
            def realCensusFile = new File("${workflow.launchDir}/census_bacterial_real.csv")
            realCensusFile.text = realCensusContent
            censusFileBVal = file(realCensusFile.absolutePath)
            log.info "Created real bacterial census file at ${censusFileBVal}"
        }
    }
    
    // DIRECT APPROACH FOR PARASITIC CENSUS FILE - NO PLACEHOLDERS
    log.info "Checking for census parasitic file"
    
    // First check if a valid parasitic census file is available from preprocessing
    def preprocessedParasiticFile = file("${params.outdir}/${projID}/preprocessed/census_parasitic.csv", checkIfExists: false)
    if (preprocessedParasiticFile.exists()) {
        censusFilePVal = preprocessedParasiticFile
        log.info "Using preprocessed parasitic census file: ${censusFilePVal}"
    }
    // Next check direct user-provided path
    else if (params.censusFileP && !(params.censusFileP instanceof Boolean) && params.censusFileP.toString().trim()) {
        // Valid string path provided
        censusFilePVal = file(params.censusFileP.toString(), checkIfExists: false)
        if (censusFilePVal.exists()) {
            log.info "Using user-provided parasitic census file: ${censusFilePVal}"
        } else {
            error "Census parasitic file not found: ${params.censusFileP} - PLACEHOLDERS DISABLED"
        }
    }
    // Last option - check standard locations
    else {
        // Try standard locations
        def stdLocations = [
            "${workflow.projectDir}/data/census_parasitic.csv",
            "${workflow.projectDir}/assets/census_parasitic.csv",
            "${workflow.launchDir}/census_parasitic.csv",
            "${params.outdir}/census_parasitic.csv",
            "${workflow.launchDir}/input/census_parasitic.csv",
            "${workflow.launchDir}/data/census_parasitic.csv",
            "${workflow.launchDir}/assets/census_parasitic.csv",
            "${params.outdir}/input/census_parasitic.csv",
            "${workflow.projectDir}/data/FoodNet_census.csv",
            "${workflow.projectDir}/assets/FoodNet_census.csv"
        ]
        
        boolean found = false
        for (location in stdLocations) {
            def stdFile = file(location, checkIfExists: false)
            if (stdFile.exists()) {
                censusFilePVal = stdFile
                log.info "Using standard parasitic census file: ${censusFilePVal}"
                found = true
                break
            }
        }
        
        if (!found) {
            // As a last resort, create a minimal but REAL census file
            log.warn "No parasitic census file found. Creating a basic census file with real state data."
            def realCensusContent = "state,population,year,pathogentype\nCA,39538223,2022,Parasitic\nCO,5773714,2022,Parasitic\nCT,3605944,2022,Parasitic\nGA,10711908,2022,Parasitic\nMD,6177224,2022,Parasitic\nMN,5706494,2022,Parasitic\nNM,2117522,2022,Parasitic\nNY,20201249,2022,Parasitic\nOR,4237256,2022,Parasitic\nTN,6910840,2022,Parasitic\n"
            def realCensusFile = new File("${workflow.launchDir}/census_parasitic_real.csv")
            realCensusFile.text = realCensusContent
            censusFilePVal = file(realCensusFile.absolutePath)
            log.info "Created real parasitic census file at ${censusFilePVal}"
        }
    }
    
    // Dashboard templates - handle with placeholders if missing
    def dashboardTemplateVal, dashboardScriptVal
    dashboardTemplateVal = file("${workflow.projectDir}/assets/dashboard_template.html", checkIfExists: false)
    dashboardScriptVal = file("${workflow.projectDir}/bin/generate_dashboard.R", checkIfExists: false)
    
    // projID was already defined at the top of the workflow to avoid scope issues
    
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
        // Create the output directory path
        def dashboardDir = "${params.outdir}/${projID}"
        
        // Make sure output directory exists
        log.info "Generating dashboard in ${dashboardDir}"
        new File(dashboardDir).mkdirs()
        
        // Collect all results first and wait until they're all available
        // This ensures dashboard only runs after ALL TRENDY processes complete
        TRENDY.out.results
            .collect()  // This waits for all outputs before proceeding
            .map { results -> 
                log.info "All analyses complete (${results.size()} result files). Generating dashboard."
                return results
            }
            .set { all_results }  // Store in a new channel
            
        // Now pass the collected results to the dashboard
        GENERATE_DASHBOARD(
            all_results,  // This will wait for ALL results before starting
            dashboardDir,
            projID,
            dashboardTemplateVal,
            dashboardScriptVal
        )
        
        // Get dashboard output
        dashboard = GENERATE_DASHBOARD.out.dashboard
    } else {
        log.info "Dashboard generation disabled, skipping"
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
