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
    def censusFileBVal, censusFilePVal
    try {
        if (params.censusFileB) {
            censusFileBVal = file(params.censusFileB, checkIfExists: false)
            if (!censusFileBVal.exists()) {
                log.warn "Census bacterial file not found: ${params.censusFileB}, will use placeholder"
                censusFileBVal = file("${workflow.projectDir}/work/placeholder_census_bact.csv")
                if (!censusFileBVal.exists()) {
                    // Create minimal placeholder file if it doesn't exist
                    def placeholder = file("${workflow.launchDir}/placeholder_census_bact.csv")
                    placeholder.text = "state,population,year,pathogentype\nCA,10000000,2020,Bacterial\n"
                    censusFileBVal = placeholder
                    log.info "Created census bacterial placeholder: ${censusFileBVal}"
                }
            } else {
                log.info "Census bacterial file found: ${censusFileBVal}"
            }
        } else {
            log.warn "No census bacterial file specified, will use placeholder"
            def placeholder = file("${workflow.launchDir}/placeholder_census_bact.csv")
            placeholder.text = "state,population,year,pathogentype\nCA,10000000,2020,Bacterial\n"
            censusFileBVal = placeholder
            log.info "Created census bacterial placeholder: ${censusFileBVal}"
        }
    } catch (Exception e) {
        log.warn "Error handling census bacterial file: ${e.message}, using placeholder"
        def placeholder = file("${workflow.launchDir}/placeholder_census_bact.csv")
        placeholder.text = "state,population,year,pathogentype\nCA,10000000,2020,Bacterial\n"
        censusFileBVal = placeholder
    }
    
    try {
        if (params.censusFileP) {
            censusFilePVal = file(params.censusFileP, checkIfExists: false)
            if (!censusFilePVal.exists()) {
                log.warn "Census parasitic file not found: ${params.censusFileP}, will use placeholder"
                censusFilePVal = file("${workflow.projectDir}/work/placeholder_census_para.csv")
                if (!censusFilePVal.exists()) {
                    // Create minimal placeholder file if it doesn't exist
                    def placeholder = file("${workflow.launchDir}/placeholder_census_para.csv")
                    placeholder.text = "state,population,year,pathogentype\nCA,10000000,2020,Parasitic\n"
                    censusFilePVal = placeholder
                    log.info "Created census parasitic placeholder: ${censusFilePVal}"
                }
            } else {
                log.info "Census parasitic file found: ${censusFilePVal}"
            }
        } else {
            log.warn "No census parasitic file specified, will use placeholder"
            def placeholder = file("${workflow.launchDir}/placeholder_census_para.csv")
            placeholder.text = "state,population,year,pathogentype\nCA,10000000,2020,Parasitic\n"
            censusFilePVal = placeholder
            log.info "Created census parasitic placeholder: ${censusFilePVal}"
        }
    } catch (Exception e) {
        log.warn "Error handling census parasitic file: ${e.message}, using placeholder"
        def placeholder = file("${workflow.launchDir}/placeholder_census_para.csv")
        placeholder.text = "state,population,year,pathogentype\nCA,10000000,2020,Parasitic\n"
        censusFilePVal = placeholder
    }
    
    // Dashboard templates - handle with placeholders if missing
    def dashboardTemplateVal, dashboardScriptVal
    dashboardTemplateVal = file("${workflow.projectDir}/assets/dashboard_template.html", checkIfExists: false)
    dashboardScriptVal = file("${workflow.projectDir}/bin/generate_dashboard.R", checkIfExists: false)
    
    // Set default projID if not specified
    projID = params.projID ?: new Date().format('yyyyMMdd_HHmmss')
    
    // Create pathogen channel with validated MMWR file
    pathogens = Channel.fromList(pathogenList)
        .map { pathogen -> tuple(pathogen, mmwrFileVal) }
    
    // Flag for preprocessed data
    isPreprocessed = params.preprocessed ?: false
    
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
    MMWR File     : ${mmwrFileVal} (exists: ${mmwrFileVal.exists()})
    Preprocessed  : ${isPreprocessed}
    Census Files  : 
      Bacterial   : ${censusFileBVal} (exists: ${censusFileBVal.exists()})
      Parasitic   : ${censusFilePVal} (exists: ${censusFilePVal.exists()})
    Travel        : ${params.travel}
    CIDT          : ${params.cidt}
    Pathogens     : ${pathogenList.join(', ')}
    States        : ${params.states ?: 'all'}
    Output Dir    : ${params.outdir}/${projID}
    Scripts Dir   : ${scripts_pathVal}
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
        censusFileBVal,
        censusFilePVal,
        projID,
        scripts_pathVal,
        params.travel,
        params.cidt,
        dashboardTemplateVal
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
            dashboardTemplateVal,
            dashboardScriptVal
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

// Legacy helper function removed - now using direct string operations instead of file operations
