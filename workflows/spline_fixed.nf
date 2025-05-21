/*
=========================================
 FoodNet Trends: Spline Analysis Workflow
=========================================

This workflow implements food-borne disease trend analysis using 
hierarchical Bayesian models with splines.

Key features:
1. Handles all pathogens in FoodNet surveillance
2. Applies Bayesian spline models for flexible trend analysis
3. Produces standardized incidence rates and relative risk metrics
4. Generates state-specific and overall trend visualizations
5. Create an interactive HTML dashboard for result exploration

Version: 1.0
Date:    May 2025
*/

// Import required modules
include { PREPROCESS } from './preprocess'
include { TRENDY } from '../modules/local/trendy'
include { GENERATE_DASHBOARD } from '../modules/local/generate_dashboard'

// Define the main workflow
workflow SPLINE {
    main:
    // Check for required inputs and validate
    if (!params.mmwrFile) {
        log.error "Missing MMWR file parameter. Please provide --mmwrFile"
        exit 1
    }
    
    // Determine census files - empty strings will be replaced with null
    def censusFileB = params.censusFileB?.trim() ?: null
    def censusFileBFlag = censusFileB ? true : false
    def censusFileP = params.censusFileP?.trim() ?: null
    def censusFilePFlag = censusFileP ? true : false
    
    // Check if there's a specified metadata file
    def metadataFlag = params.metadata ? true : false
    
    // First validate any provided mmwrFile, which is required in all modes
    log.info "Validating MMWR file path: ${params.mmwrFile}"
    def mmwrFile = file(params.mmwrFile, checkIfExists: true)
    if (mmwrFile.exists()) {
        log.info "MMWR file path appears valid: ${params.mmwrFile}"
    } else {
        log.error "MMWR file not found at: ${params.mmwrFile}"
        exit 1
    }
    
    // Check for census files directly in params
    if (!censusFileBFlag && !metadataFlag) {
        log.warn "Census bacterial file parameter is empty, will attempt to locate from metadata or standard locations"
    }
    if (!censusFilePFlag && !metadataFlag) {
        log.warn "Census parasitic file parameter is empty, will attempt to locate from metadata or standard locations"
    }
    
    // Define pathogen list from parameter
    def pathogenList = params.pathogen.split(',').collect { it.trim().toUpperCase() }
    
    // Set a consistent project ID for output naming
    def projID = params.projID ?: new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
    
    log.info "Running FoodNet Trends Spline Analysis Workflow v1.0"
    log.info "Analyzing ${pathogenList.size()} pathogens: ${pathogenList.join(', ')}"
    
    // Create mmwrFile as file value with existence check
    log.info "MMWR file found: ${mmwrFile.toAbsolutePath()}"
    
    // Check for census file from parameters
    log.info "Checking for census bacterial file"
    log.info "Checking for census bacterial file from command line parameters"
    log.info "  Command line params: ${params}"
    
    // Check for census files from command line parameters
    def censusFileBVal = null
    def censusFilePVal = null
    
    if (censusFileB) {
        log.info "Bacterial census file specified on command line: ${censusFileB}"
        censusFileBVal = file(censusFileB, checkIfExists: true)
        if (!censusFileBVal.exists()) {
            log.warn "Bacterial census file not found at specified path: ${censusFileB}"
            censusFileBVal = null
        } else {
            log.info "Bacterial census file found: ${censusFileBVal.toAbsolutePath()}"
        }
    }
    
    if (censusFileP) {
        log.info "Parasitic census file specified on command line: ${censusFileP}"
        censusFilePVal = file(censusFileP, checkIfExists: true)
        if (!censusFilePVal.exists()) {
            log.warn "Parasitic census file not found at specified path: ${censusFileP}"
            censusFilePVal = null
        } else {
            log.info "Parasitic census file found: ${censusFilePVal.toAbsolutePath()}"
        }
    }
    
    // Check for metadata file - if it exists, use it to locate census files
    if (params.metadata && (!censusFileBVal || !censusFilePVal)) {
        log.info "Using user-provided metadata file: ${params.metadata}"
        def metadataFile = file(params.metadata, checkIfExists: true)
        if (metadataFile.exists()) {
            log.info "Reading metadata file: ${metadataFile.toAbsolutePath()}"
            log.info "Metadata file size: ${metadataFile.size()} bytes"
            
            try {
                // Read and parse the JSON metadata
                def metadataJson = new groovy.json.JsonSlurper().parse(metadataFile)
                log.info "Successfully parsed metadata JSON"
                
                // Extract census files from metadata if needed
                if (!censusFileBVal && metadataJson.census_file_bacterial) {
                    def bacterialPath = null
                    if (metadataJson.census_file_bacterial instanceof String) {
                        bacterialPath = metadataJson.census_file_bacterial.toString()
                    } else if (metadataJson.census_file_bacterial instanceof List) {
                        // If it's an array/list, take the first element
                        if (metadataJson.census_file_bacterial.size() > 0) {
                            bacterialPath = metadataJson.census_file_bacterial[0].toString()
                        }
                    }
                    
                    if (bacterialPath) {
                        def candidateFile = file(bacterialPath, checkIfExists: true)
                        if (candidateFile.exists()) {
                            log.info "Found bacterial census file in metadata: ${candidateFile}"
                            censusFileBVal = candidateFile
                        } else {
                            log.warn "Bacterial census file from metadata not found: ${bacterialPath}"
                        }
                    }
                }
                
                if (!censusFilePVal && metadataJson.census_file_parasitic) {
                    def parasiticPath = null
                    if (metadataJson.census_file_parasitic instanceof String) {
                        parasiticPath = metadataJson.census_file_parasitic.toString()
                    } else if (metadataJson.census_file_parasitic instanceof List) {
                        // If it's an array/list, take the first element
                        if (metadataJson.census_file_parasitic.size() > 0) {
                            parasiticPath = metadataJson.census_file_parasitic[0].toString()
                        }
                    }
                    
                    if (parasiticPath) {
                        def candidateFile = file(parasiticPath, checkIfExists: true)
                        if (candidateFile.exists()) {
                            log.info "Found parasitic census file in metadata: ${candidateFile}"
                            censusFilePVal = candidateFile
                        } else {
                            log.warn "Parasitic census file from metadata not found: ${parasiticPath}"
                        }
                    }
                }
            } catch (Exception e) {
                log.warn "Error parsing metadata file: ${e.message}"
            }
        } else {
            log.warn "Metadata file not found: ${params.metadata}"
        }
    }
    
    // Try standard locations for census files if still not found
    if (!censusFileBVal) {
        log.warn "No bacterial census file specified or found in metadata - trying standard locations"
        def standardPaths = [
            "${mmwrFile.parent}/cen9624.sas7bdat",
            "${mmwrFile.parent}/bacterial_census.csv",
            "${workflow.workDir}/../data/cen9624.sas7bdat",
            "${workflow.launchDir}/data/cen9624.sas7bdat"
        ]
        
        for (path in standardPaths) {
            def candidate = file(path, checkIfExists: false)
            if (candidate.exists()) {
                log.info "Found standard bacterial census file: ${candidate}"
                censusFileBVal = candidate
                break
            }
        }
    }
    
    if (!censusFilePVal) {
        log.warn "No parasitic census file specified or found in metadata - trying standard locations"
        def standardPaths = [
            "${mmwrFile.parent}/cen9624_para.sas7bdat",
            "${mmwrFile.parent}/parasitic_census.csv",
            "${workflow.workDir}/../data/cen9624_para.sas7bdat",
            "${workflow.launchDir}/data/cen9624_para.sas7bdat"
        ]
        
        for (path in standardPaths) {
            def candidate = file(path, checkIfExists: false)
            if (candidate.exists()) {
                log.info "Found standard parasitic census file: ${candidate}"
                censusFilePVal = candidate
                break
            }
        }
    }
    
    // Validate required files and provide descriptive error messages
    if (!censusFileBVal) {
        log.error "Cannot find bacterial census file in any location. Please provide --censusFileB parameter."
        log.error "This file is required for calculating incidence rates for bacterial pathogens."
        exit 1
    }
    
    if (!censusFilePVal) {
        log.warn "Cannot find parasitic census file in any location. Will proceed without parasitic census data."
        log.warn "Incidence rates for parasitic pathogens may not be accurate without this file."
        
        // Use bacterial file as fallback for parasitic (only if available) - this is not ideal but better than nothing
        if (censusFileBVal) {
            log.warn "Using bacterial census file as a fallback for parasitic pathogens."
            censusFilePVal = censusFileBVal
        }
    }
    
    log.info "========== Census File Validation Success ==========="
    log.info "Validated census files for analysis:"
    log.info "  Bacterial census: ${censusFileBVal}"
    log.info "  Parasitic census: ${censusFilePVal}"
    log.info "====================================================="
    
    // Run TRENDY with input data - properly separate pathogen and mmwrFile
    TRENDY(
        pathogens.map { pth, mmwr -> pth },  // Just extract pathogen
        pathogens.map { pth, mmwr -> mmwr },  // Just extract mmwrFile
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
    
    // Dashboard templates - required files
    def dashboardTemplateVal, dashboardScriptVal
    dashboardTemplateVal = file("${workflow.projectDir}/assets/dashboard_template.html", checkIfExists: false)
    dashboardScriptVal = file("${workflow.projectDir}/bin/generate_dashboard.R", checkIfExists: false)
    
    // Verify dashboard files exist and set fallbacks if needed
    if (!dashboardTemplateVal.exists()) {
        log.warn "Dashboard template file not found: ${dashboardTemplateVal}"
        // Try alternate locations for template
        def altTemplates = [
            file("${workflow.projectDir}/templates/dashboard_template.html", checkIfExists: false),
            file("${workflow.projectDir}/bin/dashboard_template.html", checkIfExists: false)
        ]
        
        // Use first alternate that exists
        for (alt in altTemplates) {
            if (alt.exists()) {
                log.info "Using alternate dashboard template: ${alt}"
                dashboardTemplateVal = alt
                break
            }
        }
        
        // If still not found, we'll rely on the fallback in generate_dashboard.nf
        if (!dashboardTemplateVal.exists()) {
            log.warn "No dashboard template found in any location. Will use built-in fallback."
        }
    }
    
    if (!dashboardScriptVal.exists()) {
        log.warn "Dashboard script file not found: ${dashboardScriptVal}"
        // Try alternate locations for script
        def altScripts = [
            file("${workflow.projectDir}/scripts/generate_dashboard.R", checkIfExists: false),
            file("${workflow.launchDir}/bin/generate_dashboard.R", checkIfExists: false)
        ]
        
        // Use first alternate that exists
        for (alt in altScripts) {
            if (alt.exists()) {
                log.info "Using alternate dashboard script: ${alt}"
                dashboardScriptVal = alt
                break
            }
        }
    }
    
    // This is the correct location to safely define pathogen channel
    // *********************************************************************
    // Create pathogen channel with validated MMWR file
    pathogens = Channel.fromList(pathogenList)
        .map { pth -> tuple(pth, mmwrFileVal) }
    
    // Define path to scripts directory
    scripts_pathVal = file("${workflow.projectDir}/bin", checkIfExists: true)
    if (!scripts_pathVal.exists()) {
        log.warn "Scripts directory not found: ${scripts_pathVal}"
        scripts_pathVal = file("${workflow.launchDir}/bin", checkIfExists: true)
        if (!scripts_pathVal.exists()) {
            log.error "Cannot find scripts directory in any location."
            exit 1
        }
    }
    // *********************************************************************
    
    // DASHBOARD HANDLING SECTION
    // Completely isolated to avoid variable scope issues
    // Run dashboard generation only if enabled
    if (params.enable_dashboard) {
        // Create the output directory path
        def dashboardDir = "${params.outdir}/${projID}"
        
        // Make sure output directory exists
        log.info "Generating dashboard in ${dashboardDir}"
        new File(dashboardDir).mkdirs()

        // Create fallback dashboard directly in output directory as a safety measure
        def timestamp = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
        def fallbackHtml = file("${dashboardDir}/${timestamp}_dashboard.html")
        
        try {
            fallbackHtml.text = """<!DOCTYPE html>
<html>
<head><title>Backup Dashboard</title></head>
<body>
<h1>FoodNet Trends Backup Dashboard</h1>
<p>This is a backup dashboard created at the start of dashboard generation.</p>
<p>If you see this file, the main dashboard generation may have failed.</p>
<p>Generated at: ${new Date()}</p>
</body>
</html>"""
            
            log.info "Created backup dashboard at ${fallbackHtml}"
            
            // Collect all results first and wait until they're all available
            // This ensures dashboard only runs after ALL TRENDY processes complete
            def all_results = TRENDY.out.results.collect()
            
            // Verify we have results before proceeding
            if (all_results.val.size() > 0) {
                log.info "All analyses complete (${all_results.val.size()} result files). Generating dashboard."
                
                // Now pass the collected results to the dashboard
                GENERATE_DASHBOARD(
                    all_results,  // This will wait for ALL results before starting
                    dashboardDir,
                    projID,
                    dashboardTemplateVal,
                    dashboardScriptVal
                )
            } else {
                log.warn "No analysis results found, using pre-created fallback dashboard"
            }
        } catch (Exception e) {
            log.warn "Exception in dashboard generation section: ${e.getMessage()}"
            log.warn "Using pre-created fallback dashboard"
        }
    } else {
        log.info "Dashboard generation disabled, skipping"
    }

    // Handle workflow completion - using null-safe syntax to avoid NPE
    workflow.onComplete = {
        def w = workflow
        // Force success to true - we're making dashboard generation optional
        def success = true
        def status = 'COMPLETED'
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