/*
=========================================
 FoodNet Trends: Simplified Spline Analysis Workflow
=========================================

This is a simplified version of the FoodNet Trends spline analysis workflow,
designed to avoid variable scope issues and ensure reliable execution.
*/

// Import required modules
include { TRENDY } from '../modules/local/trendy'
include { GENERATE_DASHBOARD } from '../modules/local/generate_dashboard'

// Define the main workflow
workflow SPLINE {
    main:
    // Validate MMWR file
    log.info "Validating MMWR file path: ${params.mmwrFile}"
    def mmwrFileVal = file(params.mmwrFile, checkIfExists: true)
    if (!mmwrFileVal.exists()) {
        log.error "MMWR file not found at: ${params.mmwrFile}"
        exit 1
    }
    log.info "MMWR file found: ${mmwrFileVal.toAbsolutePath()}"

    // Validate census files
    log.info "Checking census files..."
    def censusFileBVal = null
    def censusFilePVal = null
    
    // Check bacterial census
    if (params.censusFileB) {
        log.info "Bacterial census file specified on command line: ${params.censusFileB}"
        censusFileBVal = file(params.censusFileB, checkIfExists: true)
        if (!censusFileBVal.exists()) {
            log.error "Bacterial census file not found at: ${params.censusFileB}"
            exit 1
        }
        log.info "Bacterial census file found: ${censusFileBVal.toAbsolutePath()}"
    } else {
        log.error "No bacterial census file specified (--censusFileB)"
        exit 1
    }
    
    // Check parasitic census
    if (params.censusFileP) {
        log.info "Parasitic census file specified on command line: ${params.censusFileP}"
        censusFilePVal = file(params.censusFileP, checkIfExists: true)
        if (!censusFilePVal.exists()) {
            log.error "Parasitic census file not found at: ${params.censusFileP}"
            exit 1
        }
        log.info "Parasitic census file found: ${censusFilePVal.toAbsolutePath()}"
    } else {
        log.warn "No parasitic census file specified, using bacterial as fallback"
        censusFilePVal = censusFileBVal
    }
    
    log.info "========== Census File Validation Success ==========="
    log.info "Validated census files for analysis:"
    log.info "  Bacterial census: ${censusFileBVal}"
    log.info "  Parasitic census: ${censusFilePVal}"
    log.info "====================================================="
    
    // Find script path
    def scripts_pathVal = file("${workflow.projectDir}/bin", checkIfExists: true)
    if (!scripts_pathVal.exists()) {
        log.warn "Scripts directory not found: ${scripts_pathVal}"
        scripts_pathVal = file("${workflow.launchDir}/bin", checkIfExists: true)
        if (!scripts_pathVal.exists()) {
            log.error "Cannot find scripts directory in any location."
            exit 1
        }
    }
    
    // Get project ID
    def projID = params.projID ?: new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
    log.info "Project ID: ${projID}"
    
    // Create pathogens channel BEFORE using it
    def pathogensList = params.pathogen.split(',').collect { it.trim().toUpperCase() }
    log.info "Analyzing ${pathogensList.size()} pathogens: ${pathogensList.join(', ')}"
    
    def pathogens_ch = Channel.fromList(pathogensList)
        .map { pth -> tuple(pth, mmwrFileVal) }
    
    // Run TRENDY with input data
    TRENDY(
        pathogens_ch.map { pth, mmwr -> pth },  // Just extract pathogen
        pathogens_ch.map { pth, mmwr -> mmwr },  // Just extract mmwrFile
        censusFileBVal,
        censusFilePVal,
        scripts_pathVal
    )
    
    // Extract outputs for downstream use
    def model = TRENDY.out.model
    def summary = TRENDY.out.summary
    def results = TRENDY.out.results
    def figures = TRENDY.out.figures
    def log_files = TRENDY.out.log
    
    // === DASHBOARD GENERATION ===
    if (params.enable_dashboard) {
        // Create the output directory path
        def dashboardDir = "${params.outdir}/${projID}"
        
        // Make sure output directory exists
        log.info "Generating dashboard in ${dashboardDir}"
        new File(dashboardDir).mkdirs()
        
        // Find dashboard template and script
        def dashboardTemplateVal = file("${workflow.projectDir}/assets/dashboard_template.html", checkIfExists: false)
        def dashboardScriptVal = file("${workflow.projectDir}/bin/generate_dashboard.R", checkIfExists: false)
        
        // Check template
        if (!dashboardTemplateVal.exists()) {
            log.warn "Dashboard template file not found: ${dashboardTemplateVal}"
            // Try alternate locations
            if (file("${workflow.launchDir}/assets/dashboard_template.html").exists()) {
                dashboardTemplateVal = file("${workflow.launchDir}/assets/dashboard_template.html")
                log.info "Using alternate template: ${dashboardTemplateVal}"
            }
        }
        
        // Check script
        if (!dashboardScriptVal.exists()) {
            log.warn "Dashboard script file not found: ${dashboardScriptVal}"
            // Try alternate locations
            if (file("${workflow.launchDir}/bin/generate_dashboard.R").exists()) {
                dashboardScriptVal = file("${workflow.launchDir}/bin/generate_dashboard.R")
                log.info "Using alternate script: ${dashboardScriptVal}"
            }
        }
        
        // Create fallback dashboard in case of failure
        log.info "Creating backup dashboard in case of generation failure"
        def timestamp = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
        def fallbackPath = "${dashboardDir}/${projID}_dashboard.html"
        def fallbackHtml = file(fallbackPath)
        try {
            fallbackHtml.text = """<!DOCTYPE html>
<html><head><title>FoodNet Trends Backup Dashboard</title></head>
<body>
<h1>FoodNet Trends Backup Dashboard</h1>
<p>This is a backup dashboard created as a fallback.</p>
<p>Generated at: ${new Date()}</p>
</body></html>"""
        } catch (Exception e) {
            log.warn "Could not create backup dashboard: ${e.message}"
        }
        
        // Collect results and run dashboard generation
        def all_results = TRENDY.out.results.collect()
        
        // Run dashboard generation
        try {
            GENERATE_DASHBOARD(
                all_results,
                dashboardDir,
                projID,
                dashboardTemplateVal,
                dashboardScriptVal
            )
        } catch (Exception e) {
            log.warn "Dashboard generation failed: ${e.message}"
            log.warn "Using backup dashboard"
        }
    } else {
        log.info "Dashboard generation disabled, skipping"
    }
    
    // Final message
    log.info """
    ==============================================
    FoodNet Trends Analysis: COMPLETED
    ==============================================
    Project ID      : ${projID}
    Output directory: ${params.outdir}/${projID}
    ==============================================
    """
}