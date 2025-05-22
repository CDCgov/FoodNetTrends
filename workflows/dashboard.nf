/*
 * ==================================================================
 * FoodNetTrends v1.0.0-rc.1 - Dashboard Generation Workflow
 * ==================================================================
 *
 * Purpose:
 *   Standalone workflow for generating interactive HTML dashboards
 *   from completed analysis results. Can be executed independently
 *   after the main modeling pipeline has finished.
 *
 * Workflow Steps:
 *   1. Locate and validate analysis result files
 *   2. Aggregate data from multiple pathogen analyses
 *   3. Generate interactive dashboard with embedded visualizations
 *
 * Use Cases:
 *   - Re-generating dashboards with updated styling
 *   - Creating dashboards from archived analysis results
 *   - Dashboard customization and branding
 *
 * Last updated: 2025-05-22
 * ==================================================================
 */

// Import required modules
include { GENERATE_DASHBOARD } from '../modules/local/dashboard'

// Define the main workflow
workflow DASHBOARD_ONLY {
    main:
    // Set project ID
    def projID = params.projID ?: new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
    
    // Set output directory
    def dashboardDir = "${params.outdir}/${projID}"
    
    log.info "Starting standalone dashboard generation"
    log.info "Project ID: ${projID}"
    log.info "Output directory: ${dashboardDir}"
    
    // Make sure the output directory exists
    new File(dashboardDir).mkdirs()
    
    // Find result files
    def resultFiles = file("${dashboardDir}/*_IRCatch.csv")
    
    // If no files found, check for a different pattern
    if (resultFiles.size() == 0) {
        resultFiles = file("${dashboardDir}/IRCatch_*.csv")
    }
    
    // Count files found
    log.info "Found ${resultFiles.size()} result files"
    
    if (resultFiles.size() == 0) {
        log.error "No result files found in ${dashboardDir}"
        log.error "Dashboard generation cannot proceed without result files"
        log.error "IMPORTANT: Make sure to run this AFTER your analysis is complete!"
        exit 1
    }
    
    // Verify job completion before proceeding
    def jobFiles = file("${dashboardDir}/.command.log")
    def trendy_logs = file("${dashboardDir}/*.log")
    
    // Check for Job Completeness by looking for sentinel files
    log.info "Verifying analysis completion..."
    if (trendy_logs.size() > 0) {
        log.info "Found ${trendy_logs.size()} log files, indicating analysis has run"
    } else {
        log.warn "No log files found. Proceeding anyway, but analysis may not be complete."
    }
    
    // Find dashboard template
    def dashboardTemplateVal = file("${workflow.projectDir}/assets/dashboard_template.html", checkIfExists: false)
    if (!dashboardTemplateVal.exists()) {
        log.warn "Dashboard template not found: ${dashboardTemplateVal}"
        // Try alternate locations
        def altTemplates = [
            file("${workflow.launchDir}/assets/dashboard_template.html", checkIfExists: false),
            file("${workflow.projectDir}/templates/dashboard_template.html", checkIfExists: false)
        ]
        
        for (alt in altTemplates) {
            if (alt.exists()) {
                log.info "Using alternate dashboard template: ${alt}"
                dashboardTemplateVal = alt
                break
            }
        }
    }
    
    // Find dashboard script
    def dashboardScriptVal = file("${workflow.projectDir}/bin/dashboard.R", checkIfExists: false)
    if (!dashboardScriptVal.exists()) {
        log.warn "Dashboard script not found: ${dashboardScriptVal}"
        def altScripts = [
            file("${workflow.launchDir}/bin/dashboard.R", checkIfExists: false)
        ]
        
        for (alt in altScripts) {
            if (alt.exists()) {
                log.info "Using alternate dashboard script: ${alt}"
                dashboardScriptVal = alt
                break
            }
        }
    }
    
    // Create channel from list of result files
    def resultChannel = Channel.fromList(resultFiles)
    
    // Run dashboard generation with proper variable names that match the module definition
    GENERATE_DASHBOARD(
        resultChannel.collect(),  // This becomes ir_outputs in the module
        dashboardDir,             // This becomes resultDir in the module
        projID,                   // This becomes projID in the module
        dashboardTemplateVal,     // This becomes dashboardTemplate in the module
        dashboardScriptVal        // This becomes dashboardScript in the module
    )
}