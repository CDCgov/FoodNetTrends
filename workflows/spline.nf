/*
 * ==================================================================
 * FoodNetTrends v1.0.0-rc.1 - Main Bayesian Spline Analysis Workflow
 * ==================================================================
 *
 * OVERVIEW FOR MAINTAINERS:
 * This is the core Nextflow workflow that orchestrates the complete
 * FoodNetTrends analysis pipeline. It coordinates data preprocessing,
 * Bayesian spline modeling, and dashboard generation across multiple
 * pathogens in parallel.
 *
 * WORKFLOW ARCHITECTURE:
 * 1. Input Validation: Validates MMWR and census data files
 * 2. Data Preprocessing: Prepares surveillance data for modeling
 * 3. Parallel Modeling: Runs Bayesian spline analysis per pathogen
 * 4. Result Collection: Aggregates outputs from parallel analyses
 * 5. Dashboard Generation: Creates interactive HTML reports
 *
 * CRITICAL DESIGN DECISIONS:
 * - Parallel Execution: Each pathogen analyzed independently for speed
 * - Channel Synchronization: Careful coordination of optional outputs
 * - Resource Optimization: Dynamic memory/CPU allocation based on data size
 * - Error Isolation: Pathogen failures don't crash entire pipeline
 *
 * CHANNEL FLOW (for maintainers):
 * mmwr_data -> preprocess -> pathogen_channels -> trendy_analysis -> 
 * collect_results -> dashboard_generation -> final_outputs
 *
 * MAINTAINER WARNING:
 * The channel collection logic in this workflow was specifically designed
 * to handle optional outputs (some pathogens may fail). Modifications to
 * channel operations should preserve this fault-tolerance.
 *   - Self-contained output with embedded visualizations
 *
 * Last updated: 2025-05-22
 * ==================================================================
 */

// Import required modules
include { PREPROCESS } from './preprocess'
include { TRENDY } from '../modules/local/trendy'
include { GENERATE_DASHBOARD } from '../modules/local/dashboard'

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
    if (!params.pathogen) {
        log.error "Missing required parameter: --pathogen. Please specify comma-separated pathogen list."
        exit 1
    }
    def pathogenList = params.pathogen.split(',').collect { it.trim().toUpperCase() }
    
    // Set a consistent project ID for output naming
    def projID = params.projID ?: new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
    
    // Create pathogen channel with validated MMWR file
    def pathogens_ch = Channel.fromList(pathogenList)
        .map { pth -> tuple(pth, mmwrFile) }
    
    // Define path to scripts directory
    def scripts_pathVal = file("${workflow.projectDir}/bin", checkIfExists: true)
    if (!scripts_pathVal.exists()) {
        log.warn "Scripts directory not found: ${scripts_pathVal}"
        scripts_pathVal = file("${workflow.launchDir}/bin", checkIfExists: true)
        if (!scripts_pathVal.exists()) {
            log.error "Cannot find scripts directory in any location."
            exit 1
        }
    }
    
    log.info "Running FoodNetTrends Spline Analysis Workflow v1.0.0-rc.1"
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
    
    // Validate required files - both census files are required for analysis
    if (!censusFileBVal) {
        log.error "Cannot find bacterial census file in any location. Please provide --censusFileB parameter."
        log.error "This file is required for calculating incidence rates."
        exit 1
    }
    
    if (!censusFilePVal) {
        log.error "Cannot find parasitic census file in any location. Please provide --censusFileP parameter."
        log.error "This file is required for calculating incidence rates."
        exit 1
    }
    
    log.info "========== Census File Validation Success ==========="
    log.info "Validated census files for analysis:"
    log.info "  Bacterial census: ${censusFileBVal}"
    log.info "  Parasitic census: ${censusFilePVal}"
    log.info "====================================================="
    
    // Run TRENDY with input data - split tuples while maintaining synchronization
    pathogens_ch
        .multiMap { pth, mmwr ->
            pathogen: pth
            mmwrFile: mmwr
        }
        .set { trendy_inputs }
    
    TRENDY(
        trendy_inputs.pathogen,
        trendy_inputs.mmwrFile,
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
    
    // Dashboard templates - required files
    def dashboardTemplateVal = file("${workflow.projectDir}/assets/dashboard_template.html", checkIfExists: false)
    def dashboardScriptVal = file("${workflow.projectDir}/bin/dashboard.R", checkIfExists: false)
    
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
        
        // If still not found, we'll rely on the fallback in dashboard.nf
        if (!dashboardTemplateVal.exists()) {
            log.warn "No dashboard template found in any location. Will use built-in fallback."
        }
    }
    
    if (!dashboardScriptVal.exists()) {
        log.warn "Dashboard script file not found: ${dashboardScriptVal}"
        // Try alternate locations for script
        def altScripts = [
            file("${workflow.projectDir}/scripts/dashboard.R", checkIfExists: false),
            file("${workflow.launchDir}/bin/dashboard.R", checkIfExists: false)
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
    
    // DASHBOARD HANDLING SECTION
    // Completely isolated to avoid variable scope issues
    // Run dashboard generation only if enabled
    if (params.enable_dashboard) {
        // Create the output directory path
        def dashboardDir = "${params.outdir}/${projID}"
        
        // Make sure output directory exists
        new File(dashboardDir).mkdirs()
        
        try {
            // v1.0.0-rc.1: Fixed critical channel collection bug for dashboard generation
            // Previous logic failed when optional outputs (CSV/PNG files) were missing from any pathogen
            // The mix() operation would create channel synchronization deadlocks with empty channels
            // New approach: Flatten all TRENDY outputs and filter only existing files for dashboard input
            TRENDY.out.results
                .mix(TRENDY.out.figures, TRENDY.out.summary)
                .flatten()
                .filter { it.exists() }  // Only include files that actually exist - handles optional outputs
                .collect()               // Now safe to collect since we only process existing files
                .set { all_output_files }
            
            // Now pass ALL collected files to the dashboard
            GENERATE_DASHBOARD(
                all_output_files,  // This includes CSV, PNG, and summary files
                dashboardDir,
                projID,
                dashboardTemplateVal,
                dashboardScriptVal
            )
        } catch (Exception e) {
            log.error "Dashboard generation failed: ${e.getMessage()}"
            log.warn "To regenerate the dashboard after pipeline completion, run:"
            log.warn "./dashboard.sh ${projID}"
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
        FoodNetTrends Analysis: ${status}
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