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
    
    // Set up census file handling with direct parameter priority
    def censusFileBVal = null
    def censusFilePVal = null
    
    // Check for census bacterial file - REQUIRED
    log.info "Checking for census bacterial file"
    
    // Direct command-line parameter for census files takes priority
    if (params.censusFileB && !(params.censusFileB instanceof Boolean) && params.censusFileB.toString().trim()) {
        // Valid string path provided
        try {
            def bPath = params.censusFileB.toString().trim()
            censusFileBVal = file(bPath)
            if (censusFileBVal.exists()) {
                log.info "Using user-provided bacterial census file: ${censusFileBVal}"
            } else {
                log.warn "User-provided census bacterial file not found: ${bPath}"
                censusFileBVal = null
            }
        } catch (Exception e) {
            log.warn "Error processing bacterial census file parameter: ${e.message}"
            censusFileBVal = null
        }
    } else {
        log.warn "No direct bacterial census file parameter provided"
    }
    
    // Second priority: Check metadata file for census file paths
    if (censusFileBVal == null) {
        // First check if user provided a metadata path
        def metadataFile = null
        if (params.metadata && !(params.metadata instanceof Boolean) && params.metadata.toString().trim()) {
            metadataFile = file(params.metadata.toString(), checkIfExists: false)
            log.info "Using user-provided metadata file: ${metadataFile}"
        } else {
            // Otherwise check default location
            metadataFile = file("${params.outdir}/${projID}/preprocessed/${projID}_metadata.json", checkIfExists: false)
        }
        
        if (metadataFile != null && metadataFile.exists()) {
            try {
                log.info "Reading metadata file: ${metadataFile}"
                def metadataContent = metadataFile.text
                log.info "Metadata file size: ${metadataContent.size()} bytes"
                
                try {
                    def slurper = new nextflow.util.JsonSlurper()
                    def metadataJson = slurper.parseText(metadataContent)
                    
                    log.info "Successfully parsed metadata JSON"
                    if (metadataJson.census_file_bacterial) {
                    def bacterialPath = metadataJson.census_file_bacterial
                    def bacterialFile = file(bacterialPath, checkIfExists: false)
                    if (bacterialFile.exists()) {
                        censusFileBVal = bacterialFile
                        log.info "Using bacterial census file from metadata: ${censusFileBVal}"
                    } else {
                        log.warn "Bacterial census file in metadata doesn't exist: ${bacterialPath}"
                    }
                } else {
                    log.warn "Metadata doesn't contain bacterial census file path"
                }
                } catch (Exception e) {
                    log.warn "Error parsing metadata JSON: ${e.message}"
                }
            } catch (Exception e) {
                log.warn "Error reading metadata file: ${e.message}"
            }
        }
    }
    
    // Third priority: Check for standard preprocessed file
    if (censusFileBVal == null) {
        def preprocessedBacterialFile = file("${params.outdir}/${projID}/preprocessed/census_bacterial.csv", checkIfExists: false)
        if (preprocessedBacterialFile.exists()) {
            censusFileBVal = preprocessedBacterialFile
            log.info "Using preprocessed bacterial census file: ${censusFileBVal}"
        }
    }
    // Final attempt: Try standard locations
    else {
        log.info "Searching standard locations for bacterial census file..."
        // Look in all possible standard locations with both CSV and SAS formats
        def stdLocations = [
            // Current project directory locations
            "${workflow.projectDir}/data/census_bacterial.csv",
            "${workflow.projectDir}/assets/census_bacterial.csv",
            "${workflow.projectDir}/census_bacterial.csv",
            "${workflow.projectDir}/data/census_bacterial.sas7bdat",
            "${workflow.projectDir}/assets/census_bacterial.sas7bdat",
            "${workflow.projectDir}/census_bacterial.sas7bdat",
            
            // User-specified output directory
            "${params.outdir}/census_bacterial.csv",
            "${params.outdir}/census_bacterial.sas7bdat",
            "${params.outdir}/data/census_bacterial.csv",
            "${params.outdir}/data/census_bacterial.sas7bdat",
            
            // Launch directory (where nextflow was started)
            "${workflow.launchDir}/census_bacterial.csv",
            "${workflow.launchDir}/census_bacterial.sas7bdat",
            "${workflow.launchDir}/data/census_bacterial.csv",
            "${workflow.launchDir}/data/census_bacterial.sas7bdat",
            "${workflow.launchDir}/input/census_bacterial.csv",
            "${workflow.launchDir}/input/census_bacterial.sas7bdat",
            "${workflow.launchDir}/assets/census_bacterial.csv",
            "${workflow.launchDir}/assets/census_bacterial.sas7bdat",
            
            // Generic census file names
            "${workflow.projectDir}/data/FoodNet_census.csv",
            "${workflow.projectDir}/assets/FoodNet_census.csv",
            "${workflow.launchDir}/data/FoodNet_census.csv",
            "${workflow.launchDir}/FoodNet_census.csv",
            "${params.outdir}/FoodNet_census.csv"
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
            error "No bacterial census file found. Please provide a valid census file using --censusFileB parameter."
        }
    }
    
    // Check for parasitic census file - REQUIRED
    log.info "Checking for census parasitic file"
    
    // Direct command-line parameter for census files takes priority
    if (params.censusFileP && !(params.censusFileP instanceof Boolean) && params.censusFileP.toString().trim()) {
        // Valid string path provided
        try {
            def pPath = params.censusFileP.toString().trim()
            censusFilePVal = file(pPath)
            if (censusFilePVal.exists()) {
                log.info "Using user-provided parasitic census file: ${censusFilePVal}"
            } else {
                log.warn "User-provided census parasitic file not found: ${pPath}"
                censusFilePVal = null
            }
        } catch (Exception e) {
            log.warn "Error processing parasitic census file parameter: ${e.message}"
            censusFilePVal = null
        }
    } else {
        log.warn "No direct parasitic census file parameter provided"
    }
    
    // Second priority: Check metadata file for census file paths
    if (censusFilePVal == null) {
        // First check if user provided a metadata path
        def metadataFile = null
        if (params.metadata && !(params.metadata instanceof Boolean) && params.metadata.toString().trim()) {
            metadataFile = file(params.metadata.toString(), checkIfExists: false)
            log.info "Using user-provided metadata file for parasitic census: ${metadataFile}"
        } else {
            // Otherwise check default location
            metadataFile = file("${params.outdir}/${projID}/preprocessed/${projID}_metadata.json", checkIfExists: false)
        }
        
        if (metadataFile != null && metadataFile.exists()) {
            try {
                log.info "Reading metadata file for parasitic census: ${metadataFile}"
                def metadataContent = metadataFile.text
                log.info "Metadata file size: ${metadataContent.size()} bytes"
                
                try {
                    def slurper = new nextflow.util.JsonSlurper()
                    def metadataJson = slurper.parseText(metadataContent)
                    
                    log.info "Successfully parsed metadata JSON for parasitic census"
                    if (metadataJson.census_file_parasitic) {
                    def parasiticPath = metadataJson.census_file_parasitic
                    def parasiticFile = file(parasiticPath, checkIfExists: false)
                    if (parasiticFile.exists()) {
                        censusFilePVal = parasiticFile
                        log.info "Using parasitic census file from metadata: ${censusFilePVal}"
                    } else {
                        log.warn "Parasitic census file in metadata doesn't exist: ${parasiticPath}"
                    }
                } else {
                    log.warn "Metadata doesn't contain parasitic census file path"
                }
                } catch (Exception e) {
                    log.warn "Error parsing metadata JSON for parasitic census: ${e.message}"
                }
            } catch (Exception e) {
                log.warn "Error reading metadata file for parasitic census: ${e.message}"
            }
        }
    }
    
    // Third priority: Check for standard preprocessed file
    if (censusFilePVal == null) {
        def preprocessedParasiticFile = file("${params.outdir}/${projID}/preprocessed/census_parasitic.csv", checkIfExists: false)
        if (preprocessedParasiticFile.exists()) {
            censusFilePVal = preprocessedParasiticFile
            log.info "Using preprocessed parasitic census file: ${censusFilePVal}"
        }
    }
    // Final attempt: Try standard locations
    else {
        log.info "Searching standard locations for parasitic census file..."
        // Look in all possible standard locations with both CSV and SAS formats
        def stdLocations = [
            // Current project directory locations
            "${workflow.projectDir}/data/census_parasitic.csv",
            "${workflow.projectDir}/assets/census_parasitic.csv",
            "${workflow.projectDir}/census_parasitic.csv",
            "${workflow.projectDir}/data/census_parasitic.sas7bdat",
            "${workflow.projectDir}/assets/census_parasitic.sas7bdat",
            "${workflow.projectDir}/census_parasitic.sas7bdat",
            
            // User-specified output directory
            "${params.outdir}/census_parasitic.csv",
            "${params.outdir}/census_parasitic.sas7bdat",
            "${params.outdir}/data/census_parasitic.csv",
            "${params.outdir}/data/census_parasitic.sas7bdat",
            
            // Launch directory (where nextflow was started)
            "${workflow.launchDir}/census_parasitic.csv",
            "${workflow.launchDir}/census_parasitic.sas7bdat",
            "${workflow.launchDir}/data/census_parasitic.csv",
            "${workflow.launchDir}/data/census_parasitic.sas7bdat",
            "${workflow.launchDir}/input/census_parasitic.csv",
            "${workflow.launchDir}/input/census_parasitic.sas7bdat",
            "${workflow.launchDir}/assets/census_parasitic.csv",
            "${workflow.launchDir}/assets/census_parasitic.sas7bdat",
            
            // Generic census file names - parasitic might use the same as bacterial in some datasets
            "${workflow.projectDir}/data/FoodNet_census_para.csv",
            "${workflow.projectDir}/assets/FoodNet_census_para.csv",
            "${workflow.launchDir}/data/FoodNet_census_para.csv",
            "${workflow.launchDir}/FoodNet_census_para.csv",
            "${params.outdir}/FoodNet_census_para.csv"
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
            error "No parasitic census file found. Please provide a valid census file using --censusFileP parameter."
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
            
            // Store metadata path for future reference
            params.metadata = metadataFromProcess.toString()
            log.info "Setting metadata parameter to: ${params.metadata}"
            
            // If we have metadata from preprocessing, extract census file paths if available
            try {
                log.info "Reading metadata from preprocessing: ${metadataFromProcess}"
                def metadataContent = metadataFromProcess.text
                log.info "Metadata file size: ${metadataContent.size()} bytes"
                
                try {
                    def slurper = new nextflow.util.JsonSlurper()
                    def metadataJson = slurper.parseText(metadataContent)
                    
                    log.info "Successfully parsed metadata JSON from preprocessing"
                    // Check for bacterial census path
                    if (metadataJson.census_file_bacterial) {
                    def bacterialPath = metadataJson.census_file_bacterial
                    def bacterialFile = file(bacterialPath, checkIfExists: false)
                    if (bacterialFile.exists()) {
                        censusFileBVal = bacterialFile
                        log.info "Updated bacterial census file from preprocessing metadata: ${censusFileBVal}"
                    }
                }
                
                // Check for parasitic census path
                if (metadataJson.census_file_parasitic) {
                    def parasiticPath = metadataJson.census_file_parasitic
                    def parasiticFile = file(parasiticPath, checkIfExists: false)
                    if (parasiticFile.exists()) {
                        censusFilePVal = parasiticFile
                        log.info "Updated parasitic census file from preprocessing metadata: ${censusFilePVal}"
                    }
                }
                } catch (Exception e) {
                    log.warn "Error parsing metadata JSON from preprocessing: ${e.message}"
                }
            } catch (Exception e) {
                log.warn "Error reading metadata from preprocessing: ${e.message}"
            }
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
    
    // Validate census files before proceeding and provide helpful error messages
    if (censusFileBVal == null || !censusFileBVal.exists()) {
        log.error """
        ========================================================================
        ERROR: Required bacterial census file not found
        
        Census files are required for accurate rate calculations.
        Please ensure the bacterial census file exists and is specified via one of:
        
        1. Command-line parameter: --censusFileB "/path/to/census_bacterial.csv"
        2. Available in metadata from preprocessing
        3. Located in a standard location (checked multiple paths)
        
        The file may be in CSV or SAS7BDAT format.
        ========================================================================
        """
        error "A valid bacterial census file (censusFileB) is required but was not found. Please check your parameters."
    }
    
    if (censusFilePVal == null || !censusFilePVal.exists()) {
        log.error """
        ========================================================================
        ERROR: Required parasitic census file not found
        
        Census files are required for accurate rate calculations.
        Please ensure the parasitic census file exists and is specified via one of:
        
        1. Command-line parameter: --censusFileP "/path/to/census_parasitic.csv"
        2. Available in metadata from preprocessing
        3. Located in a standard location (checked multiple paths)
        
        The file may be in CSV or SAS7BDAT format.
        ========================================================================
        """
        error "A valid parasitic census file (censusFileP) is required but was not found. Please check your parameters."
    }
    
    log.info "========== Census File Validation Success ==========="
    log.info "Validated census files for analysis:"
    log.info "  Bacterial census: ${censusFileBVal}"
    log.info "  Parasitic census: ${censusFilePVal}"
    log.info "====================================================="
    
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
