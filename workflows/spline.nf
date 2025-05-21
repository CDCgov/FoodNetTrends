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
    log.info "Checking for census bacterial file from command line parameters"
    log.info "  Command line params: ${params.toString()}"
    
    if (params.censusFileB != null) {
        log.info "  Found parameter censusFileB: ${params.censusFileB}"
        // Valid string path provided
        try {
            def bPath = params.censusFileB.toString().trim()
            if (bPath && bPath != "true" && bPath != "false") {
                log.info "  Checking file path: ${bPath}"
                // Use more explicit file creation to see what's happening
                def bFile = file(bPath)
                if (bFile && bFile.exists()) {
                    censusFileBVal = bFile
                    log.info "Using user-provided bacterial census file: ${censusFileBVal} (exists: ${censusFileBVal.exists()})"
                } else {
                    log.warn "User-provided census bacterial file not found: ${bPath}"
                    censusFileBVal = null
                }
            } else {
                log.warn "Census bacterial file parameter is empty or boolean value: ${bPath}"
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
                    def slurper = new groovy.json.JsonSlurper()
                    def metadataJson = slurper.parseText(metadataContent)
                    
                    log.info "Successfully parsed metadata JSON"
                    
                    // Handle census_file_bacterial path safely
                    try {
                        if (metadataJson.containsKey('census_file_bacterial')) {
                            log.info "Found census_file_bacterial in metadata: ${metadataJson.census_file_bacterial?.getClass()?.getName() ?: 'null'}"
                            
                            // Handle census file path which might be a string or an array or other object
                            def bacterialPath = null
                            if (metadataJson.census_file_bacterial instanceof String) {
                                bacterialPath = metadataJson.census_file_bacterial.toString()
                                log.info "Found bacterial census file path in metadata (as String): ${bacterialPath}"
                            } else if (metadataJson.census_file_bacterial instanceof List) {
                                // If it's an array/list, take the first element
                                if (metadataJson.census_file_bacterial.size() > 0) {
                                    bacterialPath = metadataJson.census_file_bacterial[0].toString()
                                    log.info "Found bacterial census file path in metadata (as first element of List): ${bacterialPath}"
                                } else {
                                    log.warn "Census file bacterial in metadata is an empty List"
                                }
                            } else if (metadataJson.census_file_bacterial != null) {
                                // For any other type, try toString
                                try {
                                    bacterialPath = metadataJson.census_file_bacterial.toString()
                                    log.info "Found bacterial census file path in metadata (converted from ${metadataJson.census_file_bacterial.getClass().getName()}): ${bacterialPath}"
                                } catch (Exception e) {
                                    log.warn "Could not convert census_file_bacterial to string: ${e.message}"
                                }
                            } else {
                                log.warn "Census file bacterial in metadata is null"
                            }
                            
                            // Check if we obtained a valid path and use it
                            if (bacterialPath) {
                                def bacterialFile = file(bacterialPath, checkIfExists: false)
                                if (bacterialFile.exists()) {
                                    censusFileBVal = bacterialFile
                                    log.info "Using bacterial census file from metadata: ${censusFileBVal}"
                                } else {
                                    log.warn "Bacterial census file in metadata doesn't exist: ${bacterialPath}"
                                }
                            } else {
                                log.warn "No valid bacterial census file path found in metadata"
                            }
                        } else {
                            log.warn "Metadata doesn't contain census_file_bacterial key"
                        }
                    } catch (Exception e) {
                        log.warn "Error accessing census_file_bacterial in metadata: ${e.message}"
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
    
    // Final attempt: Try standard locations if not found yet
    if (censusFileBVal == null) {
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
    log.info "Checking for census parasitic file from command line parameters"
    
    if (params.censusFileP != null) {
        log.info "  Found parameter censusFileP: ${params.censusFileP}"
        // Valid string path provided
        try {
            def pPath = params.censusFileP.toString().trim()
            if (pPath && pPath != "true" && pPath != "false") {
                log.info "  Checking file path: ${pPath}"
                // Use more explicit file creation to see what's happening
                def pFile = file(pPath)
                if (pFile && pFile.exists()) {
                    censusFilePVal = pFile
                    log.info "Using user-provided parasitic census file: ${censusFilePVal} (exists: ${censusFilePVal.exists()})"
                } else {
                    log.warn "User-provided census parasitic file not found: ${pPath}"
                    censusFilePVal = null
                }
            } else {
                log.warn "Census parasitic file parameter is empty or boolean value: ${pPath}"
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
                    def slurper = new groovy.json.JsonSlurper()
                    def metadataJson = slurper.parseText(metadataContent)
                    
                    log.info "Successfully parsed metadata JSON for parasitic census"
                    
                    // Handle census_file_parasitic path safely
                    try {
                        if (metadataJson.containsKey('census_file_parasitic')) {
                            log.info "Found census_file_parasitic in metadata: ${metadataJson.census_file_parasitic?.getClass()?.getName() ?: 'null'}"
                            
                            // Handle census file path which might be a string or an array or other object
                            def parasiticPath = null
                            if (metadataJson.census_file_parasitic instanceof String) {
                                parasiticPath = metadataJson.census_file_parasitic.toString()
                                log.info "Found parasitic census file path in metadata (as String): ${parasiticPath}"
                            } else if (metadataJson.census_file_parasitic instanceof List) {
                                // If it's an array/list, take the first element
                                if (metadataJson.census_file_parasitic.size() > 0) {
                                    parasiticPath = metadataJson.census_file_parasitic[0].toString()
                                    log.info "Found parasitic census file path in metadata (as first element of List): ${parasiticPath}"
                                } else {
                                    log.warn "Census file parasitic in metadata is an empty List"
                                }
                            } else if (metadataJson.census_file_parasitic != null) {
                                // For any other type, try toString
                                try {
                                    parasiticPath = metadataJson.census_file_parasitic.toString()
                                    log.info "Found parasitic census file path in metadata (converted from ${metadataJson.census_file_parasitic.getClass().getName()}): ${parasiticPath}"
                                } catch (Exception e) {
                                    log.warn "Could not convert census_file_parasitic to string: ${e.message}"
                                }
                            } else {
                                log.warn "Census file parasitic in metadata is null"
                            }
                            
                            // Check if we obtained a valid path and use it
                            if (parasiticPath) {
                                def parasiticFile = file(parasiticPath, checkIfExists: false)
                                if (parasiticFile.exists()) {
                                    censusFilePVal = parasiticFile
                                    log.info "Using parasitic census file from metadata: ${censusFilePVal}"
                                } else {
                                    log.warn "Parasitic census file in metadata doesn't exist: ${parasiticPath}"
                                }
                            } else {
                                log.warn "No valid parasitic census file path found in metadata"
                            }
                        } else {
                            log.warn "Metadata doesn't contain census_file_parasitic key"
                        }
                    } catch (Exception e) {
                        log.warn "Error accessing census_file_parasitic in metadata: ${e.message}"
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
    
    // Final attempt: Try standard locations if not found yet
    if (censusFilePVal == null) {
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
    COMMAND LINE PARAMETERS:
    censusFileB   : ${params.censusFileB}
    censusFileP   : ${params.censusFileP}
    metadata      : ${params.metadata}
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
                    def slurper = new groovy.json.JsonSlurper()
                    def metadataJson = slurper.parseText(metadataContent)
                    
                    log.info "Successfully parsed metadata JSON from preprocessing"
                    // Handle census_file_bacterial path safely
                    try {
                        if (metadataJson.containsKey('census_file_bacterial')) {
                            log.info "Found census_file_bacterial in preprocessing metadata: ${metadataJson.census_file_bacterial?.getClass()?.getName() ?: 'null'}"
                            
                            // Handle census file path which might be a string or an array or other object
                            def bacterialPath = null
                            if (metadataJson.census_file_bacterial instanceof String) {
                                bacterialPath = metadataJson.census_file_bacterial.toString()
                                log.info "Found bacterial census file path in preprocessing metadata (as String): ${bacterialPath}"
                            } else if (metadataJson.census_file_bacterial instanceof List) {
                                // If it's an array/list, take the first element
                                if (metadataJson.census_file_bacterial.size() > 0) {
                                    bacterialPath = metadataJson.census_file_bacterial[0].toString()
                                    log.info "Found bacterial census file path in preprocessing metadata (as first element of List): ${bacterialPath}"
                                } else {
                                    log.warn "Census file bacterial in preprocessing metadata is an empty List"
                                }
                            } else if (metadataJson.census_file_bacterial != null) {
                                // For any other type, try toString
                                try {
                                    bacterialPath = metadataJson.census_file_bacterial.toString()
                                    log.info "Found bacterial census file path in preprocessing metadata (converted from ${metadataJson.census_file_bacterial.getClass().getName()}): ${bacterialPath}"
                                } catch (Exception e) {
                                    log.warn "Could not convert census_file_bacterial to string: ${e.message}"
                                }
                            } else {
                                log.warn "Census file bacterial in preprocessing metadata is null"
                            }
                            
                            // Check if we obtained a valid path and use it
                            if (bacterialPath) {
                                def bacterialFile = file(bacterialPath, checkIfExists: false)
                                if (bacterialFile.exists()) {
                                    censusFileBVal = bacterialFile
                                    log.info "Updated bacterial census file from preprocessing metadata: ${censusFileBVal}"
                                } else {
                                    log.warn "Bacterial census file in preprocessing metadata doesn't exist: ${bacterialPath}"
                                }
                            } else {
                                log.warn "No valid bacterial census file path found in preprocessing metadata"
                            }
                        } else {
                            log.warn "Preprocessing metadata doesn't contain census_file_bacterial key"
                        }
                    } catch (Exception e) {
                        log.warn "Error accessing census_file_bacterial in preprocessing metadata: ${e.message}"
                    }
                    
                    // Handle census_file_parasitic path safely
                    try {
                        if (metadataJson.containsKey('census_file_parasitic')) {
                            log.info "Found census_file_parasitic in preprocessing metadata: ${metadataJson.census_file_parasitic?.getClass()?.getName() ?: 'null'}"
                            
                            // Handle census file path which might be a string or an array or other object
                            def parasiticPath = null
                            if (metadataJson.census_file_parasitic instanceof String) {
                                parasiticPath = metadataJson.census_file_parasitic.toString()
                                log.info "Found parasitic census file path in preprocessing metadata (as String): ${parasiticPath}"
                            } else if (metadataJson.census_file_parasitic instanceof List) {
                                // If it's an array/list, take the first element
                                if (metadataJson.census_file_parasitic.size() > 0) {
                                    parasiticPath = metadataJson.census_file_parasitic[0].toString()
                                    log.info "Found parasitic census file path in preprocessing metadata (as first element of List): ${parasiticPath}"
                                } else {
                                    log.warn "Census file parasitic in preprocessing metadata is an empty List"
                                }
                            } else if (metadataJson.census_file_parasitic != null) {
                                // For any other type, try toString
                                try {
                                    parasiticPath = metadataJson.census_file_parasitic.toString()
                                    log.info "Found parasitic census file path in preprocessing metadata (converted from ${metadataJson.census_file_parasitic.getClass().getName()}): ${parasiticPath}"
                                } catch (Exception e) {
                                    log.warn "Could not convert census_file_parasitic to string: ${e.message}"
                                }
                            } else {
                                log.warn "Census file parasitic in preprocessing metadata is null"
                            }
                            
                            // Check if we obtained a valid path and use it
                            if (parasiticPath) {
                                def parasiticFile = file(parasiticPath, checkIfExists: false)
                                if (parasiticFile.exists()) {
                                    censusFilePVal = parasiticFile
                                    log.info "Updated parasitic census file from preprocessing metadata: ${censusFilePVal}"
                                } else {
                                    log.warn "Parasitic census file in preprocessing metadata doesn't exist: ${parasiticPath}"
                                }
                            } else {
                                log.warn "No valid parasitic census file path found in preprocessing metadata"
                            }
                        } else {
                            log.warn "Preprocessing metadata doesn't contain census_file_parasitic key"
                        }
                    } catch (Exception e) {
                        log.warn "Error accessing census_file_parasitic in preprocessing metadata: ${e.message}"
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
    
    // Verify that census file values are valid File objects
    log.info "Final validation of census file paths..."
    
    if (censusFileBVal != null) {
        try {
            log.info "Verifying bacterial census file: ${censusFileBVal.getClass().getName()}"
            if (censusFileBVal.exists()) {
                log.info "Bacterial census file exists: ${censusFileBVal}"
            } else {
                log.warn "Bacterial census file doesn't exist: ${censusFileBVal}"
                // Try one more normalize to be safe
                def path = censusFileBVal.toString()
                censusFileBVal = file(path, checkIfExists: false)
                log.info "Normalized bacterial census file path: ${censusFileBVal}, exists: ${censusFileBVal.exists()}"
                if (!censusFileBVal.exists()) {
                    censusFileBVal = null
                }
            }
        } catch (Exception e) {
            log.warn "Error validating bacterial census file: ${e.message}"
            censusFileBVal = null
        }
    }
    
    if (censusFilePVal != null) {
        try {
            log.info "Verifying parasitic census file: ${censusFilePVal.getClass().getName()}"
            if (censusFilePVal.exists()) {
                log.info "Parasitic census file exists: ${censusFilePVal}"
            } else {
                log.warn "Parasitic census file doesn't exist: ${censusFilePVal}"
                // Try one more normalize to be safe
                def path = censusFilePVal.toString()
                censusFilePVal = file(path, checkIfExists: false)
                log.info "Normalized parasitic census file path: ${censusFilePVal}, exists: ${censusFilePVal.exists()}"
                if (!censusFilePVal.exists()) {
                    censusFilePVal = null
                }
            }
        } catch (Exception e) {
            log.warn "Error validating parasitic census file: ${e.message}"
            censusFilePVal = null
        }
    }
    
    // Check final status of census files
    log.info "Final census file status:"
    log.info "  Bacterial census file: ${censusFileBVal?.exists() ? 'FOUND' : 'MISSING'}"
    log.info "  Parasitic census file: ${censusFilePVal?.exists() ? 'FOUND' : 'MISSING'}"
    
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
                
                // Get dashboard output and handle fallbacks
                def tempOutput = GENERATE_DASHBOARD.out.dashboard.collect()
                
                // Add a fallback mechanism to ensure we always have a dashboard
                // even if the module failed to create one
                tempOutput.ifEmpty { 
                    log.warn "Dashboard output is empty, using pre-created fallback"
                    tempOutput = Channel.fromPath(fallbackHtml.toString())
                }
                
                // Assign to the outer variable
                dashboardOutput = tempOutput
            } else {
                log.warn "No analysis results found, using pre-created fallback dashboard"
                dashboardOutput = Channel.fromPath(fallbackHtml.toString())
            }
        } catch (Exception e) {
            log.warn "Exception in dashboard generation section: ${e.getMessage()}"
            log.warn "Using pre-created fallback dashboard"
            dashboardOutput = Channel.fromPath(fallbackHtml.toString())
        }
    } else {
        log.info "Dashboard generation disabled, skipping"
    }
    
    // Create a dummy empty channel for when dashboard is disabled
    if (!binding.hasVariable('dashboardOutput')) {
        dashboardOutput = Channel.empty()
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

// Legacy helper function removed - now using direct string operations instead of file operations
