nextflow.enable.dsl = 2

// Handle help parameter
if (params.help) {
    println """
    ============================================
    FoodNet Trends Pipeline
    ============================================
    Usage:
        nextflow run main.nf [options]

    Options:
        --help                  Show this help message
        --mmwrFile              Path to the MMWR data file (raw SAS or preprocessed CSV)
        --censusFileB          Path to the bacterial census data file
        --censusFileP          Path to the parasitic census data file
        --travel                Travel types (e.g., NO,UNKNOWN)
        --cidt                  CIDT types (e.g., CIDT+,CX+,PARASITIC)
        --projID                Project ID (e.g., 20240705)
        --outdir                Base output directory for pipeline reports and results
        --outputBase            Base name for output files (preprocessing mode)
        --preprocessed          TRUE/FALSE indicating if using preprocessed CSV data
        --metadata              Path to metadata JSON file (for preprocessed data)
        --pathogen              Comma-separated list of pathogens to analyze (e.g., CAMPYLOBACTER,SALMONELLA)
        --states                Comma-separated list of states to analyze (e.g., CA,NY,GA)
        --salmonella_serotypes  Comma-separated list of Salmonella serotypes to analyze
        --chains                Number of MCMC chains
        --iterations            Number of MCMC iterations
        --adapt_delta           Adaptation parameter for MCMC
        --max_treedepth         Maximum tree depth for MCMC
        --seed                  Random seed for reproducibility
    """
    System.exit(0)
}

// Handle common file errors
def validateFilePaths() {
    // Validate mmwrFile first - always required
    if (!params.mmwrFile) {
        error "Missing required parameter: --mmwrFile"
    }
    
    // Check if the path is a valid string that can be used safely
    try {
        def mmwrFilePath = params.mmwrFile.toString()
        log.info "MMWR file path appears valid: ${mmwrFilePath}"
    } catch (Exception e) {
        error "Invalid MMWR file path format: ${params.mmwrFile}"
    }
    
    // Census files are now optional with warnings
    if (params.censusFileB && params.censusFileB != "") {
        try {
            def censusPath = params.censusFileB.toString()
            log.info "Census bacterial file path appears valid: ${censusPath}"
        } catch (Exception e) {
            log.warn "Invalid census bacterial file path format: ${params.censusFileB}"
        }
    }
    
    if (params.censusFileP && params.censusFileP != "") {
        try {
            def censusPath = params.censusFileP.toString()
            log.info "Census parasitic file path appears valid: ${censusPath}"
        } catch (Exception e) {
            log.warn "Invalid census parasitic file path format: ${params.censusFileP}"
        }
    }
}

// Include the workflows
include { SPLINE } from './workflows/spline.nf'
include { PREPROCESS_WORKFLOW } from './workflows/preprocess.nf'

// Default workflow
workflow {
    try {
        // Validate file paths first
        validateFilePaths()
        
        // Then run main workflow
        SPLINE()
    } catch (Exception e) {
        log.error "Error in SPLINE workflow: ${e.message}"
        if (e.message =~ /(?i).*filesystem.*/ || e.message =~ /(?i).*getFileSystem.*/) {
            log.error """
            ==========================================
            File System Error Detected
            ==========================================
            This appears to be a file path handling issue. Check that:
            
            1. Your input paths exist and are accessible
            2. You're using proper path syntax for your system
            3. File names don't contain special characters
            4. For preprocessed files, ensure they were generated with compatible versions
            
            Try running with the -entry PREPROCESS_ONLY flag first to generate 
            fresh preprocessed files, then run the main workflow with --preprocessed true.
            ==========================================
            """
        }
        System.exit(1)
    }
}

// PREPROCESS_ONLY workflow entry point
workflow PREPROCESS_ONLY {
    try {
        // Validate file paths first
        validateFilePaths()
        
        // Then run preprocessing workflow
        PREPROCESS_WORKFLOW()
    } catch (Exception e) {
        log.error "Error in PREPROCESS_WORKFLOW: ${e.message}"
        if (e.message =~ /(?i).*filesystem.*/ || e.message =~ /(?i).*getFileSystem.*/) {
            log.error """
            ==========================================
            File System Error Detected
            ==========================================
            This appears to be a file path handling issue. Check that:
            
            1. Your input paths exist and are accessible
            2. You're using proper path syntax for your system
            3. File names don't contain special characters
            
            Additional debugging tips:
            - Try using absolute paths instead of relative ones
            - Verify paths don't have spaces or special characters
            - Check if your HPC environment has any path restrictions
            ==========================================
            """
        }
        System.exit(1)
    }
}
