/*
========================================
  FoodNet Trends v1.0
========================================

  A Nextflow pipeline for Bayesian modeling of foodborne disease surveillance data
  for the CDC Foodborne Diseases Active Surveillance Network (FoodNet).

----------------------------------------
  Main workflow
----------------------------------------
*/

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
        --cores                 Number of CPU cores to use per pathogen analysis
        --enable_dashboard      TRUE/FALSE to enable/disable dashboard generation
        --dashboard_title       Title for dashboard
    ============================================
    """
    exit 0
}

// Include the workflows
include { SPLINE } from './workflows/spline_fixed.nf'
include { PREPROCESS_WORKFLOW } from './workflows/preprocess.nf'

// Determine workflow entry point based on params
workflow {
    if (params.mmwrFile) {
        log.info "MMWR file provided: ${params.mmwrFile}"
        
        // Add validation of key parameters
        validateParams()
        
        // Determine entry point
        if (params.skip_preprocessing) {
            log.info "Skipping preprocessing as requested."
            SPLINE()
        } else if (params.preprocessing_only) {
            log.info "Running preprocessing workflow only."
            PREPROCESS_WORKFLOW()
        } else {
            log.info "Running full workflow (preprocessing + analysis)."
            SPLINE()
        }
    } else {
        error "No MMWR file provided! Use --mmwrFile to specify an input dataset."
    }
}

def validateParams() {
    // Essential parameter validation for the FoodNet Trends pipeline
    if (!params.mmwrFile) {
        error "Missing required parameter: --mmwrFile"
    }
    
    if (!params.pathogen) {
        error "Missing required parameter: --pathogen (comma-separated list of pathogens to analyze)"
    }
    
    // Check for valid preprocessing flag settings
    if (params.skip_preprocessing && params.preprocessing_only) {
        error "Incompatible options: Cannot specify both --skip_preprocessing and --preprocessing_only"
    }
    
    if (params.preprocessing_only && !params.outputBase) {
        log.warn "No output base name specified for preprocessing. Using default derived from input filename."
    }
    
    // Validate workflow-specific parameters
    if (!params.skip_preprocessing && !params.preprocessing_only) {
        log.info "Running complete workflow (preprocessing + analysis)"
        // No additional validation needed
    } else if (params.skip_preprocessing) {
        log.info "Skipping preprocessing, running analysis directly"
        // Ensure we have a preprocessed data source
        if (!params.cleanFile && !params.preprocessed) {
            error "When skipping preprocessing, you must specify either --cleanFile or set --preprocessed=true"
        }
    } else if (params.preprocessing_only) {
        log.info "Running preprocessing only"
        // No additional validation needed
    }
    
    log.info "Parameter validation complete"
}