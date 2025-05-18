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

// Include the workflows
include { SPLINE } from './workflows/spline.nf'
include { PREPROCESS_WORKFLOW } from './workflows/preprocess.nf'

// Default workflow
workflow {
    SPLINE()
}

// PREPROCESS_ONLY workflow entry point
workflow PREPROCESS_ONLY {
    PREPROCESS_WORKFLOW()
}
