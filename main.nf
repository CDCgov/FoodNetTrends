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
        --mmwrFile              Path to the MMWR data file (or cleaned CSV if preprocessed is TRUE)
        --censusFile_B          Path to the bacterial census data file
        --censusFile_P          Path to the parasitic census data file
        --travel                Travel types (e.g., NO,UNKNOWN)
        --cidt                  CIDT types (e.g., CIDT+,CX+,PARASITIC)
        --projID                Project ID (e.g., 20240705)
        --outdir                Base output directory for pipeline reports and results
        --preprocessed          TRUE/FALSE indicating if using preprocessed CSV data
        --cleanFile             Path to cleaned CSV file (if preprocessed is TRUE)
        --pathogen              Comma-separated list of pathogens to analyze (e.g., CAMPYLOBACTER,SALMONELLA)
        --pathogen_grouping     Pipe-separated pathogen groupings (e.g., STEC:O157|STEC:nonO157)
        --chains                Number of MCMC chains
        --iterations            Number of MCMC iterations
        --adapt_delta           Adaptation parameter for MCMC
        --max_treedepth         Maximum tree depth for MCMC
        --seed                  Random seed for reproducibility
    """
    System.exit(0)
}

// Include the SPLINE workflow
include { SPLINE } from './workflows/spline.nf'

workflow FoodNetTrends {
    SPLINE()
}
