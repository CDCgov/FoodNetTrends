nextflow.enable.dsl = 2

// Handle help parameter
if (params.help) {
    println """
    ============================================
    cdc/spline Pipeline
    ============================================
    Usage:
        nextflow run main.nf [options]

    Options:
        --help                  Show this help message
        --mmwrFile              Path to the MMWR data file
        --censusFile_B          Path to the bacterial census data file
        --censusFile_P          Path to the parasitic census data file
        --travel                Travel types (e.g., NO,UNKNOWN)
        --cidt                  CIDT types (e.g., CIDT+,CX+,PARASITIC)
        --projID                Project ID
    """
    System.exit(0)
}

// Include the SPLINE workflow
include { SPLINE } from './workflows/spline.nf'

workflow CDC_SPLINE {
    SPLINE()
}

