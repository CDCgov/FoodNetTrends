nextflow.enable.dsl=2

process TRENDY {
    tag "trendy"
    label 'process_large'

    // conda "${baseDir}/assets/Renv.yaml"

    input:
    path whichScript
    path whichFunctions
    path mmwrFile
    path censusFile_B
    path censusFile_P
    val travel
    val cidt
    val projID

    output:
    path 'SplineResults/*.Rds', emit: Rds
    path 'EstIRRCatch_summary.csv', emit: csv
    path 'SplineResults/*.png', emit: png

    script:
    """
    # Run the R script with the provided arguments
    Rscript "${whichScript}" \
        --mmwrFile "${mmwrFile}" \
        --censusFile_B "${censusFile_B}" \
        --censusFile_P "${censusFile_P}" \
        --travel "${travel}" \
        --cidt "${cidt}" \
        --projID "${projID}" \
        --whichFunctions "${whichFunctions}"

    # Combine CSV files into a summary
    paste SplineResults/EstIRRCatch*.csv > EstIRRCatch_summary.csv
    """
}

workflow {
    TRENDY(
        file(params.whichScript),
        file(params.whichFunctions),
        file(params.mmwrFile),
        file(params.censusFile_B),
        file(params.censusFile_P),
        params.travel,
        params.cidt,
        params.projID
    )
}

