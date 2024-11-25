process TRENDY {
    tag "trendy"
    label 'process_large'

    conda "${baseDir}/assets/Renv.yaml"

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
    Rscript "trendy.R" \
        --mmwrFile "$mmwrFile" \
        --censusFile_B "$censusFile_B" \
        --censusFile_P "$censusFile_P" \
        --travel "${travel}" \
        --cidt "${cidt}" \
        --projID "${projID}" \
        --whichFunctions "functions.R" \
        --debug FALSE

    # Combine CSV files into a summary
    paste SplineResults/EstIRRCatch*.csv > EstIRRCatch_summary.csv
    """
}

