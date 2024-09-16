process TRENDY {
    tag "trendy"
    label 'process_large'

    conda "assets/Renv.yaml"

    input:
    path whichScript
    path whichFunctions
    path mmwrFile
    path censusFile_B
    path censusFile_P
    val(travel)
    val(cidt)
    val(projID)

    output:
    path 'SplineResults/*.Rds'       , emit: Rds
    path 'EstIRRCatch_summary.csv'       , emit: csv
    path 'SplineResults/*.png'       , emit: png
    
    script:
    """
    Rscript "${whichScript}" \
        --debug FALSE \
        --mmwrFile $mmwrFile \
        --censusFile_B $censusFile_B \
        --censusFile_P $censusFile_P \
        --travel "${travel}" \
        --cidt "${cidt}" \
        --projID $projID

    paste SplineResults/EstIRRCatch*.csv > EstIRRCatch_summary.csv
    """
}
