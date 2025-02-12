process TRENDY {
    tag "$pathogen"
    label 'process_large'
    shell '/bin/bash'
    container "foodnet.sif"

    input:
      val pathogen
      path mmwrFile
      path censusFile_B
      path censusFile_P
      val travel
      val cidt
      val projID
      val whichScript
      val preprocessed
      val cleanFile

    output:
      // path "${projID}/SplineResults/${pathogen}_brm.Rds", emit: rds
      // path "${projID}/SplineResults/${pathogen}_IRCatch.csv", emit: csv
      // path "${projID}/SplineResults/*.png", emit: png

    script:
    """
    mkdir -p ${projID}/SplineResults

    # Now run Rscript using the default invocation; the updated PATH and R_LIBS should take effect.
    Rscript ${whichScript} \\
      --mmwrFile "$mmwrFile" \\
      --censusFile_B "$censusFile_B" \\
      --censusFile_P "$censusFile_P" \\
      --travel "${travel}" \\
      --cidt "${cidt}" \\
      --projID "${projID}" \\
      --outDir "${projID}/SplineResults" \\
      --pathogen "${pathogen}" \\
      --preprocessed ${preprocessed} \\
      --cleanFile "$cleanFile" \\
      --debug FALSE

    # paste ${projID}/SplineResults/EstIRRCatch*.csv > EstIRRCatch_summary.csv
    """
}

