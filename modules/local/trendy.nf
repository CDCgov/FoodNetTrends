process TRENDY {
    tag "$pathogen"
    label 'process_large'
    shell '/bin/bash'
    container "foodnet.sif"
    publishDir "${params.outdir}/${projID}/spline_results", mode: 'copy'

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
      path "${pathogen}_brm.Rds", emit: rds
      path "${pathogen}_IRCatch.csv", emit: csv
      path "${pathogen}_*.png", emit: png

    script:
    """
    # Run Rscript using the default invocation
    Rscript ${whichScript} \\
      --mmwrFile "$mmwrFile" \\
      --censusFile_B "$censusFile_B" \\
      --censusFile_P "$censusFile_P" \\
      --travel "${travel}" \\
      --cidt "${cidt}" \\
      --projID "${projID}" \\
      --outDir "./" \\
      --pathogen "${pathogen}" \\
      --preprocessed ${preprocessed} \\
      --cleanFile "$cleanFile" \\
      --debug FALSE
    """
}

