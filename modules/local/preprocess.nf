process PREPROCESS {
    tag "Preprocessing MMWR data"
    label 'process_medium'
    shell '/bin/bash'
    container "foodnet.sif"
    publishDir "${params.outdir}/${params.projID}/preprocess_results", mode: 'copy'

    input:
      path mmwrFile
      val projID

    output:
      path "${projID}_clean_mmwr.csv", emit: clean_csv

    script:
    """
    Rscript ${launchDir}/bin/calcIR.R \\
      --mmwrFile "$mmwrFile" \\
      --outputFile "${projID}_clean_mmwr.csv"
    """
}

