process PREPROCESS {
    tag "Preprocessing MMWR data"
    label 'process_medium'
    shell "/bin/bash"
    container 'foodnet.sif'

    publishDir "${params.outdir}/${projID}/preprocessed", mode: 'copy'

    input:
    path mmwrFile
    val projID

    output:
    path "clean_mmwr.csv", emit: cleanFile

    script:
    // Use absolute path to the script or a relative path from the current directory
    def scriptPath = "${workflow.projectDir}/bin/calcIR.R"

    """
    Rscript ${scriptPath} \\
      --mmwrFile ${mmwrFile} \\
      --outputFile clean_mmwr.csv
    """
}
