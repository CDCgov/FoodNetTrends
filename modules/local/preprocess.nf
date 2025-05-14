process PREPROCESS {
    tag "Preprocessing MMWR data"
    label 'process_medium'
    shell "/bin/bash"
    container 'foodnet.sif'

    publishDir "${params.outdir}/preprocessed", mode: 'copy'

    input:
    path mmwrFile
    val outputBase
    val generateMetadata

    output:
    path "${outputBase}.csv", emit: cleanedData
    path "${outputBase}_metadata.json", optional: true, emit: metadata

    script:
    """
    # Use calcIR.R with metadata generation enabled if requested
    Rscript ${workflow.projectDir}/bin/calcIR.R \
      --mmwrFile ${mmwrFile} \
      --outputFile ${outputBase}.csv \
      --generate_metadata ${generateMetadata}
    """
}
