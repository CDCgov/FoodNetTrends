process GENERATE_DASHBOARD {
    tag "Generate dashboard"
    label 'process_medium'
    shell "/bin/bash"
    container 'foodnet.sif'

    input:
    path ir_outputs
    val resultDir
    val projID
    path dashboardTemplate
    val dashboardScript

    output:
    path "${projID}_dashboard.html", emit: dashboard
    publishDir "${params.outdir}/${params.projID}", mode: params.publish_dir_mode

    script:
    """
    Rscript \
      ${dashboardScript} \
      --outDir=${resultDir} \
      --resultDir=${resultDir} \
      --outputFile=${projID}_dashboard.html \
      --title="FoodNet Trends Analysis: ${projID}" \
      --templateFile=${dashboardTemplate}
    """
} 