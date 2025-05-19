process GENERATE_DASHBOARD {
    tag "Generate dashboard"
    label 'process_medium'
    shell "/bin/bash"
    container 'foodnet.sif'

    input:
    path ir_outputs
    val outDir
    val projID
    path dashboardTemplate
    val dashboardScript

    output:
    path "${projID}_dashboard.html", emit: dashboard

    script:
    """
    Rscript \
      ${dashboardScript} \
      --outDir=${outDir} \
      --resultDir=${outDir} \
      --outputFile=${projID}_dashboard.html \
      --title="FoodNet Trends Analysis: ${projID}" \
      --templateFile=${dashboardTemplate}
    """
} 