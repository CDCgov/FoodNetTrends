process GENERATE_DASHBOARD {
    tag "Generate dashboard"
    label 'process_medium'
    shell "/bin/bash"
    container 'foodnet.sif'

    input:
    val outDir
    val projID
    path dashboardTemplate

    output:
    path "${projID}_dashboard.html", emit: dashboard

    script:
    """
    Rscript \
      {workflow.projectDir}/bin/generate_dashboard.R \
      --outDir=${outDir} \
      --resultDir=${outDir} \
      --outputFile=${projID}_dashboard.html \
      --title="FoodNet Trends Analysis: ${projID}" \
      --templateFile=${dashboardTemplate}
    """
} 