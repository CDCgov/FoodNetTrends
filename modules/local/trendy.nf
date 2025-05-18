/*
 * ==================================================================
 * FoodNet Trends - TRENDY Process Module
 * ==================================================================
 *
 * Purpose:
 *   This process runs the main Bayesian modeling for each pathogen.
 *   It handles input data preparation, executes the R modeling script,
 *   captures detailed logs, and manages potential errors.
 *
 * Inputs:
 *   - Pathogen name (for targeting specific organism)
 *   - MMWR data file (preprocessed or raw)
 *   - Census files (bacterial and parasitic)
 *   - Travel and CIDT filtering parameters
 *   - Project ID and script paths
 *
 * Outputs:
 *   - Bayesian model RDS files
 *   - Incidence rate estimates CSV
 *   - Visualizations (PNG files)
 *   - Summary statistics and logs
 *   - Interactive HTML dashboard
 *
 * Error handling:
 *   - Retries on memory/resource errors
 *   - Detailed logging for diagnostics
 *
 * Last updated: 2025-05-18
 * ==================================================================
 */

process TRENDY {
    tag "$pathogen"
    label 'process_high_memory'
    label 'error_retry'
    
    container 'foodnet.sif'
    
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "*.log"
    publishDir "${params.outdir}/${params.projID}", mode: params.publish_dir_mode, pattern: "{*.png,*_summary.txt,*_IRCatch.csv,*.Rds,*_EstIRRCatch_*.csv}"
    publishDir "${params.outdir}/${params.projID}", mode: params.publish_dir_mode, pattern: "dashboard.html", saveAs: { "${params.projID}_" + new java.text.SimpleDateFormat("yyyyMMdd").format(new Date()) + "_dashboard.html" }
    
    input:
    tuple val(pathogen), path(mmwrFile)
    path censusFileBact
    path censusFileParas
    val projID
    path scripts_path
    val filter_travel
    val filter_cidt
    path template_path
    
    output:
    path "${pathogen}_${params.output_suffix}.Rds", emit: model
    path "${pathogen}_summary.txt", emit: summary
    path "${pathogen}_IRCatch.csv", emit: ir_outputs
    path "${pathogen}_*.png", emit: plots
    path "${pathogen}_EstIRRCatch_*.csv", emit: irr_outputs
    path "${params.projID}_dashboard.html", optional: true, emit: dashboard
    path "*.log", emit: logs
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    echo "Starting analysis for pathogen: ${pathogen}" > ${pathogen}_trendy.log
    echo "Using MMWR data file: ${mmwrFile}" >> ${pathogen}_trendy.log
    echo "Using Census bacterial file: ${censusFileBact}" >> ${pathogen}_trendy.log
    echo "Using Census parasitic file: ${censusFileParas}" >> ${pathogen}_trendy.log
    
    # Run the main trend analysis
    Rscript ${scripts_path}/trendy.R \
        --pathogen=${pathogen} \
        --mmwrFile=${mmwrFile} \
        --censusBact=${censusFileBact} \
        --censusParas=${censusFileParas} \
        --projID=${projID} \
        --filter_travel=${filter_travel} \
        --filter_cidt=${filter_cidt} \
        --outDir=./
    
    # If this is the last pathogen (determined by task ID), generate the dashboard
    if [[ \$TASK_ID == \$TASK_MAX ]]; then
        echo "Generating interactive dashboard..." >> ${pathogen}_trendy.log
        Rscript ${scripts_path}/generate_dashboard.R \
            --outDir=./ \
            --resultDir=./ \
            --outputFile=${params.projID}_dashboard.html \
            --title="FoodNet Trends Analysis: ${params.projID}" \
            --templateFile=${template_path}/dashboard_template.html || \
            echo "Warning: Dashboard generation failed, proceeding with standard outputs" >> ${pathogen}_trendy.log
    fi
    
    echo "Analysis completed successfully" >> ${pathogen}_trendy.log
    """
}
