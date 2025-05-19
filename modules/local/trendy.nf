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
    publishDir "${params.outdir}/${params.projID}", mode: params.publish_dir_mode, pattern: "dashboard.html", saveAs: { "${params.projID}_dashboard.html" }
    
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
    set -e
    echo "Starting analysis for pathogen: ${pathogen}" > ${pathogen}_trendy.log
    echo "Using MMWR data file: ${mmwrFile}" >> ${pathogen}_trendy.log
    
    # Create empty placeholder files if needed
    if [ ! -f "${censusFileBact}" ] || [ ! -s "${censusFileBact}" ]; then
        echo "Census bacterial file missing or empty, creating placeholder" >> ${pathogen}_trendy.log
        echo "state,population,year,pathogentype" > empty_census_bact.csv
        echo "CA,10000000,2020,Bacterial" >> empty_census_bact.csv
        CENSUS_B_ARG="--censusFileB=empty_census_bact.csv"
    else
        echo "Census bacterial file exists: ${censusFileBact}" >> ${pathogen}_trendy.log
        CENSUS_B_ARG="--censusFileB=\"${censusFileBact}\""
    fi
    
    if [ ! -f "${censusFileParas}" ] || [ ! -s "${censusFileParas}" ]; then
        echo "Census parasitic file missing or empty, creating placeholder" >> ${pathogen}_trendy.log
        echo "state,population,year,pathogentype" > empty_census_para.csv
        echo "CA,10000000,2020,Parasitic" >> empty_census_para.csv
        CENSUS_P_ARG="--censusFileP=empty_census_para.csv"
    else
        echo "Census parasitic file exists: ${censusFileParas}" >> ${pathogen}_trendy.log
        CENSUS_P_ARG="--censusFileP=\"${censusFileParas}\""
    fi
    
    # Run the main trend analysis with explicit path handling for everything
    echo "Using scripts path: ${scripts_path}" >> ${pathogen}_trendy.log
    
    # Create explicit path to R script
    SCRIPT_PATH="${scripts_path}/trendy.R"
    echo "Full script path: \${SCRIPT_PATH}" >> ${pathogen}_trendy.log
    
    # Check that R script exists
    if [ ! -f "\${SCRIPT_PATH}" ]; then
        echo "ERROR: R script not found at \${SCRIPT_PATH}" >> ${pathogen}_trendy.log
        echo "Directory contents of ${scripts_path}:" >> ${pathogen}_trendy.log
        ls -la "${scripts_path}" >> ${pathogen}_trendy.log
        exit 1
    fi
    
    # Execute with carefully quoted arguments
    Rscript "\${SCRIPT_PATH}" \\
        --pathogen="${pathogen}" \\
        --mmwrFile="${mmwrFile}" \\
        \${CENSUS_B_ARG} \\
        \${CENSUS_P_ARG} \\
        --projID="${projID}" \\
        --travel="${filter_travel}" \\
        --cidt="${filter_cidt}" \\
        --outDir="./"
    
    echo "Analysis completed successfully" >> ${pathogen}_trendy.log
    """
}
