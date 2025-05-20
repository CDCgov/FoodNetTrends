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
 * Last updated: 2025-05-20
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
    // Get file extension to handle CSV vs SAS files properly
    def mmwrExt = mmwrFile.toString().toLowerCase().endsWith('.csv') ? 'csv' : 'sas7bdat'
    """
    set -e
    echo "Starting analysis for pathogen: ${pathogen}" > ${pathogen}_trendy.log
    echo "Using MMWR data file: ${mmwrFile} (${mmwrExt} format)" >> ${pathogen}_trendy.log
    
    # Make a local copy of the MMWR file to handle path issues
    cp -v "${mmwrFile}" ./input_data.${mmwrExt}
    echo "Created local copy of MMWR file as: input_data.${mmwrExt}" >> ${pathogen}_trendy.log
    
    # Create empty placeholder files if needed
    if [ ! -f "${censusFileBact}" ] || [ ! -s "${censusFileBact}" ]; then
        echo "Census bacterial file missing or empty, creating placeholder" >> ${pathogen}_trendy.log
        echo "state,population,year,pathogentype" > empty_census_bact.csv
        echo "CA,10000000,2020,Bacterial" >> empty_census_bact.csv
        CENSUS_B_ARG="--censusFileB=empty_census_bact.csv"
    else
        # Create local copy of census file for consistent handling
        cp -v "${censusFileBact}" ./census_bact.sas7bdat
        echo "Created local copy of census bacterial file as: census_bact.sas7bdat" >> ${pathogen}_trendy.log
        CENSUS_B_ARG="--censusFileB=census_bact.sas7bdat"
    fi
    
    if [ ! -f "${censusFileParas}" ] || [ ! -s "${censusFileParas}" ]; then
        echo "Census parasitic file missing or empty, creating placeholder" >> ${pathogen}_trendy.log
        echo "state,population,year,pathogentype" > empty_census_para.csv
        echo "CA,10000000,2020,Parasitic" >> empty_census_para.csv
        CENSUS_P_ARG="--censusFileP=empty_census_para.csv"
    else
        # Create local copy of census file for consistent handling
        cp -v "${censusFileParas}" ./census_para.sas7bdat
        echo "Created local copy of census parasitic file as: census_para.sas7bdat" >> ${pathogen}_trendy.log
        CENSUS_P_ARG="--censusFileP=census_para.sas7bdat"
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
    
    # Check that our local data file copy exists and is readable
    if [ ! -f "./input_data.${mmwrExt}" ] || [ ! -r "./input_data.${mmwrExt}" ]; then
        echo "ERROR: Local MMWR data file copy not found or not readable" >> ${pathogen}_trendy.log
        echo "Original file: ${mmwrFile}" >> ${pathogen}_trendy.log
        echo "Local copy attempt: ./input_data.${mmwrExt}" >> ${pathogen}_trendy.log
        echo "Current directory contents:" >> ${pathogen}_trendy.log
        ls -la ./ >> ${pathogen}_trendy.log
        exit 1
    fi
    
    # Add preprocessed flag based on file extension
    if [ "${mmwrExt}" == "csv" ]; then
        PREPROC_ARG="--preprocessed=TRUE --cleanFile=./input_data.csv"
        echo "Using preprocessed mode for CSV file" >> ${pathogen}_trendy.log
    else
        PREPROC_ARG="--preprocessed=FALSE"
        echo "Using raw data mode for SAS file" >> ${pathogen}_trendy.log
    fi
    
    # Execute with carefully quoted arguments
    Rscript "\${SCRIPT_PATH}" \\
        --pathogen="${pathogen}" \\
        --mmwrFile="./input_data.${mmwrExt}" \\
        \${CENSUS_B_ARG} \\
        \${CENSUS_P_ARG} \\
        \${PREPROC_ARG} \\
        --projID="${projID}" \\
        --travel="${filter_travel}" \\
        --cidt="${filter_cidt}" \\
        --outDir="./"
    
    # Check return code from R script
    R_STATUS=\$?
    if [ \$R_STATUS -ne 0 ]; then
        echo "ERROR: R script failed with exit code \$R_STATUS" >> ${pathogen}_trendy.log
        exit \$R_STATUS
    fi
    
    echo "Analysis completed successfully" >> ${pathogen}_trendy.log
    """
}
