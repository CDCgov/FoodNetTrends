/*
 * ==================================================================
 * FoodNetTrends v1.0 - Bayesian Modeling Module
 * ==================================================================
 *
 * Purpose:
 *   Executes Bayesian hierarchical spline models for individual pathogens.
 *   Handles data preparation, model fitting, and result generation with
 *   specialized processing for different pathogen types.
 *
 * Inputs:
 *   - Target pathogen identifier
 *   - Preprocessed MMWR surveillance data
 *   - Population census data (bacterial/parasitic)
 *   - Analysis parameters and filtering criteria
 *
 * Outputs:
 *   - Fitted Bayesian model objects (RDS)
 *   - Incidence rate estimates with confidence intervals
 *   - Trend visualizations and summary statistics
 *   - Detailed execution logs
 *
 * Last updated: 2025-05-22
 * ==================================================================
 */

process TRENDY {
    tag "${pathogen}"
    label 'process_high_memory'
    label 'error_retry'
    
    container 'foodnet.sif'
    
    publishDir "${params.outdir}/results/${pathogen}", 
        pattern: "${pathogen}*.{csv,txt,log}", 
        mode: 'copy'
    publishDir "${params.outdir}/figures/${pathogen}", 
        pattern: "${pathogen}*.{png,pdf}", 
        mode: 'copy'
    publishDir "${params.outdir}/models/${pathogen}", 
        pattern: "${pathogen}*.{Rds,rds}", 
        mode: 'copy'
    
    input:
    val pathogen
    path mmwrFile
    path censusBFile
    path censusPFile
    path bin_dir
    
    output:
    path "${pathogen}_summary.txt", emit: summary
    path "${pathogen}_brm.Rds", emit: model, optional: true
    path "${pathogen}_*.csv", emit: results, optional: true
    path "${pathogen}_*.png", emit: figures, optional: true
    path "${pathogen}_*.log", emit: log
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def mmwrExt = mmwrFile.toString().endsWith('.csv') ? 'csv' : 'sas7bdat'
    def mmwrIsPreprocessed = mmwrExt == 'csv'
    def localMmwrFile = "input_data.${mmwrExt}"
    """
    #!/usr/bin/env bash
    set -e  # Exit immediately if a command exits with non-zero status
    
    # Create log file
    echo "Starting analysis for pathogen: ${pathogen}" > ${pathogen}_trendy.log
    echo "CPU cores: ${task.cpus}, Memory: ${task.memory}" >> ${pathogen}_trendy.log
    
    # Create error function
    error_exit() {
        # Log error message
        echo "ERROR: \$1" | tee -a "${pathogen}_trendy.log"
        echo "Error in TRENDY process for ${pathogen}: \$1" >> "${pathogen}_error_summary.txt"
        
        # Create a basic summary file
        echo "Error processing ${pathogen}" > "${pathogen}_summary.txt"
        echo "Error occurred at `date`" >> "${pathogen}_summary.txt"
        echo "Error message: \$1" >> "${pathogen}_summary.txt"
        exit 1
    }
    
    # Set error trap
    trap 'error_exit "Command failed with exit code \$?: \$BASH_COMMAND"' ERR
    
    # Verify MMWR file exists before proceeding
    if [ ! -f "${mmwrFile}" ]; then
        error_exit "MMWR data file does not exist: ${mmwrFile}"
    fi
    
    echo "Using MMWR data file: ${mmwrFile} (${mmwrExt} format)" | tee -a ${pathogen}_trendy.log
    
    # Make a local copy of the MMWR file to handle path issues
    if ! cp -v "${mmwrFile}" "${localMmwrFile}"; then
        error_exit "Failed to copy MMWR file - check file permissions and path"
    fi
    
    echo "Created local copy of MMWR file as: ${localMmwrFile}" | tee -a ${pathogen}_trendy.log
    
    # Validate file is not empty
    if [ ! -s "${localMmwrFile}" ]; then
        error_exit "MMWR data file is empty"
    fi
    
    # Initialize variable defaults (empty but will be set if required files exist)
    CENSUS_B_ARG=""
    CENSUS_P_ARG=""
    PREPROC_ARG=""
    
    # Handle census bacterial file
    if [ -f "${censusBFile}" ] && [ -s "${censusBFile}" ]; then
        # File exists and is not empty
        echo "Census bacterial file exists: ${censusBFile}" >> ${pathogen}_trendy.log
        
        # Get file extension from basename
        CENSUS_B_BASENAME=`basename "${censusBFile}"`
        CENSUS_B_EXT="\${CENSUS_B_BASENAME##*.}"
        echo "Census bacterial file extension: \${CENSUS_B_EXT}" >> ${pathogen}_trendy.log
        
        # Create appropriate local copy based on extension
        if [ "\${CENSUS_B_EXT}" = "csv" ]; then
            cp -v "${censusBFile}" ./census_bact.csv || error_exit "Failed to copy census bacterial CSV file"
            echo "Created local copy of census bacterial file as: census_bact.csv" >> ${pathogen}_trendy.log
            chmod 644 ./census_bact.csv
            CENSUS_B_ARG="--censusFileB=census_bact.csv"
        elif [ "\${CENSUS_B_EXT}" = "sas7bdat" ]; then
            cp -v "${censusBFile}" ./census_bact.sas7bdat || error_exit "Failed to copy census bacterial SAS file"
            echo "Created local copy of census bacterial file as: census_bact.sas7bdat" >> ${pathogen}_trendy.log
            chmod 644 ./census_bact.sas7bdat
            CENSUS_B_ARG="--censusFileB=census_bact.sas7bdat"
        else
            echo "Unknown census bacterial file extension, defaulting to CSV" >> ${pathogen}_trendy.log
            cp -v "${censusBFile}" ./census_bact.csv || error_exit "Failed to copy census bacterial file"
            echo "Created local copy of census bacterial file as: census_bact.csv" >> ${pathogen}_trendy.log
            chmod 644 ./census_bact.csv
            CENSUS_B_ARG="--censusFileB=census_bact.csv"
        fi
    else
        # Census file is required - exit with error
        error_exit "ERROR: Census bacterial file is required but was missing or empty"
    fi
    
    # Handle census parasitic file
    if [ -f "${censusPFile}" ] && [ -s "${censusPFile}" ]; then
        # File exists and is not empty
        echo "Census parasitic file exists: ${censusPFile}" >> ${pathogen}_trendy.log
        
        # Get file extension from basename
        CENSUS_P_BASENAME=`basename "${censusPFile}"`
        CENSUS_P_EXT="\${CENSUS_P_BASENAME##*.}"
        echo "Census parasitic file extension: \${CENSUS_P_EXT}" >> ${pathogen}_trendy.log
        
        # Create appropriate local copy based on extension
        if [ "\${CENSUS_P_EXT}" = "csv" ]; then
            cp -v "${censusPFile}" ./census_para.csv || error_exit "Failed to copy census parasitic CSV file"
            echo "Created local copy of census parasitic file as: census_para.csv" >> ${pathogen}_trendy.log
            chmod 644 ./census_para.csv
            CENSUS_P_ARG="--censusFileP=census_para.csv"
        elif [ "\${CENSUS_P_EXT}" = "sas7bdat" ]; then
            cp -v "${censusPFile}" ./census_para.sas7bdat || error_exit "Failed to copy census parasitic SAS file"
            echo "Created local copy of census parasitic file as: census_para.sas7bdat" >> ${pathogen}_trendy.log
            chmod 644 ./census_para.sas7bdat
            CENSUS_P_ARG="--censusFileP=census_para.sas7bdat"
        else
            echo "Unknown census parasitic file extension, defaulting to CSV" >> ${pathogen}_trendy.log
            cp -v "${censusPFile}" ./census_para.csv || error_exit "Failed to copy census parasitic file"
            echo "Created local copy of census parasitic file as: census_para.csv" >> ${pathogen}_trendy.log
            chmod 644 ./census_para.csv
            CENSUS_P_ARG="--censusFileP=census_para.csv"
        fi
    else
        # Census parasitic file is required - exit with error
        error_exit "ERROR: Census parasitic file is required but was missing or empty"
    fi
    
    # Check data files and log their status
    echo "===== DATA FILES VERIFICATION =====" >> ${pathogen}_trendy.log
    if [ -f "./census_bact.sas7bdat" ] && [ -s "./census_bact.sas7bdat" ]; then
        echo "VALID: Using real bacterial census data (SAS format)" >> ${pathogen}_trendy.log
    elif [ -f "./census_bact.csv" ] && [ -s "./census_bact.csv" ]; then
        echo "VALID: Using real bacterial census data (CSV format)" >> ${pathogen}_trendy.log
    else 
        echo "ERROR: No valid bacterial census data found" >> ${pathogen}_trendy.log
        error_exit "ERROR: Bacterial census file required but not found or empty"
    fi
    
    if [ -f "./census_para.sas7bdat" ] && [ -s "./census_para.sas7bdat" ]; then
        echo "VALID: Using real parasitic census data (SAS format)" >> ${pathogen}_trendy.log
    elif [ -f "./census_para.csv" ] && [ -s "./census_para.csv" ]; then
        echo "VALID: Using real parasitic census data (CSV format)" >> ${pathogen}_trendy.log
    else
        echo "ERROR: No valid parasitic census data found" >> ${pathogen}_trendy.log
        error_exit "ERROR: Parasitic census file required but not found or empty"
    fi
    echo "===================================" >> ${pathogen}_trendy.log
    
    # Check local data file copy exists and is readable
    if [ ! -f "${localMmwrFile}" ] || [ ! -r "${localMmwrFile}" ]; then
        # List directory contents for debugging
        ls -la ./ > dir_contents.txt
        error_exit "Local MMWR data file copy not found or not readable. Original file: ${mmwrFile}, Local copy attempt: ${localMmwrFile}"
    fi
    
    # Add preprocessed flag based on file type
    if [ "${mmwrIsPreprocessed}" = "true" ]; then
        PREPROC_ARG="--preprocessed=TRUE --cleanFile=./${localMmwrFile}"
        echo "Using preprocessed mode for CSV file" >> ${pathogen}_trendy.log
    else
        PREPROC_ARG="--preprocessed=FALSE --rawFile=./${localMmwrFile}"
        echo "Using raw data mode for SAS file" >> ${pathogen}_trendy.log
    fi
    
    # Log census arguments
    echo "Census bacterial arg: \${CENSUS_B_ARG}" >> ${pathogen}_trendy.log
    echo "Census parasitic arg: \${CENSUS_P_ARG}" >> ${pathogen}_trendy.log
    echo "Preprocessing arg: \${PREPROC_ARG}" >> ${pathogen}_trendy.log
    
    # Create a results summary file
    echo "====================================================" > ${pathogen}_data_summary.txt
    echo "     ANALYSIS SETUP FOR PATHOGEN: ${pathogen}" >> ${pathogen}_data_summary.txt
    echo "====================================================" >> ${pathogen}_data_summary.txt
    echo "Data Sources:" >> ${pathogen}_data_summary.txt
    echo "  MMWR Data:       ${localMmwrFile}" >> ${pathogen}_data_summary.txt
    echo "  Census Data:" >> ${pathogen}_data_summary.txt
    
    # Check census data files and record their status
    echo "  Census Bacterial Status:" >> ${pathogen}_data_summary.txt
    if [ -f "./census_bact.sas7bdat" ] && [ -s "./census_bact.sas7bdat" ]; then
        echo "    REAL DATA (SAS format)" >> ${pathogen}_data_summary.txt
        ls -la "./census_bact.sas7bdat" >> ${pathogen}_data_summary.txt
    elif [ -f "./census_bact.csv" ] && [ -s "./census_bact.csv" ]; then
        echo "    REAL DATA (CSV format)" >> ${pathogen}_data_summary.txt
        ls -la "./census_bact.csv" >> ${pathogen}_data_summary.txt
    elif [ -f "./empty_census_bact.csv" ] && [ -s "./empty_census_bact.csv" ]; then
        echo "    *** ERROR: Invalid census file *** (Real census data is required)" >> ${pathogen}_data_summary.txt
        ls -la "./empty_census_bact.csv" >> ${pathogen}_data_summary.txt
    else
        echo "    *** NO CENSUS FILE FOUND ***" >> ${pathogen}_data_summary.txt
    fi
    
    echo "  Census Parasitic Status:" >> ${pathogen}_data_summary.txt
    if [ -f "./census_para.sas7bdat" ] && [ -s "./census_para.sas7bdat" ]; then
        echo "    REAL DATA (SAS format)" >> ${pathogen}_data_summary.txt
        ls -la "./census_para.sas7bdat" >> ${pathogen}_data_summary.txt
    elif [ -f "./census_para.csv" ] && [ -s "./census_para.csv" ]; then
        echo "    REAL DATA (CSV format)" >> ${pathogen}_data_summary.txt
        ls -la "./census_para.csv" >> ${pathogen}_data_summary.txt
    elif [ -f "./empty_census_para.csv" ] && [ -s "./empty_census_para.csv" ]; then
        echo "    *** ERROR: Invalid census file *** (Real census data is required)" >> ${pathogen}_data_summary.txt
        ls -la "./empty_census_para.csv" >> ${pathogen}_data_summary.txt
    else
        echo "    *** NO CENSUS FILE FOUND ***" >> ${pathogen}_data_summary.txt
    fi
    
    echo "" >> ${pathogen}_data_summary.txt
    echo "Analysis Parameters:" >> ${pathogen}_data_summary.txt
    echo "  Pathogen:        ${pathogen}" >> ${pathogen}_data_summary.txt
    echo "  Travel Types:    ${params.travel}" >> ${pathogen}_data_summary.txt
    echo "  CIDT Types:      ${params.cidt}" >> ${pathogen}_data_summary.txt
    echo "  MCMC Chains:     ${params.chains}" >> ${pathogen}_data_summary.txt
    echo "  Iterations:      ${params.iterations}" >> ${pathogen}_data_summary.txt
    echo "  Cores:           ${params.cores}" >> ${pathogen}_data_summary.txt
    echo "====================================================" >> ${pathogen}_data_summary.txt
    echo "" >> ${pathogen}_data_summary.txt
    
    # Copy to the log file as well
    cat ${pathogen}_data_summary.txt >> ${pathogen}_trendy.log
    
    # Find the R script
    SCRIPT_PATH="${bin_dir}/trendy.R"
    echo "Checking script path: \${SCRIPT_PATH}" >> ${pathogen}_trendy.log
    
    # Check that R script exists
    if [ ! -f "\${SCRIPT_PATH}" ]; then
        ls -la ${bin_dir} > script_dir_contents.txt
        error_exit "R script not found: \${SCRIPT_PATH}"
    fi
    
    # Copy progress tracking file if it exists
    PROGRESS_FILE="${bin_dir}/progress.R"
    if [ -f "\${PROGRESS_FILE}" ]; then
        echo "Copying progress tracking module..." >> ${pathogen}_trendy.log
        cp -v "\${PROGRESS_FILE}" ./ || echo "Warning: Could not copy progress tracking" >> ${pathogen}_trendy.log
    else
        echo "Progress tracking module not found (this is okay)" >> ${pathogen}_trendy.log
    fi
    
    # Run the pathogen analysis script
    echo "Running pathogen analysis for ${pathogen}" >> ${pathogen}_trendy.log
    echo "Using script: \${SCRIPT_PATH}" >> ${pathogen}_trendy.log
    
    # Run the R script with both stderr and stdout captured to a dedicated log file
    Rscript "\${SCRIPT_PATH}" \\
        --pathogen="${pathogen}" \\
        --mmwrFile="./${localMmwrFile}" \\
        \${CENSUS_B_ARG} \\
        \${CENSUS_P_ARG} \\
        \${PREPROC_ARG} \\
        --projID="${params.projID ?: 'foodnet'}" \\
        --travel="${params.travel}" \\
        --cidt="${params.cidt}" \\
        --outDir="./" \\
        --cores=${params.cores} \\
        --chains=${params.chains} \\
        --iterations=${params.iterations} \\
        --adapt_delta=${params.adapt_delta} \\
        --max_treedepth=${params.max_treedepth} \\
        --seed=${params.seed} \\
        ${params.debug ? '--debug' : ''} 2>&1 | tee ${pathogen}_R_output.log
    
    # Check return code from R script
    R_STATUS=\$?
    if [ \$R_STATUS -ne 0 ]; then
        echo "ERROR: R script failed with exit code \$R_STATUS" >> ${pathogen}_trendy.log
        
        # Create a more visible error summary
        echo "EXECUTION ERROR" >> ${pathogen}_data_summary.txt
        echo "----------------" >> ${pathogen}_data_summary.txt
        echo "The R script failed with exit code \$R_STATUS" >> ${pathogen}_data_summary.txt
        echo "Check ${pathogen}_R_output.log for details" >> ${pathogen}_data_summary.txt
        
        error_exit "Analysis failed for pathogen ${pathogen}"
    fi
    
    # Check for the presence of key output files
    if [ -f "${pathogen}_IRCatch.csv" ]; then
        echo "SUCCESS: Generated incidence rate file ${pathogen}_IRCatch.csv" >> ${pathogen}_data_summary.txt
        wc -l "${pathogen}_IRCatch.csv" >> ${pathogen}_data_summary.txt
    else
        echo "WARNING: No incidence rate file was generated!" >> ${pathogen}_data_summary.txt
    fi
    
    if [ -f "${pathogen}_brm.Rds" ]; then
        echo "SUCCESS: Generated model file ${pathogen}_brm.Rds" >> ${pathogen}_data_summary.txt
        ls -la "${pathogen}_brm.Rds" >> ${pathogen}_data_summary.txt
    else
        echo "WARNING: No model file was generated!" >> ${pathogen}_data_summary.txt
    fi
    
    # Check for figures
    PNG_COUNT=`ls -1 ${pathogen}_*.png 2>/dev/null | wc -l`
    if [ \$PNG_COUNT -gt 0 ]; then
        echo "SUCCESS: Generated \$PNG_COUNT visualization files" >> ${pathogen}_data_summary.txt
        ls -la ${pathogen}_*.png >> ${pathogen}_data_summary.txt
    else
        echo "WARNING: No visualization files were generated!" >> ${pathogen}_data_summary.txt
    fi
    
    echo "Analysis completed successfully" >> ${pathogen}_trendy.log
    echo "ANALYSIS COMPLETED SUCCESSFULLY" >> ${pathogen}_data_summary.txt
    
    # Copy the data summary to the log path as well to ensure it's published
    cp ${pathogen}_data_summary.txt ${pathogen}_data_summary.log
    """
}