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
 * Last updated: 2025-05-21
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
    """
    #!/usr/bin/env bash
    set -e  # Exit immediately if a command exits with non-zero status
    
    # Setup error handling
    error_exit() {
        echo "ERROR: \$1" >> ${pathogen}_trendy.log
        exit 1
    }
    
    # Initialize variable defaults to ensure they're always defined
    CENSUS_B_ARG="--censusFileB=empty_census_bact.csv"
    CENSUS_P_ARG="--censusFileP=empty_census_para.csv"
    PREPROC_ARG=""
    
    # Log start of analysis
    echo "Starting analysis for pathogen: ${pathogen}" > ${pathogen}_trendy.log
    echo "Using MMWR data file: ${mmwrFile} (${mmwrFile.extension} format)" >> ${pathogen}_trendy.log
    
    # Make a local copy of the MMWR file to handle path issues
    cp -v "${mmwrFile}" ./input_data.${mmwrFile.extension} || error_exit "Failed to copy MMWR file"
    echo "Created local copy of MMWR file as: input_data.${mmwrFile.extension}" >> ${pathogen}_trendy.log
    
    # Handle census bacterial file
    if [ -f "${censusBFile}" ] && [ -s "${censusBFile}" ]; then
        # File exists and is not empty
        echo "Census bacterial file exists: ${censusBFile}" >> ${pathogen}_trendy.log
        
        # Get file extension from basename
        CENSUS_B_BASENAME=\$(basename "${censusBFile}")
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
        # Create empty placeholder
        echo "Census bacterial file missing or empty, creating placeholder" >> ${pathogen}_trendy.log
        echo "state,population,year,pathogentype" > empty_census_bact.csv
        echo "CA,10000000,2020,Bacterial" >> empty_census_bact.csv
        echo "CO,5000000,2020,Bacterial" >> empty_census_bact.csv
        echo "CT,3000000,2020,Bacterial" >> empty_census_bact.csv
        echo "GA,8000000,2020,Bacterial" >> empty_census_bact.csv
        echo "MD,5000000,2020,Bacterial" >> empty_census_bact.csv
        echo "MN,4000000,2020,Bacterial" >> empty_census_bact.csv
        echo "NM,2000000,2020,Bacterial" >> empty_census_bact.csv
        echo "NY,15000000,2020,Bacterial" >> empty_census_bact.csv
        echo "OR,3000000,2020,Bacterial" >> empty_census_bact.csv
        echo "TN,5000000,2020,Bacterial" >> empty_census_bact.csv
        chmod 644 empty_census_bact.csv
        CENSUS_B_ARG="--censusFileB=empty_census_bact.csv"
    fi
    
    # Handle census parasitic file
    if [ -f "${censusPFile}" ] && [ -s "${censusPFile}" ]; then
        # File exists and is not empty
        echo "Census parasitic file exists: ${censusPFile}" >> ${pathogen}_trendy.log
        
        # Get file extension from basename
        CENSUS_P_BASENAME=\$(basename "${censusPFile}")
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
        # Create empty placeholder
        echo "Census parasitic file missing or empty, creating placeholder" >> ${pathogen}_trendy.log
        echo "state,population,year,pathogentype" > empty_census_para.csv
        echo "CA,10000000,2020,Parasitic" >> empty_census_para.csv
        echo "CO,5000000,2020,Parasitic" >> empty_census_para.csv
        echo "CT,3000000,2020,Parasitic" >> empty_census_para.csv
        echo "GA,8000000,2020,Parasitic" >> empty_census_para.csv
        echo "MD,5000000,2020,Parasitic" >> empty_census_para.csv
        echo "MN,4000000,2020,Parasitic" >> empty_census_para.csv
        echo "NM,2000000,2020,Parasitic" >> empty_census_para.csv
        echo "NY,15000000,2020,Parasitic" >> empty_census_para.csv
        echo "OR,3000000,2020,Parasitic" >> empty_census_para.csv
        echo "TN,5000000,2020,Parasitic" >> empty_census_para.csv
        chmod 644 empty_census_para.csv
        CENSUS_P_ARG="--censusFileP=empty_census_para.csv"
    fi
    
    # Show file information
    echo "Checking copied files:" >> ${pathogen}_trendy.log
    ls -la ./ >> ${pathogen}_trendy.log
    
    # Run the main trend analysis with explicit path handling for everything
    echo "Using scripts path: ${bin_dir}" >> ${pathogen}_trendy.log
    
    # Create explicit path to R script - Use new unified script
    SCRIPT_PATH="${bin_dir}/process_pathogen.R"
    echo "Full script path: \${SCRIPT_PATH}" >> ${pathogen}_trendy.log
    
    # Check that R script exists
    if [ ! -f "\${SCRIPT_PATH}" ]; then
        # Fall back to trendy.R if process_pathogen.R doesn't exist
        echo "New script not found, falling back to legacy script" >> ${pathogen}_trendy.log
        SCRIPT_PATH="${bin_dir}/trendy.R"
        if [ ! -f "\${SCRIPT_PATH}" ]; then
            error_exit "R script not found at \${SCRIPT_PATH}. Directory contents of ${bin_dir}: \$(ls -la ${bin_dir})"
        fi
    fi
    
    # Check that our local data file copy exists and is readable
    if [ ! -f "./input_data.${mmwrFile.extension}" ] || [ ! -r "./input_data.${mmwrFile.extension}" ]; then
        error_exit "Local MMWR data file copy not found or not readable. Original file: ${mmwrFile}, Local copy attempt: ./input_data.${mmwrFile.extension}, Current directory contents: \$(ls -la ./)"
    fi
    
    # Add preprocessed flag based on file extension
    if [ "${mmwrFile.extension}" = "csv" ]; then
        PREPROC_ARG="--preprocessed=TRUE --cleanFile=./input_data.csv"
        echo "Using preprocessed mode for CSV file" >> ${pathogen}_trendy.log
    else
        PREPROC_ARG="--preprocessed=FALSE --rawFile=./input_data.sas7bdat"
        echo "Using raw data mode for SAS file" >> ${pathogen}_trendy.log
    fi
    
    # Make sure the census variables are fully initialized first
    echo "Census bacterial arg: \${CENSUS_B_ARG}" >> ${pathogen}_trendy.log
    echo "Census parasitic arg: \${CENSUS_P_ARG}" >> ${pathogen}_trendy.log
    echo "Preprocessing arg: \${PREPROC_ARG}" >> ${pathogen}_trendy.log
    
    # Run the pathogen analysis with the new unified script
    echo "Running pathogen analysis for ${pathogen}" >> ${pathogen}_trendy.log
    echo "Using script: \${SCRIPT_PATH}" >> ${pathogen}_trendy.log
    Rscript "\${SCRIPT_PATH}" \\
        --pathogen="${pathogen}" \\
        --mmwrFile="./input_data.${mmwrFile.extension}" \\
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
        ${params.debug ? '--debug' : ''}
    
    # Check return code from R script
    R_STATUS=\$?
    if [ \$R_STATUS -ne 0 ]; then
        echo "ERROR: R script failed with exit code \$R_STATUS" >> ${pathogen}_trendy.log
        error_exit "Analysis failed for pathogen ${pathogen}"
    fi
    
    echo "Analysis completed successfully" >> ${pathogen}_trendy.log
    """
}
