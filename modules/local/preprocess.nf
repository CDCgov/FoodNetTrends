/*
 * ==================================================================
 * FoodNet Trends - PREPROCESS Process Module
 * ==================================================================
 *
 * Purpose:
 *   This process handles the preprocessing of raw MMWR data files.
 *   It standardizes formats, cleans data, and generates optional
 *   metadata for downstream discovery and filtering.
 *
 * Inputs:
 *   - MMWR data file in SAS format
 *   - Census bacterial file in SAS format
 *   - Census parasitic file in SAS format
 *   - Output base name
 *   - Flag to generate metadata
 *
 * Outputs:
 *   - Cleaned CSV file with standardized format
 *   - Optional JSON metadata about dataset contents
 *   - Process logs for troubleshooting
 *
 * Error handling:
 *   - Input file validation
 *   - Output verification
 *   - Detailed logging
 *
 * Last updated: 2025-05-18
 * ==================================================================
 */

process PREPROCESS {
    tag "Preprocessing MMWR data"
    label 'process_medium'
    shell "/bin/bash"
    container 'foodnet.sif'

    // Organize outputs with better structure
    publishDir "${params.outdir}/preprocessed", mode: params.publish_dir_mode, saveAs: { filename ->
        if (filename.endsWith('.log')) {
            return "logs/$filename"
        } else {
            // Keep metadata files at the root level to match workflow expectations
            return filename
        }
    }

    input:
    path mmwrFile
    path censusFileB
    path censusFileP
    val outputBase
    val generateMetadata

    output:
    // Using consistent naming pattern across all outputs
    path "${outputBase}.csv", emit: cleanedData
    path "${outputBase}_metadata.json", optional: true, emit: metadata
    path "${outputBase}_*.log", emit: logs, optional: true

    script:
    """
    set -e  # Exit on error to prevent silent failures
    
    current_date=$(date)
    echo "Starting preprocessing at $current_date" | tee ${outputBase}_process.log
    echo "Input file: ${mmwrFile}" | tee -a ${outputBase}_process.log
    echo "Census file (bacterial): ${censusFileB}" | tee -a ${outputBase}_process.log
    echo "Census file (parasitic): ${censusFileP}" | tee -a ${outputBase}_process.log
    echo "Output base: ${outputBase}" | tee -a ${outputBase}_process.log
    echo "Generate metadata: ${generateMetadata}" | tee -a ${outputBase}_process.log
    
    # Additional validation for input file - prevents cryptic R errors later
    if [ ! -f "${mmwrFile}" ]; then
        echo "ERROR: Input file does not exist: ${mmwrFile}" > ${outputBase}_error.log
        exit 1
    fi
    
    # Define reusable placeholder function directly - avoid heredoc that causes syntax issues
    handle_census_file() {
        local file_path="\$1"
        local pathogen_type="\$2"
        local placeholder_file="\$3"
        local log_file="\$4"
        local output_path=""

        if [ ! -f "\${file_path}" ] || [ ! -s "\${file_path}" ]; then
            echo "WARNING: Census \${pathogen_type} file missing or empty: \${file_path}" >> \${log_file}
            echo "Creating standardized placeholder for census \${pathogen_type} file" >> \${log_file}
            
            # Create directory for placeholder file safely
            placeholder_dir=\${placeholder_file%/*}
            mkdir -p "\${placeholder_dir}"
            
            # Create standardized placeholder file
            echo "state,population,year,pathogentype" > "\${placeholder_file}"
            for state in CA CO CT GA MD MN NM NY OR TN; do
                for year in {2016..2023}; do
                    echo "\${state},5000000,\${year},\${pathogen_type}" >> "\${placeholder_file}"
                done
            done
            
            echo "CRITICAL WARNING: Using PLACEHOLDER \${pathogen_type} census data" | tee -a \${log_file}
            echo "                  Results will NOT be valid for production use!" | tee -a \${log_file}
            echo "                  Placeholder file created at: \${placeholder_file}" | tee -a \${log_file}
            
            output_path="\${placeholder_file}"
        else
            echo "Census \${pathogen_type} file exists: \${file_path}" >> \${log_file}
            output_path="\${file_path}"
            
            # Verify file format
            if [[ "\${file_path}" == *.csv ]]; then
                echo "Census \${pathogen_type} file format: CSV" >> \${log_file}
                # Verify file has required columns
                if ! head -1 "\${file_path}" | grep -i -q "state" || ! head -1 "\${file_path}" | grep -i -q "year"; then
                    echo "WARNING: Census \${pathogen_type} file may be missing required columns" >> \${log_file}
                fi
            elif [[ "\${file_path}" == *.sas7bdat ]]; then
                echo "Census \${pathogen_type} file format: SAS" >> \${log_file}
            else
                echo "WARNING: Census \${pathogen_type} file has unknown format: \${file_path}" >> \${log_file}
            fi
        fi
        
        echo "\${output_path}"
    }
    
    # Handle census bacterial file
    CENSUS_B=\$(handle_census_file "${censusFileB}" "bacterial" "${PWD}/placeholder_census_bacterial.csv" "${outputBase}_warnings.log")
    
    # Handle census parasitic file 
    CENSUS_P=\$(handle_census_file "${censusFileP}" "parasitic" "${PWD}/placeholder_census_parasitic.csv" "${outputBase}_warnings.log")
    
    # Execute the R preprocessing script with output capturing and proper argument handling
    echo "Using census bacterial file: \${CENSUS_B}" | tee -a ${outputBase}_process.log
    echo "Using census parasitic file: \${CENSUS_P}" | tee -a ${outputBase}_process.log
    
    # Note: tee command duplicates output to both console and log file
    Rscript ${workflow.projectDir}/bin/calcIR.R \\
      --mmwrFile="${mmwrFile}" \\
      --censusFileB="\${CENSUS_B}" \\
      --censusFileP="\${CENSUS_P}" \\
      --outputFile="${outputBase}.csv" \\
      --generate_metadata=${generateMetadata} \\
      2>&1 | tee ${outputBase}_R.log
    
    # Verify script created expected output before proceeding
    if [ ! -f "${outputBase}.csv" ]; then
        echo "ERROR: Expected output file not created: ${outputBase}.csv" > ${outputBase}_error.log
        exit 1
    fi
    
    end_date=$(date)
    echo "Preprocessing completed at $end_date" | tee -a ${outputBase}_process.log
    """
}