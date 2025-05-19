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
    
    echo "Starting preprocessing at \$(date)" | tee ${outputBase}_process.log
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
    if [ ! -f "${censusFileB}" ]; then
        echo "WARNING: Census bacterial file does not exist: ${censusFileB}" >> ${outputBase}_warnings.log
        echo "Creating empty placeholder for census bacterial file" >> ${outputBase}_warnings.log
        touch empty_census_bacterial.csv
        censusFileB="empty_census_bacterial.csv"
    fi
    if [ ! -f "${censusFileP}" ]; then
        echo "WARNING: Census parasitic file does not exist: ${censusFileP}" >> ${outputBase}_warnings.log
        echo "Creating empty placeholder for census parasitic file" >> ${outputBase}_warnings.log
        touch empty_census_parasitic.csv
        censusFileP="empty_census_parasitic.csv"
    fi
    
    # Execute the R preprocessing script with output capturing
    # Note: tee command duplicates output to both console and log file
    Rscript ${workflow.projectDir}/bin/calcIR.R \
      --mmwrFile ${mmwrFile} \
      --censusFileB "${censusFileB}" \
      --censusFileP "${censusFileP}" \
      --outputFile ${outputBase}.csv \
      --generate_metadata ${generateMetadata} \
      2>&1 | tee ${outputBase}_R.log
    
    # Verify script created expected output before proceeding
    if [ ! -f "${outputBase}.csv" ]; then
        echo "ERROR: Expected output file not created: ${outputBase}.csv" > ${outputBase}_error.log
        exit 1
    fi
    
    echo "Preprocessing completed at \$(date)" | tee -a ${outputBase}_process.log
    """
}
