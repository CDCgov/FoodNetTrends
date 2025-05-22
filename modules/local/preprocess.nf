/*
 * ==================================================================
 * FoodNetTrends v1.0.0-rc.1 - Data Preprocessing Module
 * ==================================================================
 * 
 * Purpose:
 *   Standardizes raw MMWR surveillance data for downstream analysis.
 *   Performs data cleaning, format validation, and metadata generation
 *   to ensure consistent input for Bayesian modeling processes.
 *
 * Inputs:
 *   - Raw MMWR data file (SAS format)
 *   - Census files (bacterial and parasitic populations)
 *   - Output naming parameters
 *   - Metadata generation flags
 * 
 * Outputs:
 *   - Standardized CSV data file
 *   - JSON metadata with dataset characteristics
 *   - Process execution logs
 * 
 * Last updated: 2025-05-22
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
    # Create process log
    echo "===============================================" > ${outputBase}_process.log
    echo "FoodNetTrends Preprocessing" >> ${outputBase}_process.log
    echo "===============================================" >> ${outputBase}_process.log
    echo "Starting preprocessing" >> ${outputBase}_process.log
    echo "Input file: ${mmwrFile}" >> ${outputBase}_process.log
    echo "Census file (bacterial): ${censusFileB}" >> ${outputBase}_process.log
    echo "Census file (parasitic): ${censusFileP}" >> ${outputBase}_process.log
    echo "Output base: ${outputBase}" >> ${outputBase}_process.log
    echo "Generate metadata: ${generateMetadata}" >> ${outputBase}_process.log
    
    # Create warnings log
    touch ${outputBase}_warnings.log
    
    # Validate MMWR file
    if [ ! -f "${mmwrFile}" ]; then
        echo "ERROR: Input file does not exist: ${mmwrFile}" > ${outputBase}_error.log
        exit 1
    fi
    
    # Handle census bacterial file - get absolute path
    if [ ! -f "${censusFileB}" ] || [ ! -s "${censusFileB}" ]; then
        echo "ERROR: Census bacterial file missing or empty" >> ${outputBase}_error.log
        echo "Census bacterial file is required for analysis" >> ${outputBase}_error.log
        exit 1
    else
        echo "Using provided bacterial census data: ${censusFileB}" >> ${outputBase}_process.log
        CENSUS_B_ABS=\$(readlink -f "${censusFileB}")
        echo "Absolute path to bacterial census: \$CENSUS_B_ABS" >> ${outputBase}_process.log
    fi
    
    # Handle census parasitic file - get absolute path
    if [ ! -f "${censusFileP}" ] || [ ! -s "${censusFileP}" ]; then
        echo "ERROR: Census parasitic file missing or empty" >> ${outputBase}_error.log
        echo "Census parasitic file is required for analysis" >> ${outputBase}_error.log
        exit 1
    else
        echo "Using provided parasitic census data: ${censusFileP}" >> ${outputBase}_process.log
        CENSUS_P_ABS=\$(readlink -f "${censusFileP}")
        echo "Absolute path to parasitic census: \$CENSUS_P_ABS" >> ${outputBase}_process.log
    fi
    
    # Create symbolic links with stable names for downstream steps
    ln -sf "\$CENSUS_B_ABS" census_bacterial.csv
    ln -sf "\$CENSUS_P_ABS" census_parasitic.csv
    
    # Execute the R preprocessing script with absolute paths
    echo "Running preprocessing script with:" >> ${outputBase}_process.log
    echo "  MMWR file: ${mmwrFile}" >> ${outputBase}_process.log
    echo "  Census bacterial file: \$CENSUS_B_ABS" >> ${outputBase}_process.log
    echo "  Census parasitic file: \$CENSUS_P_ABS" >> ${outputBase}_process.log
    echo "  Output file: ${outputBase}.csv" >> ${outputBase}_process.log
    echo "  Generate metadata: ${generateMetadata}" >> ${outputBase}_process.log
    
    Rscript ${workflow.projectDir}/bin/preprocess.R \\
      --mmwrFile="${mmwrFile}" \\
      --censusFileB="\$CENSUS_B_ABS" \\
      --censusFileP="\$CENSUS_P_ABS" \\
      --outputFile="${outputBase}.csv" \\
      --generate_metadata=${generateMetadata} \\
      2>&1 | tee ${outputBase}_R.log
    
    # Verify output was created
    if [ ! -f "${outputBase}.csv" ]; then
        echo "ERROR: Expected output file not created: ${outputBase}.csv" > ${outputBase}_error.log
        exit 1
    fi
    
    echo "Preprocessing completed successfully" >> ${outputBase}_process.log
    echo "===============================================" >> ${outputBase}_process.log
    """
}