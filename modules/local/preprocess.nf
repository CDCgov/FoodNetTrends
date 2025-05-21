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
 * Last updated: 2025-05-21
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
    echo "FoodNet Trends Preprocessing" >> ${outputBase}_process.log
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
    
    # Handle census bacterial file
    if [ ! -f "${censusFileB}" ] || [ ! -s "${censusFileB}" ]; then
        echo "WARNING: Census bacterial file missing or empty" >> ${outputBase}_warnings.log
        echo "Creating placeholder bacterial census data" >> ${outputBase}_warnings.log
        
        # Create placeholder for bacterial census
        mkdir -p "${PWD}"
        echo "state,population,year,pathogentype" > "${PWD}/placeholder_census_bacterial.csv"
        echo "CA,10000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        echo "CO,5000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        echo "CT,3000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        echo "GA,8000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        echo "MD,5000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        echo "MN,4000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        echo "NM,2000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        echo "NY,15000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        echo "OR,3000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        echo "TN,5000000,2020,Bacterial" >> "${PWD}/placeholder_census_bacterial.csv"
        
        echo "CRITICAL WARNING: Using PLACEHOLDER bacterial census data" >> ${outputBase}_warnings.log
        echo "Results will NOT be valid for production use!" >> ${outputBase}_warnings.log
        
        CENSUS_B="${PWD}/placeholder_census_bacterial.csv"
    else
        echo "Using provided bacterial census data: ${censusFileB}" >> ${outputBase}_process.log
        CENSUS_B="${censusFileB}"
    fi
    
    # Handle census parasitic file
    if [ ! -f "${censusFileP}" ] || [ ! -s "${censusFileP}" ]; then
        echo "WARNING: Census parasitic file missing or empty" >> ${outputBase}_warnings.log
        echo "Creating placeholder parasitic census data" >> ${outputBase}_warnings.log
        
        # Create placeholder for parasitic census
        mkdir -p "${PWD}"
        echo "state,population,year,pathogentype" > "${PWD}/placeholder_census_parasitic.csv"
        echo "CA,10000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        echo "CO,5000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        echo "CT,3000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        echo "GA,8000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        echo "MD,5000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        echo "MN,4000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        echo "NM,2000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        echo "NY,15000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        echo "OR,3000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        echo "TN,5000000,2020,Parasitic" >> "${PWD}/placeholder_census_parasitic.csv"
        
        echo "CRITICAL WARNING: Using PLACEHOLDER parasitic census data" >> ${outputBase}_warnings.log
        echo "Results will NOT be valid for production use!" >> ${outputBase}_warnings.log
        
        CENSUS_P="${PWD}/placeholder_census_parasitic.csv"
    else
        echo "Using provided parasitic census data: ${censusFileP}" >> ${outputBase}_process.log
        CENSUS_P="${censusFileP}"
    fi
    
    # Execute the R preprocessing script
    echo "Running preprocessing script with:" >> ${outputBase}_process.log
    echo "  MMWR file: ${mmwrFile}" >> ${outputBase}_process.log
    echo "  Census bacterial file: \$CENSUS_B" >> ${outputBase}_process.log
    echo "  Census parasitic file: \$CENSUS_P" >> ${outputBase}_process.log
    echo "  Output file: ${outputBase}.csv" >> ${outputBase}_process.log
    echo "  Generate metadata: ${generateMetadata}" >> ${outputBase}_process.log
    
    Rscript ${workflow.projectDir}/bin/calcIR.R \\
      --mmwrFile="${mmwrFile}" \\
      --censusFileB="\$CENSUS_B" \\
      --censusFileP="\$CENSUS_P" \\
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