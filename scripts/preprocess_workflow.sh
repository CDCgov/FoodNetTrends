#!/bin/bash

# Preprocess workflow module for FoodNet Trends pipeline
# Mode 1: Preprocess data and optionally continue to analysis
# Version: 1.0 (2025-05)

handle_preprocess_workflow() {
    echo ""
    echo -e "${BLUE}======== Input Files ========${NC}"
    
    # Default data file
    defaultMmwrFile="${DEFAULT_DATA_DIR}/mmwr9624_May2025.sas7bdat"
    
    read -p "MMWR data file [${defaultMmwrFile}]: " mmwrFile
    mmwrFile=${mmwrFile:-$defaultMmwrFile}
    
    # Validate file exists
    validate_file "${mmwrFile}" "MMWR" true
    
    # Set output location
    echo ""
    echo -e "${BLUE}======== Output Settings ========${NC}"
    
    # Default output directory for preprocessed data
    defaultPreprocessedDir="preprocessed_$(date +%Y%m%d_%H%M%S)"
    read -p "Output directory for preprocessed data [${defaultPreprocessedDir}]: " preprocessedDir
    preprocessedDir=${preprocessedDir:-$defaultPreprocessedDir}
    
    # Default output base name (derived from input filename)
    fileBasename=$(basename "${mmwrFile}" | sed 's/\.[^.]*$//')
    defaultOutputBase="foodnet_data_${fileBasename}"
    read -p "Base name for output files [${defaultOutputBase}]: " outputBase
    outputBase=${outputBase:-$defaultOutputBase}
    
    # Ask about metadata generation
    echo ""
    read -p "Generate metadata JSON? (y/n) [y]: " generate_metadata
    generate_metadata=${generate_metadata:-y}
    if [[ "$generate_metadata" =~ ^[Yy]$ ]]; then
        metadata_param="--generateMetadata true"
    else
        metadata_param="--generateMetadata false"
    fi
    
    echo -e "${GREEN}Running preprocessing...${NC}"
    
    # Run the preprocessing workflow
    preprocess_cmd="nextflow run main.nf -profile singularity -entry PREPROCESS_WORKFLOW \
      --mmwrFile \"${mmwrFile}\" \
      --outdir \"${preprocessedDir}\" \
      --outputBase \"${outputBase}\" \
      ${metadata_param}"
    
    echo "$(date): Running preprocessing command: ${preprocess_cmd}" >> "$error_log"
    
    if ! eval $preprocess_cmd; then
        echo -e "${RED}Error: Preprocessing failed.${NC}"
        echo -e "${RED}Check .nextflow.log for details.${NC}"
        echo "$(date): Preprocessing failed, check .nextflow.log" >> "$error_log"
        return 1
    fi
    
    # Set paths to preprocessed data and metadata
    preprocessed_data="${preprocessedDir}/preprocessed/${outputBase}.csv"
    preprocessed_metadata="${preprocessedDir}/preprocessed/${outputBase}_metadata.json"
    
    # Check if files were created
    if [ ! -f "${preprocessed_data}" ]; then
        echo -e "${RED}Error: Preprocessing failed to create expected CSV file: ${preprocessed_data}${NC}"
        echo "$(date): Missing expected output CSV: ${preprocessed_data}" >> "$error_log"
        echo -e "${RED}Check logs for details.${NC}"
        return 1
    fi
    
    echo -e "${GREEN}Preprocessing complete!${NC}"
    echo -e "- Cleaned data: ${GREEN}${preprocessed_data}${NC}"
    
    if [[ "$generate_metadata" =~ ^[Yy]$ ]] && [ -f "${preprocessed_metadata}" ]; then
        # Verify metadata is valid JSON
        if ! jq '.' "${preprocessed_metadata}" > /dev/null 2>&1; then
            echo -e "${YELLOW}Warning: Generated metadata file is not valid JSON: ${preprocessed_metadata}${NC}"
            echo "$(date): Invalid JSON metadata: ${preprocessed_metadata}" >> "$error_log"
        else
            echo -e "- Metadata: ${GREEN}${preprocessed_metadata}${NC}"
        fi
    fi
    
    echo ""
    read -p "Proceed to analysis with this preprocessed data? (y/n) [y]: " proceed_to_analysis
    proceed_to_analysis=${proceed_to_analysis:-y}
    
    if [[ ! "$proceed_to_analysis" =~ ^[Yy]$ ]]; then
        echo -e "${GREEN}Preprocessing complete. You can run analysis later using mode 3.${NC}"
        return 0
    fi
    
    # Continue to analysis using the preprocessed data
    echo -e "${GREEN}Proceeding to analysis...${NC}"
    
    # Set mmwrFile to the preprocessed CSV for downstream use
    mmwrFile=$preprocessed_data
    
    # Continue with census files, pathogens, and other settings
    return 2  # Special return code to indicate continuing to analysis
} 