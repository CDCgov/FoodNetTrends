#!/bin/bash

# Analysis workflow module for FoodNet Trends pipeline
# Mode 2: Run complete analysis with raw data files
# Version: 1.0 (2025-05)

handle_analysis_workflow() {
    echo ""
    echo -e "${BLUE}======== Input Files ========${NC}"
    
    # Default data files
    defaultMmwrFile="${DEFAULT_DATA_DIR}/mmwr9624_May2025.sas7bdat"
    defaultCensusFileB="${DEFAULT_DATA_DIR}/cen9624.sas7bdat"
    defaultCensusFileP="${DEFAULT_DATA_DIR}/cen9624_para.sas7bdat"
    
    read -p "MMWR data file [${defaultMmwrFile}]: " mmwrFile
    mmwrFile=${mmwrFile:-$defaultMmwrFile}
    
    read -p "Census file (bacterial) [${defaultCensusFileB}]: " censusFileB
    censusFileB=${censusFileB:-$defaultCensusFileB}
    
    read -p "Census file (parasitic) [${defaultCensusFileP}]: " censusFileP
    censusFileP=${censusFileP:-$defaultCensusFileP}
    
    # Validate that files exist
    for file in "${mmwrFile}" "${censusFileB}" "${censusFileP}"; do
        validate_file "$file" "Input" true
    done
    
    return 0
} 