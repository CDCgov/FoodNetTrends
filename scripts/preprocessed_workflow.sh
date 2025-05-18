#!/bin/bash

# Preprocessed workflow module for FoodNet Trends pipeline
# Mode 3: Use existing preprocessed data for analysis
# Version: 1.0 (2025-05)

handle_preprocessed_workflow() {
    echo ""
    echo -e "${BLUE}======== Preprocessed Data ========${NC}"
    
    read -p "Path to preprocessed CSV file: " preprocessed_data
    
    # Validate path exists
    validate_file "${preprocessed_data}" "Preprocessed" true
    
    # Check for metadata file in same directory with _metadata.json suffix
    base_path=${preprocessed_data%.csv}
    auto_metadata="${base_path}_metadata.json"
    
    if [ -f "${auto_metadata}" ]; then
        # Verify metadata is valid JSON
        if ! jq '.' "${auto_metadata}" > /dev/null 2>&1; then
            echo -e "${YELLOW}Warning: Metadata file is not valid JSON: ${auto_metadata}${NC}"
            echo "$(date): Invalid JSON metadata: ${auto_metadata}" >> "$error_log"
            echo -e "${YELLOW}Continuing without metadata. Discovery will be limited.${NC}"
        else
            echo -e "${GREEN}Found valid metadata: ${auto_metadata}${NC}"
            preprocessed_metadata="${auto_metadata}"
        fi
    else
        echo -e "${YELLOW}Warning: No metadata file found with naming pattern ${base_path}_metadata.json${NC}"
        read -p "Path to metadata JSON file (leave empty to skip): " user_metadata
        
        if [ -n "$user_metadata" ]; then
            if [ -f "${user_metadata}" ]; then
                # Verify user-provided metadata is valid JSON
                if ! jq '.' "${user_metadata}" > /dev/null 2>&1; then
                    echo -e "${YELLOW}Warning: Metadata file is not valid JSON: ${user_metadata}${NC}"
                    echo "$(date): Invalid JSON metadata: ${user_metadata}" >> "$error_log"
                    echo -e "${YELLOW}Continuing without metadata. Discovery will be limited.${NC}"
                else
                    preprocessed_metadata="${user_metadata}"
                    echo -e "${GREEN}Using metadata: ${preprocessed_metadata}${NC}"
                fi
            else
                echo -e "${RED}Error: Metadata file does not exist: ${user_metadata}${NC}"
                echo "$(date): Missing metadata file: ${user_metadata}" >> "$error_log"
                echo -e "${YELLOW}Continuing without metadata. Discovery will be limited.${NC}"
            fi
        else
            echo -e "${YELLOW}Continuing without metadata. Discovery will be limited.${NC}"
        fi
    fi
    
    # Set mmwrFile to the preprocessed CSV
    mmwrFile=$preprocessed_data
    
    return 0
} 