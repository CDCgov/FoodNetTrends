#!/bin/bash

# Pathogen and state selection functions for FoodNet Trends pipeline
# Version: 1.0 (2025-05)

# Load metadata from file
load_metadata() {
    local preprocessed_metadata=$1
    local has_metadata=false
    local ALL_PATHOGENS="CAMPYLOBACTER,CYCLOSPORA,SALMONELLA,SHIGELLA,STEC,VIBRIO,YERSINIA"
    local ALL_STATES="CA,CO,CT,GA,MD,MN,NM,NY,OR,TN"
    local DEFAULT_PATHOGENS="CAMPYLOBACTER,CYCLOSPORA"
    local has_serotypes=false
    local metadata_values=()
    
    if [[ -n "${preprocessed_metadata}" && -f "${preprocessed_metadata}" ]]; then
        echo "Found metadata file: ${preprocessed_metadata}"
        echo "Metadata will be used by the workflow directly."
        
        # Set flag to indicate metadata is available
        has_metadata=true
        
        # Check if file contains serotypes
        if has_serotypes "${preprocessed_metadata}"; then
            has_serotypes=true
            echo "Metadata contains Salmonella serotype information."
        fi
        
        # Extract key information from metadata file
        echo "Extracting key information from metadata:"
        
        # Extract record count using basic parsing
        record_count=$(parse_json_value "${preprocessed_metadata}" "record_count")
        if [ -n "$record_count" ]; then
            echo "Dataset contains approximately $record_count records"
        fi
        
        # Try to extract and display pathogens from metadata
        if [[ "$have_jq" == true ]]; then
            # Use jq if available
            metadata_pathogens=$(jq -r '.pathogens | join(",")' "${preprocessed_metadata}" 2>/dev/null)
        else
            # Fallback to basic parsing
            metadata_pathogens=$(parse_json_array "${preprocessed_metadata}" "pathogens")
        fi
        
        if [ -n "$metadata_pathogens" ]; then
            echo "Pathogens in dataset: $metadata_pathogens"
            # Update ALL_PATHOGENS if we found valid ones in metadata
            ALL_PATHOGENS="$metadata_pathogens"
        fi
        
        # Try to extract and display states from metadata
        if [[ "$have_jq" == true ]]; then
            # Use jq if available
            metadata_states=$(jq -r '.states | join(",")' "${preprocessed_metadata}" 2>/dev/null)
        else
            # Fallback to basic parsing
            metadata_states=$(parse_json_array "${preprocessed_metadata}" "states")
        fi
        
        if [ -n "$metadata_states" ]; then
            echo "States in dataset: $metadata_states"
            # Update ALL_STATES if we found valid ones in metadata
            ALL_STATES="$metadata_states"
        fi
    else
        echo "Warning: No metadata file found. Using default values."
        has_metadata=false
    fi
    
    # Package values into an array
    metadata_values=("$has_metadata" "$has_serotypes" "$ALL_PATHOGENS" "$ALL_STATES" "$DEFAULT_PATHOGENS")
    echo "${metadata_values[@]}"
}

# Get pathogen selection from user
select_pathogens() {
    local all_pathogens=$1
    local default_pathogens=$2
    
    echo ""
    echo "======== Pathogen Selection ========"
    echo "Available pathogens in this dataset:"
    
    # Parse the comma-separated list and display each pathogen
    IFS=',' read -ra PATHOGEN_ARRAY <<< "$all_pathogens"
    for p in "${PATHOGEN_ARRAY[@]}"; do
        echo "- $p"
    done
    echo ""
    
    echo "1) Run ALL available pathogens"
    echo "2) Select specific pathogens"
    read -p "Enter selection [2]: " pathogen_mode
    pathogen_mode=${pathogen_mode:-2}
    
    if [[ "$pathogen_mode" == "1" ]]; then
        # Use all pathogens
        echo "Selected: ALL pathogens (${all_pathogens})"
        echo "$all_pathogens"
    else
        # Ask for specific pathogens
        echo ""
        echo "Enter pathogens to analyze (comma-separated with NO spaces)"
        read -p "Leave blank for default (${default_pathogens}): " pathogens
        pathogens=${pathogens:-"$default_pathogens"}
        
        # Use the validation function
        validate_list "$pathogens" "$all_pathogens" "pathogen"
        
        echo "$pathogens"
    fi
}

# Get state selection from user
select_states() {
    local all_states=$1
    
    echo ""
    echo "======== State Selection ========"
    echo "Available states in this dataset:"
    
    # Parse the comma-separated list and display each state
    IFS=',' read -ra STATE_ARRAY <<< "$all_states"
    for s in "${STATE_ARRAY[@]}"; do
        echo "- $s"
    done
    echo ""
    
    echo "1) Use ALL available states"
    echo "2) Select specific states"
    read -p "Enter selection [1]: " state_mode
    state_mode=${state_mode:-1}
    
    if [[ "$state_mode" == "1" ]]; then
        # Use all states
        echo "Selected: ALL states (${all_states})"
        echo "$all_states"
    else
        # Ask for specific states
        echo ""
        echo "Enter states to analyze (comma-separated with NO spaces)"
        read -p "Leave blank for all states: " states
        states=${states:-"$all_states"}
        
        # Use the validation function
        validate_list "$states" "$all_states" "state"
        
        echo "$states"
    fi
}

# Get Salmonella serotype selection if applicable
select_serotypes() {
    local pathogens=$1
    local has_serotypes=$2
    local preprocessed_metadata=$3
    local serotype_param=""
    
    if [[ ",$pathogens," == *",SALMONELLA,"* ]] && [[ "$has_serotypes" == "true" ]] && [[ -n "$preprocessed_metadata" ]]; then
        # Check if serotypes are available
        echo ""
        echo "======== Salmonella Serotype Analysis ========"
        echo "Salmonella was selected. Do you want to:"
        echo "1) Analyze ALL Salmonella serotypes together"
        echo "2) Analyze specific serotypes separately"
        echo "3) Focus on a single top serotype only"
        read -p "Enter selection [1]: " serotype_mode
        serotype_mode=${serotype_mode:-1}
        
        if [[ "$serotype_mode" == "2" ]]; then
            # Display available serotypes (top 20 to keep it manageable)
            echo ""
            echo "Top Salmonella serotypes in dataset:"
            
            if [[ "$have_jq" == true ]]; then
                # Use jq for more advanced display if available
                jq_cmd='.salmonella_serotypes | to_entries | sort_by(.value) | reverse | .[0:20] | .[] | "\(.key): \(.value) isolates"'
                top_serotypes=$(jq -r "$jq_cmd" "${preprocessed_metadata}" 2>/dev/null)
                if [[ -n "$top_serotypes" ]]; then
                    echo "$top_serotypes" | sed 's/^/- /'
                else
                    echo "Warning: Could not process serotype information with jq."
                    # Fallback to basic parsing
                    basic_serotypes=$(grep -o "\"[^\"]*\":[0-9]*" "${preprocessed_metadata}" | sort -t':' -k2,2nr | head -n 20)
                    echo "$basic_serotypes" | sed 's/"//g' | sed 's/:/: /' | sed 's/^/- /'
                fi
            else
                # Use basic parsing without jq
                basic_serotypes=$(grep -o "\"[^\"]*\":[0-9]*" "${preprocessed_metadata}" | sort -t':' -k2,2nr | head -n 20)
                echo "$basic_serotypes" | sed 's/"//g' | sed 's/:/: /' | sed 's/^/- /'
            fi
            
            echo ""
            echo "Enter serotypes to analyze (comma-separated with NO spaces)"
            read -p "Serotypes: " serotypes
            
            # Validate serotypes
            if [[ -n "$serotypes" ]]; then
                # Add parameter for serotypes
                serotype_param="--salmonella_serotypes \"$serotypes\""
            else
                echo "Warning: No serotypes specified. Analyzing all Salmonella together."
                serotype_param=""
            fi
        elif [[ "$serotype_mode" == "3" ]]; then
            # Get top serotype 
            if [[ "$have_jq" == true ]]; then
                # Use jq if available
                top_serotype=$(jq -r '.salmonella_serotypes | to_entries | sort_by(.value) | reverse | .[0].key' "${preprocessed_metadata}" 2>/dev/null)
            else
                # Fallback to basic parsing - get the first serotype with highest count
                top_serotype=$(grep -o "\"[^\"]*\":[0-9]*" "${preprocessed_metadata}" | sort -t':' -k2,2nr | head -n 1 | sed 's/"//g' | cut -d':' -f1)
            fi
            
            if [[ -n "$top_serotype" ]]; then
                echo "Will focus on top serotype: $top_serotype"
                serotype_param="--salmonella_serotypes \"$top_serotype\""
            else
                echo "Warning: Could not determine top serotype. Analyzing all Salmonella together."
                echo "$(date): Failed to extract top serotype" >> "$error_log"
                serotype_param=""
            fi
        else
            # Analyze all serotypes together (default)
            serotype_param=""
        fi
    fi
    
    echo "$serotype_param"
} 