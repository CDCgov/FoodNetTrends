#!/bin/bash

# Utility functions for FoodNet Trends pipeline
# This includes file validation, JSON parsing, and logging functions
# Version: 1.0 (2025-05)

# Define validation function to reduce code duplication
validate_list() {
    local input_list=$1
    local valid_values=$2
    local item_type=$3
    local invalid_found=false
    
    IFS=',' read -ra INPUT_ARRAY <<< "$input_list"
    IFS=',' read -ra VALID_ARRAY <<< "$valid_values"
    
    for item in "${INPUT_ARRAY[@]}"; do
        valid=false
        for valid_item in "${VALID_ARRAY[@]}"; do
            if [[ "$item" == "$valid_item" ]]; then
                valid=true
                break
            fi
        done
        
        if [[ "$valid" == false ]]; then
            echo -e "${YELLOW}Warning: '$item' is not in the discovered $item_type list and may cause errors.${NC}"
            echo "$(date): Invalid $item_type: $item" >> "$error_log"
            invalid_found=true
        fi
    done
    
    if [[ "$invalid_found" == true ]]; then
        echo ""
        read -p "Continue anyway? (y/n) [n]: " continue_choice
        continue_choice=${continue_choice:-n}
        if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
            echo -e "${RED}Exiting.${NC}"
            echo "$(date): User canceled due to invalid $item_type" >> "$error_log"
            exit 1
        fi
    fi
}

# Validate that input files exist
validate_file() {
    local file_path=$1
    local file_type=$2
    local required=$3
    
    if [ ! -f "${file_path}" ]; then
        if [ "$required" = true ]; then
            echo -e "${RED}Error: $file_type file does not exist: ${file_path}${NC}"
            echo "$(date): Missing required $file_type file: ${file_path}" >> "$error_log"
            echo -e "${RED}Exiting.${NC}"
            exit 1
        else
            echo -e "${YELLOW}Warning: $file_type file does not exist: ${file_path}${NC}"
            echo "$(date): Missing $file_type file: ${file_path}" >> "$error_log"
            read -p "Continue anyway? (y/n) [n]: " continue_choice
            continue_choice=${continue_choice:-n}
            if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
                echo -e "${RED}Exiting.${NC}"
                exit 1
            fi
        fi
    fi
}

# Function for basic JSON parsing without jq
parse_json_value() {
    local json_file=$1
    local key=$2
    
    # Extract simple key-value pairs with grep and sed
    grep -o "\"$key\":[^,}]*" "$json_file" | sed 's/.*://' | sed 's/^[ \t]*//;s/[ \t]*$//' | sed 's/^"//;s/"$//'
}

# Function to extract array values from JSON without jq
parse_json_array() {
    local json_file=$1
    local key=$2
    local max_items=$3
    
    # Find the array in the JSON
    local array_text=$(grep -o "\"$key\":\[[^]]*\]" "$json_file")
    
    # Extract items from the array
    if [[ -n "$array_text" ]]; then
        # Remove the key and brackets
        array_text=${array_text#*\[}
        array_text=${array_text%\]*}
        
        # Split by comma and extract values
        local items=()
        local count=0
        
        # Split string by commas and extract values
        IFS=',' read -ra raw_items <<< "$array_text"
        for item in "${raw_items[@]}"; do
            # Clean up quotes and whitespace
            clean_item=$(echo "$item" | sed 's/^[ \t]*"//;s/"[ \t]*$//')
            items+=("$clean_item")
            ((count++))
            
            # Limit to max_items if specified
            if [[ -n "$max_items" && $count -ge $max_items ]]; then
                break
            fi
        done
        
        # Return as comma-separated list
        local result=$(IFS=,; echo "${items[*]}")
        echo "$result"
    else
        # Try alternative approach for complex nested arrays
        grep -o "\"[^\"]*\":[0-9]*" "$json_file" | 
            head -n ${max_items:-10} | 
            sed 's/"//g' | 
            cut -d':' -f1 | 
            paste -sd,
    fi
}

# Function to extract Salmonella serotypes from metadata JSON
get_serotypes() {
    local json_file=$1
    local count=${2:-5}  # Default to top 5 if not specified
    
    if [[ "$have_jq" == true ]]; then
        # Use jq for advanced extraction and sorting
        jq_cmd='.salmonella_serotypes | to_entries | sort_by(.value) | reverse | .[0:'$count'] | .[] | .key'
        jq -r "$jq_cmd" "$json_file" 2>/dev/null | paste -sd, || echo ""
    else
        # Fallback to grep/sed
        grep -o "\"[^\"]*\":[0-9]*" "$json_file" | 
            sort -t':' -k2,2nr | 
            head -n "$count" | 
            sed 's/"//g' | 
            cut -d':' -f1 | 
            paste -sd,
    fi
}

# Helper function to format file sizes
format_size() {
    local size=$1
    if [ "$size" -lt 1024 ]; then 
        echo "${size} B"
    elif [ "$size" -lt 1048576 ]; then 
        echo "$(( (size * 10 + 5) / 1024 / 10 )).$(( (size * 100) / 1024 % 10 )) KB"
    elif [ "$size" -lt 1073741824 ]; then 
        echo "$(( (size * 10 + 5) / 1048576 / 10 )).$(( (size * 100) / 1048576 % 10 )) MB"
    else 
        echo "$(( (size * 10 + 5) / 1073741824 / 10 )).$(( (size * 100) / 1073741824 % 10 )) GB"
    fi
}

# Helper function to check if metadata contains serotype information
has_serotypes() {
    local json_file=$1
    if [ -f "$json_file" ] && grep -q "salmonella_serotypes" "$json_file"; then
        return 0  # True in bash
    else
        return 1  # False in bash
    fi
} 