#!/bin/bash

# Execution functions for FoodNet Trends pipeline
# This includes functions for executing the command and handling the result

# Execute the command and handle results
execute_command() {
    local final_cmd=$1
    local outDir=$2
    local background=$3
    local log_file=$4
    
    # Create the output directory if it doesn't exist
    mkdir -p "$outDir"
    echo "$(date): Created output directory: ${outDir}" >> "$error_log"
    
    # Run the command
    echo "$(date): Executing command: ${final_cmd}" >> "$error_log"
    if eval $final_cmd; then
        if [[ $background == true ]]; then
            # Check if background process started successfully
            sleep 2
            # Check for presence of process with specific parameters as signature
            if pgrep -f "nextflow.*${outDir}" > /dev/null; then
                echo -e "${GREEN}Process started in background. Check status with:${NC}"
                echo -e "${YELLOW}tail -f ${log_file}${NC}"
                echo "$(date): Background process started successfully" >> "$error_log"
                return 0
            else
                echo -e "${YELLOW}Warning: Background process may not have started correctly.${NC}"
                echo -e "${YELLOW}Check ${log_file} for details.${NC}"
                echo "$(date): Background process start verification failed" >> "$error_log"
                return 1
            fi
        else
            echo -e "${GREEN}Analysis completed successfully.${NC}"
            echo -e "${GREEN}Results are available in: $outDir${NC}"
            echo "$(date): Analysis completed successfully" >> "$error_log"
            return 0
        fi
    else
        echo -e "${RED}Error running analysis command.${NC}"
        echo -e "${RED}Check .nextflow.log for details.${NC}"
        echo "$(date): Command execution failed. Exit code: $?" >> "$error_log"
        return 1
    fi
}

# Function to get confirmation before executing
get_execution_confirmation() {
    read -p "Proceed with analysis? (y/n) [y]: " proceed
    proceed=${proceed:-y}
    
    if [[ "$proceed" =~ ^[Yy]$ ]]; then
        return 0
    else
        echo -e "${RED}Analysis cancelled.${NC}"
        echo "$(date): Analysis cancelled by user" >> "$error_log"
        return 1
    fi
} 