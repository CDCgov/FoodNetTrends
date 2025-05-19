#!/bin/bash

# Execution functions for FoodNet Trends pipeline
# This includes functions for executing the command and handling the result

# Execute the command and handle results
execute_command() {
    local final_cmd=$1
    local outDir=$2
    local background=$3
    local log_file=$4
    
    # Trim any leading/trailing spaces from paths
    outDir=$(echo "$outDir" | xargs)
    
    # Create the output directory if it doesn't exist
    mkdir -p "$outDir"
    echo "$(date): Created output directory: ${outDir}" >> "$error_log"
    
    # Run the command
    echo "$(date): Executing command: ${final_cmd}" >> "$error_log"
    echo "Starting analysis..."
    if eval $final_cmd; then
        if [[ $background == true ]]; then
            # Check if background process started successfully
            sleep 2
            # Check for presence of process with specific parameters as signature
            if pgrep -f "nextflow.*${outDir}" > /dev/null; then
                echo "Process started in background. Check status with:"
                echo "tail -f ${log_file}"
                echo "$(date): Background process started successfully" >> "$error_log"
                return 0
            else
                echo "Warning: Background process may not have started correctly."
                echo "Check ${log_file} for details."
                echo "$(date): Background process start verification failed" >> "$error_log"
                return 1
            fi
        else
            echo "Analysis completed successfully."
            echo "Results are available in: $outDir"
            echo "$(date): Analysis completed successfully" >> "$error_log"
            return 0
        fi
    else
        echo "Error running analysis command."
        echo "Check .nextflow.log for details."
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
        echo "Analysis cancelled."
        echo "$(date): Analysis cancelled by user" >> "$error_log"
        return 1
    fi
} 