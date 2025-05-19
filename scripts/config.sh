#!/bin/bash

# Configuration functions for FoodNet Trends pipeline
# This includes functions for handling MCMC parameters and command construction
# Version: 1.0 (2025-05)

# Function to get MCMC parameters for full mode
get_mcmc_params() {
    local params=()
    
    # Number of chains
    echo ""
    read -p "Number of chains [2]: " chains
    chains=${chains:-2}
    
    # Validate chains
    if ! [[ "$chains" =~ ^[0-9]+$ ]]; then
        echo "Error: Invalid input. Using default value (2)."
        echo "$(date): Invalid chains value: ${chains}" >> "$error_log"
        chains=2
    elif [ "$chains" -lt 2 ]; then
        echo "Warning: At least 2 chains recommended for convergence diagnostics."
        echo "Continuing with $chains chain(s)."
        echo "$(date): Using only ${chains} chain(s) (not recommended)" >> "$error_log"
    elif [ "$chains" -gt 8 ]; then
        echo "Warning: Large number of chains may significantly increase runtime."
        echo "Continuing with $chains chains."
        echo "$(date): Using high chain count: ${chains}" >> "$error_log"
    fi
    
    # Number of iterations
    echo ""
    read -p "Number of iterations [500]: " iterations
    iterations=${iterations:-500}
    
    # Validate iterations
    if ! [[ "$iterations" =~ ^[0-9]+$ ]]; then
        echo "Error: Invalid input. Using default value (500)."
        echo "$(date): Invalid iterations value: ${iterations}" >> "$error_log"
        iterations=500
    elif [ "$iterations" -lt 200 ]; then
        echo "Warning: Low iteration count may lead to poor convergence."
        echo "Continuing with $iterations iterations."
        echo "$(date): Low iteration count: ${iterations}" >> "$error_log"
    elif [ "$iterations" -gt 2000 ]; then
        echo "Warning: High iteration count will significantly increase runtime."
        echo "Continuing with $iterations iterations."
        echo "$(date): High iteration count: ${iterations}" >> "$error_log"
    fi
    
    # Adapt delta
    echo ""
    read -p "Adapt delta (0.0-1.0) [0.95]: " adapt_delta
    adapt_delta=${adapt_delta:-0.95}
    
    # Validate adapt delta without bc
    if [[ ! "$adapt_delta" =~ ^0?\.[0-9]+$ ]]; then
        echo "Error: Invalid input. Using default value (0.95)."
        echo "$(date): Invalid adapt_delta value: ${adapt_delta}" >> "$error_log"
        adapt_delta=0.95
    elif [[ "$adapt_delta" == "0.7"* ]] || [[ "$adapt_delta" == "0.6"* ]] || [[ "$adapt_delta" == "0.5"* ]] || [[ "$adapt_delta" == "0.4"* ]] || [[ "$adapt_delta" == "0.3"* ]] || [[ "$adapt_delta" == "0.2"* ]] || [[ "$adapt_delta" == "0.1"* ]] || [[ "$adapt_delta" == "0.0"* ]]; then
        echo "Warning: Low adapt_delta may cause algorithm issues."
        echo "Continuing with adapt_delta=$adapt_delta."
        echo "$(date): Low adapt_delta value: ${adapt_delta}" >> "$error_log"
    elif [[ "$adapt_delta" == "0.99"* ]] || [[ "$adapt_delta" == "1.0"* ]]; then
        echo "Warning: Very high adapt_delta may significantly increase runtime."
        echo "Continuing with adapt_delta=$adapt_delta."
        echo "$(date): High adapt_delta value: ${adapt_delta}" >> "$error_log"
    fi
    
    # Max treedepth
    echo ""
    read -p "Max treedepth [10]: " max_treedepth
    max_treedepth=${max_treedepth:-10}
    
    # Validate max treedepth
    if ! [[ "$max_treedepth" =~ ^[0-9]+$ ]]; then
        echo "Error: Invalid input. Using default value (10)."
        echo "$(date): Invalid max_treedepth: ${max_treedepth}" >> "$error_log"
        max_treedepth=10
    elif [ "$max_treedepth" -lt 8 ]; then
        echo "Warning: Low max_treedepth may cause truncated trajectories."
        echo "Continuing with max_treedepth=$max_treedepth."
        echo "$(date): Low max_treedepth value: ${max_treedepth}" >> "$error_log"
    elif [ "$max_treedepth" -gt 15 ]; then
        echo "Warning: High max_treedepth will significantly increase runtime."
        echo "Continuing with max_treedepth=$max_treedepth."
        echo "$(date): High max_treedepth value: ${max_treedepth}" >> "$error_log"
    fi
    
    # Return as an array
    params=("$chains" "$iterations" "$adapt_delta" "$max_treedepth")
    echo "${params[@]}"
}

# Function to get MCMC parameters for test mode
get_test_params() {
    local chains=1
    local iterations=100
    local adapt_delta=0.8
    local max_treedepth=8
    
    echo ""
    echo "Test Mode: Using minimal settings"
    echo "Chains: $chains"
    echo "Iterations: $iterations"
    echo "Adapt delta: $adapt_delta"
    echo "Max treedepth: $max_treedepth"
    
    # Return as an array
    local params=("$chains" "$iterations" "$adapt_delta" "$max_treedepth")
    echo "${params[@]}"
}

# Function to build the Nextflow command
build_command() {
    local resume=$1
    local censusFileB=$2
    local censusFileP=$3
    local travel=$4
    local cidt=$5
    local iterations=$6
    local chains=$7
    local adapt_delta=$8
    local max_treedepth=$9
    local outDir=${10}
    local pathogens=${11}
    local states=${12}
    local serotype_param=${13}
    
    # Trim whitespace from all path variables
    censusFileB=$(echo "$censusFileB" | xargs)
    censusFileP=$(echo "$censusFileP" | xargs)
    outDir=$(echo "$outDir" | xargs)
    
    local resume_flag=""
    if [[ "$resume" == "true" ]]; then
        resume_flag="-resume"
    fi
    
    cmd="nextflow run main.nf -profile singularity ${resume_flag} -entry SPLINE"
    cmd="${cmd} --censusFileB \"${censusFileB}\""
    cmd="${cmd} --censusFileP \"${censusFileP}\""
    cmd="${cmd} --travel \"${travel}\""
    cmd="${cmd} --cidt \"${cidt}\""
    cmd="${cmd} --iterations ${iterations}"
    cmd="${cmd} --chains ${chains}"
    cmd="${cmd} --adapt_delta ${adapt_delta}"
    cmd="${cmd} --max_treedepth ${max_treedepth}"
    cmd="${cmd} --seed 123"
    cmd="${cmd} --outdir \"${outDir}\""
    cmd="${cmd} --pathogen \"${pathogens}\""
    cmd="${cmd} --states \"${states}\""
    
    if [[ -n "$serotype_param" ]]; then
        cmd="${cmd} ${serotype_param}"
    fi
    
    echo "$cmd"
}

# Function to finalize command with additional parameters
finalize_command() {
    local cmd=$1
    local mmwrFile=$2
    local workflow_mode=$3
    local preprocessed_metadata=$4
    local background=$5
    local timestamp=$6
    
    # Trim whitespace from all path variables
    mmwrFile=$(echo "$mmwrFile" | xargs)
    preprocessed_metadata=$(echo "$preprocessed_metadata" | xargs)
    
    # Add mmwrFile parameter
    if [[ -n "$mmwrFile" ]]; then
        cmd="${cmd} --mmwrFile \"${mmwrFile}\""
    fi
    
    # Add preprocessed flag if using preprocessed data
    if [[ "$workflow_mode" == "1" || "$workflow_mode" == "3" ]]; then
        cmd="${cmd} --preprocessed true"
    fi
    
    # Add metadata parameter if available
    if [[ -n "$preprocessed_metadata" && -f "$preprocessed_metadata" ]]; then
        cmd="${cmd} --metadata \"${preprocessed_metadata}\""
    fi
    
    # Add background option if requested
    local final_cmd="$cmd"
    if [[ $background == true ]]; then
        # Check if output directory exists for log, create if needed
        logDir="logs"
        mkdir -p "$logDir"
        log_file="${logDir}/foodnet_run_${timestamp}.log"
        echo "$(date): Starting background process, log file at ${log_file}" >> "$error_log"
        bg_cmd="nohup ${cmd} > \"${log_file}\" 2>&1 &"
        final_cmd="${bg_cmd}"
        echo "Process will run in background with log: ${log_file}"
    fi
    
    echo "$final_cmd"
} 