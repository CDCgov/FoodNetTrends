#!/bin/bash
#==============================================================================
# FoodNetTrends Pipeline v1.0.0-rc.1 - Resume Previous Run
#==============================================================================
#
# This script allows users to easily resume a previous FoodNetTrends analysis
# by loading saved configurations and automatically applying the -resume flag.
#
# Usage:
#   ./resume_pipeline.sh
#
# Dependencies:
#   - Nextflow ≥ 24.10.4
#   - Singularity ≥ 4.1.4
#   - Previous run configuration files in ./.nextflow_configs/
#
# Last updated: 2025-05-27
#==============================================================================

# Initialize log file
error_log="foodnet_errors.log"
echo "$(date): Starting FoodNetTrends Resume Script v1.0.0-rc.1" >> "$error_log"

# Set up environment (same as run_pipeline.sh)
if [ -d "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/" ]; then
    DEFAULT_DATA_DIR="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data"
elif [ -d "/project/foodnet/data" ]; then
    DEFAULT_DATA_DIR="/project/foodnet/data"
else
    DEFAULT_DATA_DIR="./data"
    mkdir -p ./data
fi

# Load required modules if available
if command -v module >/dev/null 2>&1; then
    echo "Loading required modules..."
    module load singularity/3.10.0 2>/dev/null || module load singularity 2>/dev/null
    module load nextflow/24.10.4 2>/dev/null || module load nextflow 2>/dev/null
else
    echo "Warning: Module system not detected. Dependencies must be in PATH."
fi

# Set temporary directory
TMPDIR=/scicomp/scratch/$(whoami)
mkdir -p $TMPDIR/nextflow

echo "==========================================="
echo "   FoodNetTrends Pipeline Resume Tool      "
echo "           v1.0.0-rc.1                     "
echo "==========================================="
echo ""

# Function to find and display available configurations
find_configurations() {
    local config_dir="./.nextflow_configs"
    
    # Check if configs directory exists
    if [[ ! -d "$config_dir" ]]; then
        echo "Error: No saved configurations found."
        echo "The configuration directory '$config_dir' does not exist."
        echo ""
        echo "To create a resumable run, use ./run_pipeline.sh first."
        exit 1
    fi
    
    # Find all config files (excluding history)
    mapfile -t config_files < <(find "$config_dir" -maxdepth 1 -name "foodnet_config_*.sh" -type f 2>/dev/null | sort -r)
    
    if [[ ${#config_files[@]} -eq 0 ]]; then
        echo "Error: No saved configurations found in $config_dir"
        echo ""
        echo "To create a resumable run, use ./run_pipeline.sh first."
        exit 1
    fi
    
    echo "Available configurations to resume:"
    echo ""
    
    # Display available configurations with details
    for i in "${!config_files[@]}"; do
        config_file="${config_files[$i]}"
        # Extract project ID from filename
        proj_id=$(basename "$config_file" | sed 's/foodnet_config_//;s/\.sh$//')
        
        # Get configuration details
        if [[ -f "$config_file" ]]; then
            config_date=$(grep "# Generated" "$config_file" | cut -d' ' -f4-)
            config_pathogens=$(grep "export pathogens=" "$config_file" | cut -d'"' -f2)
            config_states=$(grep "export states=" "$config_file" | cut -d'"' -f2)
            config_chains=$(grep "export chains=" "$config_file" | cut -d'"' -f2)
            config_iterations=$(grep "export iterations=" "$config_file" | cut -d'"' -f2)
            
            # Check if Nextflow log exists for this project
            nextflow_status="Unknown"
            if [[ -f ".nextflow.log" ]]; then
                if grep -q "$proj_id" .nextflow.log 2>/dev/null; then
                    # Try to find completion status
                    if grep -q "Succeeded.*$proj_id" .nextflow.log 2>/dev/null; then
                        nextflow_status="Completed"
                    elif grep -q "WARN.*$proj_id" .nextflow.log 2>/dev/null; then
                        nextflow_status="Failed/Incomplete"
                    else
                        nextflow_status="In Progress"
                    fi
                fi
            fi
            
            echo "$((i+1))) Project: $proj_id"
            echo "   Date: $config_date"
            echo "   Status: $nextflow_status"
            echo "   Pathogens: $config_pathogens"
            echo "   States: ${config_states:-ALL}"
            echo "   Chains: $config_chains | Iterations: $config_iterations"
            echo ""
        fi
    done
    
    echo "0) Exit"
    echo ""
}

# Main script logic
echo "Searching for previous run configurations..."
echo ""

find_configurations

read -p "Select configuration to resume [1]: " config_choice
config_choice=${config_choice:-1}

if [[ "$config_choice" == "0" ]]; then
    echo "Exiting."
    exit 0
fi

if [[ "$config_choice" -ge 1 && "$config_choice" -le ${#config_files[@]} ]]; then
    selected_config="${config_files[$((config_choice-1))]}"
    
    echo ""
    echo "Loading configuration: $selected_config"
    
    # Source the configuration
    source "$selected_config"
    
    # Verify critical variables are set
    if [[ -z "$projID" || -z "$mmwrFile" || -z "$pathogens" ]]; then
        echo "Error: Configuration file is missing critical parameters."
        echo "Cannot resume this run."
        exit 1
    fi
    
    echo ""
    echo "================= RESUME CONFIGURATION ================="
    echo ""
    echo "Project ID: $projID"
    echo "Pathogens: $pathogens"
    echo "States: ${states:-ALL}"
    echo "Travel: $travel"
    echo "CIDT: $cidt"
    echo "Chains: $chains | Iterations: $iterations"
    echo "Output Directory: $outDir/$projID"
    echo ""
    echo "This will RESUME the analysis from the last successful step."
    echo "========================================================"
    echo ""
    
    read -p "Resume this analysis? (y/n) [y]: " confirm
    confirm=${confirm:-y}
    
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
        echo "Resume cancelled."
        exit 0
    fi
    
    # Build the nextflow command with -resume flag
    cmd="nextflow run main.nf -profile singularity,production -resume"
    
    # Add resource parameters
    cmd="$cmd -process.memory ${memory:-16.GB}"
    cmd="$cmd -process.cpus ${cores:-8}"
    cmd="$cmd -executor.queueSize 100"
    cmd="$cmd -executor.submitRateLimit '10/1min'"
    
    # Add all the analysis parameters
    cmd="$cmd --mmwrFile \"$mmwrFile\""
    
    if [[ -n "$censusFileB" && "$censusFileB" != "true" ]]; then
        cmd="$cmd --censusFileB \"$censusFileB\""
    fi
    
    if [[ -n "$censusFileP" && "$censusFileP" != "true" ]]; then
        cmd="$cmd --censusFileP \"$censusFileP\""
    fi
    
    cmd="$cmd --travel \"$travel\""
    cmd="$cmd --cidt \"$cidt\""
    cmd="$cmd --iterations $iterations"
    cmd="$cmd --chains $chains"
    cmd="$cmd --adapt_delta ${adapt_delta:-0.95}"
    cmd="$cmd --max_treedepth ${max_treedepth:-10}"
    cmd="$cmd --cores ${cores:-8}"
    cmd="$cmd --seed 123"
    cmd="$cmd --outdir \"$outDir\""
    cmd="$cmd --projID \"$projID\""
    cmd="$cmd --pathogen \"$pathogens\""
    cmd="$cmd --states \"$states\""
    
    # Add preprocessed data parameters if applicable
    if [[ -n "$preprocessed_data" ]]; then
        cmd="$cmd --preprocessed true --cleanFile \"$preprocessed_data\""
    fi
    
    if [[ -n "$preprocessed_metadata" ]]; then
        cmd="$cmd --metadata \"$preprocessed_metadata\""
    fi
    
    # Add dashboard parameters
    if [[ -n "$dashboard_params" ]]; then
        cmd="$cmd $dashboard_params"
    fi
    
    # Add serotype/serogroup parameters
    if [[ -n "$stec_serogroups" ]]; then
        cmd="$cmd --stec_serogroups \"$stec_serogroups\""
    fi
    if [[ -n "$salmonella_serotypes" ]]; then
        cmd="$cmd --salmonella_serotypes \"$salmonella_serotypes\""
    fi
    
    # Execute the command
    echo ""
    echo "Executing resume command..."
    echo ""
    echo "Command: $cmd"
    echo ""
    
    # Log the command
    echo "$(date): Resuming analysis with command: $cmd" >> "$error_log"
    
    # Execute
    eval "$cmd"
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "✓ Analysis resumed successfully!"
    else
        echo ""
        echo "✗ Error resuming analysis. Check .nextflow.log for details."
        exit 1
    fi
    
else
    echo "Invalid selection."
    exit 1
fi