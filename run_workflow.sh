#!/bin/bash
#==============================================================================
# FoodNet Trends Pipeline v1.0 - Main Execution Script
#==============================================================================
#
# Purpose:
#   This interactive script serves as the main entry point for the FoodNet
#   Trends pipeline, guiding users through pipeline configuration and execution.
#   It provides a user-friendly interface to the underlying Nextflow workflow.
#
# Features:
#   - Interactive mode for user-guided parameter selection
#   - Support for three workflow modes: preprocess, full analysis, or using 
#     preprocessed data
#   - Parameter validation and error handling
#   - Options for background execution with logging
#   - Interactive dashboard generation for result visualization
#
# Usage:
#   ./run_workflow.sh
#
# Dependencies:
#   - Nextflow must be installed and available in PATH
#   - Properly configured Singularity container with all required R packages
#   - jq for JSON processing (will check and warn if not available)
#   - Properly configured module environment
#
# Output:
#   - Logs and results in user-specified output directory
#   - Background process logs in logs/ directory if running in background
#   - Interactive HTML dashboard for visualizing results
#
# Last updated: 2025-05-18
#==============================================================================

# Cleaner variables (no color codes)
GREEN=""
YELLOW=""
RED=""
BLUE=""
NC=""

# Check for jq but make it optional
have_jq=true
if ! command -v jq &> /dev/null; then
    echo "Warning: jq is not installed. Basic functionality will work, but advanced serotype filtering will be limited."
    have_jq=false
fi

# Set up paths, files
# Make data directory path environment-aware
if [ -d "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/" ]; then
  DEFAULT_DATA_DIR="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/"
elif [ -d "/project/foodnet/data" ]; then
  DEFAULT_DATA_DIR="/project/foodnet/data"
else
  DEFAULT_DATA_DIR="./data"
  # Create local data directory if it doesn't exist
  mkdir -p "$DEFAULT_DATA_DIR"
fi
outDir="output"  # Default to "output" directory in current location

# Initialize log file for error tracking
error_log="foodnet_errors.log"
echo "$(date): Starting FoodNet Trends Analysis Pipeline" > "$error_log"

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
            # Strip any trailing periods which may be causing confusion
            cleaned_item="${item%.}"
            # Check again with cleaned item
            for valid_item in "${VALID_ARRAY[@]}"; do
                if [[ "$cleaned_item" == "$valid_item" ]]; then
                    valid=true
                    break
                fi
            done
            
            if [[ "$valid" == false ]]; then
                echo "Warning: '$item' is not in the discovered $item_type list and may cause errors."
                echo "$(date): Invalid $item_type: $item" >> "$error_log"
                invalid_found=true
            fi
        fi
    done
    
    if [[ "$invalid_found" == true ]]; then
        echo ""
        read -p "Continue anyway? (y/n) [n]: " continue_choice
        continue_choice=${continue_choice:-n}
        if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
            echo "Exiting."
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
            echo "Error: $file_type file does not exist: ${file_path}"
            echo "$(date): Missing required $file_type file: ${file_path}" >> "$error_log"
            echo "Exiting."
            exit 1
        else
            echo "Warning: $file_type file does not exist: ${file_path}"
            echo "$(date): Missing $file_type file: ${file_path}" >> "$error_log"
            read -p "Continue anyway? (y/n) [n]: " continue_choice
            continue_choice=${continue_choice:-n}
            if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
                echo "Exiting."
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

# Set up modules if we're in an HPC environment
if command -v module &> /dev/null; then
    module purge
    module load nextflow/24.10.4
    module load singularity/4.1.4
    module load java/17.0.6
else
    echo "Warning: Module system not detected. Assuming dependencies are available in PATH."
    echo "$(date): Module system not detected" >> "$error_log"
fi

# Make sure TMPDIR is set
TMPDIR=${TMPDIR:-/scicomp/scratch/$(whoami)}
mkdir -p "$TMPDIR/nextflow" 2>/dev/null

# Create timestamp for automatic project ID
timestamp=$(date +%Y%m%d_%H%M%S)

# Display welcome banner
echo "Your Nextflow temporary/cache files will be placed in $TMPDIR/nextflow/ by default"
echo "========================================="
echo "   FoodNet Trends Analysis Pipeline      "
echo "========================================="
echo ""

# Ask for workflow mode
echo "Select mode:"
echo "1) Preprocess data (clean raw data files and generate metadata)"
echo "2) Run analysis (with complete pipeline)"
echo "3) Use existing preprocessed data"
read -p "Enter selection [1]: " workflow_mode
workflow_mode=${workflow_mode:-1}

# Handle preprocessed data
preprocessed_data=""
preprocessed_metadata=""

# Rest of the script follows...
# (This is just a starter to verify the menu displays properly)

# Debug output - initial script execution
echo "DEBUG: Starting run_workflow.sh script"

# Source all module scripts
SCRIPT_DIR="$(dirname "$0")/scripts"

# Debug output - script directory
echo "DEBUG: Script directory is $SCRIPT_DIR"

# Ensure scripts directory exists
if [ ! -d "$SCRIPT_DIR" ]; then
  echo "Error: Scripts directory not found at $SCRIPT_DIR"
  exit 1
fi

# Debug output - before sourcing modules
echo "DEBUG: About to source modules"

# Source all modules (order matters for some dependencies)
source "$SCRIPT_DIR/ui.sh"
echo "DEBUG: Loaded ui.sh"
source "$SCRIPT_DIR/utils.sh"
echo "DEBUG: Loaded utils.sh"
source "$SCRIPT_DIR/environment.sh"
echo "DEBUG: Loaded environment.sh"
source "$SCRIPT_DIR/config.sh"
echo "DEBUG: Loaded config.sh"
source "$SCRIPT_DIR/pathogen.sh"
echo "DEBUG: Loaded pathogen.sh"
source "$SCRIPT_DIR/preprocess_workflow.sh"
echo "DEBUG: Loaded preprocess_workflow.sh"
source "$SCRIPT_DIR/analysis_workflow.sh"
echo "DEBUG: Loaded analysis_workflow.sh"
source "$SCRIPT_DIR/preprocessed_workflow.sh"
echo "DEBUG: Loaded preprocessed_workflow.sh"
source "$SCRIPT_DIR/execution.sh"
echo "DEBUG: Loaded execution.sh"

# Debug output - calling main function
echo "DEBUG: Starting main function"

# Get dashboard preferences
get_dashboard_preference() {
  echo ""
  echo "======== Dashboard Generation ========"
  echo "An interactive HTML dashboard can be generated to visualize results."
  read -p "Generate interactive dashboard? (y/n) [y]: " enable_dashboard
  enable_dashboard=${enable_dashboard:-y}
  
  if [[ "$enable_dashboard" =~ ^[Yy]$ ]]; then
    # Get custom dashboard title
    read -p "Custom dashboard title [FoodNet Trends Analysis]: " dashboard_title
    dashboard_title=${dashboard_title:-"FoodNet Trends Analysis"}
    
    # Get logo path (optional)
    read -p "Path to logo image for dashboard (optional): " dashboard_logo
    
    # Return dashboard parameters
    echo "--enable_dashboard true --dashboard_title \"$dashboard_title\" ${dashboard_logo:+--dashboard_logo \"$dashboard_logo\"}"
  else
    echo "--enable_dashboard false"
  fi
}

# Update finalize_command function to include dashboard parameters
finalize_command() {
  local cmd="$1"
  local mmwrFile="$2"
  local workflow_mode="$3"
  local preprocessed_metadata="$4"
  local background="$5"
  local timestamp="$6"
  
  # Add mmwrFile parameter based on workflow mode
  if [[ "$workflow_mode" == "1" ]]; then
    # We're preprocessing and analyzing, use the cleaned file
    cmd="$cmd"
  elif [[ "$workflow_mode" == "2" ]]; then
    # Using raw data directly
    cmd="$cmd --mmwrFile \"$mmwrFile\""
  else
    # Using preprocessed data
    cmd="$cmd --mmwrFile \"$mmwrFile\" --preprocessed true"
    
    # Add metadata file if available
    if [[ -n "$preprocessed_metadata" ]]; then
      cmd="$cmd --metadata \"$preprocessed_metadata\""
    fi
  fi
  
  # Get dashboard parameters
  dashboard_params=$(get_dashboard_preference)
  cmd="$cmd $dashboard_params"
  
  # Add background option if needed
  if [[ "$background" == "yes" ]]; then
    # Create logs directory if it doesn't exist
    mkdir -p logs
    cmd="$cmd -bg"
  fi
  
  # Add proper error handling for resuming
  [[ "$flag" == "resume" ]] && cmd="$cmd -resume"
  
  # Return complete command
  echo "$cmd"
}

# Main function that orchestrates the entire workflow
main() {
  # Debug output - inside main function
  echo "DEBUG: Inside main function"

  # Initialize environment and paths
  echo "DEBUG: About to call setup_environment"
  setup_environment
  echo "DEBUG: Finished setup_environment"
  
  # Check for jq availability
  echo "DEBUG: About to call check_jq"
  check_jq
  echo "DEBUG: Finished check_jq"
  
  # Display welcome banner
  echo "DEBUG: About to call display_welcome"
  display_welcome
  echo "DEBUG: Finished display_welcome"
  
  # Get workflow mode from user
  echo "DEBUG: About to call get_workflow_mode"
  workflow_mode=$(get_workflow_mode)
  echo "DEBUG: get_workflow_mode returned: $workflow_mode"
  
  # Initialize variables
  preprocessed_data=""
  preprocessed_metadata=""
  censusFileB=""
  censusFileP=""
  mmwrFile=""
  
  # Handle workflow mode-specific initialization
  echo "DEBUG: Starting workflow mode-specific initialization for mode $workflow_mode"
  case "$workflow_mode" in
    1) 
      # Preprocess data
      echo "DEBUG: About to call handle_preprocess_workflow"
      handle_preprocess_workflow
      preprocess_result=$?
      echo "DEBUG: handle_preprocess_workflow returned: $preprocess_result"
      
      if [ $preprocess_result -eq 1 ]; then
        # Preprocessing failed
        echo "Preprocessing failed. Exiting."
        exit 1
      elif [ $preprocess_result -eq 0 ]; then
        # User chose not to continue to analysis
        exit 0
      fi
      # Otherwise preprocess_result is 2, continue to analysis
      ;;
    2) 
      # Run analysis with raw data
      echo "DEBUG: About to call handle_analysis_workflow"
      handle_analysis_workflow
      echo "DEBUG: Finished handle_analysis_workflow"
      ;;
    3) 
      # Use existing preprocessed data
      echo "DEBUG: About to call handle_preprocessed_workflow"
      handle_preprocessed_workflow
      echo "DEBUG: Finished handle_preprocessed_workflow"
      ;;
  esac
  
  # For modes 1 and 3, we need to ask for census files
  if [[ "$workflow_mode" == "1" || "$workflow_mode" == "3" ]]; then
    echo ""
    echo "======== Census Files ========"
    
    # Default census data files
    defaultCensusFileB="${DEFAULT_DATA_DIR}/cen9624.sas7bdat"
    defaultCensusFileP="${DEFAULT_DATA_DIR}/cen9624_para.sas7bdat"
    
    read -p "Census file (bacterial) [${defaultCensusFileB}]: " censusFileB
    censusFileB=${censusFileB:-$defaultCensusFileB}
    
    read -p "Census file (parasitic) [${defaultCensusFileP}]: " censusFileP
    censusFileP=${censusFileP:-$defaultCensusFileP}
    
    # Validate census files exist
    validate_file "${censusFileB}" "Census bacterial" true
    validate_file "${censusFileP}" "Census parasitic" true
  fi
  
  # Load metadata if available
  metadata_info=($(load_metadata "$preprocessed_metadata"))
  has_metadata=${metadata_info[0]}
  has_serotypes=${metadata_info[1]}
  ALL_PATHOGENS=${metadata_info[2]}
  ALL_STATES=${metadata_info[3]}
  DEFAULT_PATHOGENS=${metadata_info[4]}
  
  # Get travel status filter
  travel=$(get_travel_status)
  
  # Get CIDT/culture method filter
  cidt=$(get_cidt_method)
  
  # Get pathogen selection
  pathogens=$(select_pathogens "$ALL_PATHOGENS" "$DEFAULT_PATHOGENS")
  
  # Get state selection
  states=$(select_states "$ALL_STATES")
  
  # Get serotype selection if applicable
  serotype_param=$(select_serotypes "$pathogens" "$has_serotypes" "$preprocessed_metadata")
  
  # Get run mode
  flag=$(get_run_mode)
  
  # Get MCMC parameters based on run mode
  if [[ "$flag" == "test" ]]; then
    mcmc_params=($(get_test_params))
  elif [[ "$flag" == "full" ]]; then
    mcmc_params=($(get_mcmc_params))
  elif [[ "$flag" == "resume" ]]; then
    # Use default values for resumed runs
    mcmc_params=(2 500 0.95 10)
    echo ""
    echo "Resume Mode: Using parameters from previous run"
  fi
  
  # Extract MCMC parameters
  chains=${mcmc_params[0]}
  iterations=${mcmc_params[1]}
  adapt_delta=${mcmc_params[2]}
  max_treedepth=${mcmc_params[3]}
  
  # Get background execution preference
  background=$(get_background_choice)
  
  # Get output directory
  outDir=$(get_output_directory "$outDir")
  # Trim any leading/trailing spaces
  outDir=$(echo "$outDir" | xargs)
  
  # Build the command
  resume_flag="false"
  [[ "$flag" == "resume" ]] && resume_flag="true"
  cmd=$(build_command "$resume_flag" "$censusFileB" "$censusFileP" "$travel" "$cidt" "$iterations" \
         "$chains" "$adapt_delta" "$max_treedepth" "$outDir" "$pathogens" "$states" "$serotype_param")
  
  # Finalize command with additional parameters
  final_cmd=$(finalize_command "$cmd" "$mmwrFile" "$workflow_mode" "$preprocessed_metadata" "$background" "$timestamp")
  
  # Determine mode string for summary
  if [ "$flag" == "test" ]; then
    mode_string="Test run"
  elif [ "$flag" == "full" ]; then
    mode_string="Full analysis"
  else
    mode_string="Resume previous run"
  fi
  
  # Display summary
  display_summary "$mode_string" "$workflow_mode" "$pathogens" "$states" "$travel" "$cidt" \
                 "$serotype_param" "$mmwrFile" "$censusFileB" "$censusFileP" "$preprocessed_metadata" \
                 "$flag" "$chains" "$iterations" "$adapt_delta" "$max_treedepth" "$outDir" \
                 "$background" "$final_cmd"
  
  # Get confirmation from user
  if get_execution_confirmation; then
    echo "Starting analysis..."
    
    # Log file for background jobs
    log_file="logs/foodnet_run_${timestamp}.log"
    
    # Execute the command
    if ! execute_command "$final_cmd" "$outDir" "$background" "$log_file"; then
      echo "Error in execution."
      exit 1
    fi
  fi
}

# Run the main function
echo "DEBUG: Calling main function"
main "$@"
echo "DEBUG: Finished running main function" 