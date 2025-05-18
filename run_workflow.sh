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

# Source all module scripts
SCRIPT_DIR="$(dirname "$0")/scripts"

# Ensure scripts directory exists
if [ ! -d "$SCRIPT_DIR" ]; then
  echo "Error: Scripts directory not found at $SCRIPT_DIR"
  exit 1
fi

# Source all modules (order matters for some dependencies)
source "$SCRIPT_DIR/ui.sh"
source "$SCRIPT_DIR/utils.sh"
source "$SCRIPT_DIR/environment.sh"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/pathogen.sh"
source "$SCRIPT_DIR/preprocess_workflow.sh"
source "$SCRIPT_DIR/analysis_workflow.sh"
source "$SCRIPT_DIR/preprocessed_workflow.sh"
source "$SCRIPT_DIR/execution.sh"

# Get dashboard preferences
get_dashboard_preference() {
  echo ""
  echo -e "${BLUE}======== Dashboard Generation ========${NC}"
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
  # Initialize environment and paths
  setup_environment
  
  # Check for jq availability
  check_jq
  
  # Display welcome banner
  display_welcome
  
  # Get workflow mode from user
  workflow_mode=$(get_workflow_mode)
  
  # Initialize variables
  preprocessed_data=""
  preprocessed_metadata=""
  censusFileB=""
  censusFileP=""
  mmwrFile=""
  
  # Handle workflow mode-specific initialization
  case "$workflow_mode" in
    1) 
      # Preprocess data
      handle_preprocess_workflow
      preprocess_result=$?
      
      if [ $preprocess_result -eq 1 ]; then
        # Preprocessing failed
        echo -e "${RED}Preprocessing failed. Exiting.${NC}"
        exit 1
      elif [ $preprocess_result -eq 0 ]; then
        # User chose not to continue to analysis
        exit 0
      fi
      # Otherwise preprocess_result is 2, continue to analysis
      ;;
    2) 
      # Run analysis with raw data
      handle_analysis_workflow
      ;;
    3) 
      # Use existing preprocessed data
      handle_preprocessed_workflow
      ;;
  esac
  
  # For modes 1 and 3, we need to ask for census files
  if [[ "$workflow_mode" == "1" || "$workflow_mode" == "3" ]]; then
    echo ""
    echo -e "${BLUE}======== Census Files ========${NC}"
    
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
    echo -e "${BLUE}Resume Mode: Using parameters from previous run${NC}"
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
    echo -e "${GREEN}Starting analysis...${NC}"
    
    # Log file for background jobs
    log_file="logs/foodnet_run_${timestamp}.log"
    
    # Execute the command
    if ! execute_command "$final_cmd" "$outDir" "$background" "$log_file"; then
      echo -e "${RED}Error in execution.${NC}"
      exit 1
    fi
  fi
}

# Run the main function
main "$@" 