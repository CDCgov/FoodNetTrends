#!/bin/bash

# UI functions for FoodNet Trends pipeline
# This includes functions for user interaction and display
# Version: 1.0 (2025-05)

# Define color codes for better user experience
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to display the welcome banner
display_welcome() {
  # Make sure TMPDIR is set before referencing it
  TMPDIR=${TMPDIR:-/scicomp/scratch/$(whoami)}
  
  # Ensure nextflow directory exists
  mkdir -p "$TMPDIR/nextflow" 2>/dev/null
  
  echo "Your Nextflow temporary/cache files will be placed in ${TMPDIR}/nextflow/ by default"
  echo "========================================="
  echo "   FoodNet Trends Analysis Pipeline      "
  echo "========================================="
  echo ""
}

# Function to get workflow mode from the user
get_workflow_mode() {
  # Explicitly print each option with plain text formatting
  echo "Select mode:"
  echo ""
  echo "1) Preprocess data (clean raw data files and generate metadata)"
  echo "2) Run analysis (with complete pipeline)"
  echo "3) Use existing preprocessed data"
  echo ""
  # Add debug output
  echo "Waiting for your selection (enter 1, 2, or 3)..."
  read -p "Enter selection [1]: " workflow_mode
  workflow_mode=${workflow_mode:-1}
  
  # Echo the selection for debugging
  echo "You selected: $workflow_mode"
  
  # Validate workflow mode
  if [[ ! "$workflow_mode" =~ ^[1-3]$ ]]; then
    echo "Error: Invalid mode selection. Using default (Preprocess data)."
    echo "$(date): Invalid workflow_mode: ${workflow_mode}" >> "$error_log"
    workflow_mode=1
  fi
  
  echo "$workflow_mode"
}

# Function to get run mode (test/full/resume)
get_run_mode() {
  echo ""
  echo -e "Select run mode:"
  echo "1) Test run (minimal settings - automatically sets: chains=1, iterations=100)"
  echo "2) Full analysis (custom settings - you'll specify parameters)"
  echo "3) Resume previous run (continues from last execution)"
  read -p "Enter selection [1]: " run_mode
  run_mode=${run_mode:-1}
  
  # Validate run mode
  if [[ ! "$run_mode" =~ ^[1-3]$ ]]; then
    echo -e "${RED}Error: Invalid selection. Using default (Test run).${NC}"
    echo "$(date): Invalid run_mode selection: ${run_mode}" >> "$error_log"
    run_mode=1
  fi
  
  # Convert run mode to flag for backward compatibility
  case $run_mode in
    1) echo "test" ;;
    2) echo "full" ;;
    3) echo "resume" ;;
  esac
}

# Function to prompt for background execution
get_background_choice() {
  echo ""
  read -p "Run in background? (y/n) [n]: " bg_choice
  bg_choice=${bg_choice:-n}
  if [[ "$bg_choice" =~ ^[Yy]$ ]]; then
    echo "true"
  else
    echo "false"
  fi
}

# Function to get output directory
get_output_directory() {
  local default_dir=$1
  echo ""
  read -p "Output directory [${default_dir}]: " user_outdir
  echo "${user_outdir:-$default_dir}"
}

# Function to get travel status selection
get_travel_status() {
  echo ""
  echo -e "${BLUE}======== Travel Status Filter ========${NC}"
  echo "Select travel status to include:"
  echo "1) All travel statuses (NO, UNKNOWN, YES)"
  echo "2) Only non-travel related (NO only)"
  echo "3) Non-travel and unknown (NO, UNKNOWN)"
  echo "4) Travel-related only (YES only)"
  echo "5) Custom selection"
  read -p "Enter selection [1]: " travel_choice
  travel_choice=${travel_choice:-1}
  
  case $travel_choice in
    1) echo "NO,UNKNOWN,YES" ;;
    2) echo "NO" ;;
    3) echo "NO,UNKNOWN" ;;
    4) echo "YES" ;;
    5)
      echo "Enter comma-separated travel statuses (NO,UNKNOWN,YES):"
      read -p "Travel statuses: " travel
      # Default if empty
      echo "${travel:-NO,UNKNOWN,YES}"
      ;;
    *)
      echo -e "${RED}Error: Invalid selection. Using default (All travel statuses).${NC}"
      echo "$(date): Invalid travel_choice: ${travel_choice}" >> "$error_log"
      echo "NO,UNKNOWN,YES"
      ;;
  esac
}

# Function to get CIDT/culture method selection
get_cidt_method() {
  echo ""
  echo -e "${BLUE}======== CIDT/Culture Method Filter ========${NC}"
  echo "Select CIDT/culture methods to include:"
  echo "1) All methods (CIDT+, CX+, PARASITIC)"
  echo "2) Culture positive only (CX+)"
  echo "3) CIDT positive only (CIDT+)"
  echo "4) Custom selection"
  read -p "Enter selection [1]: " cidt_choice
  cidt_choice=${cidt_choice:-1}
  
  case $cidt_choice in
    1) echo "CIDT+,CX+,PARASITIC" ;;
    2) echo "CX+" ;;
    3) echo "CIDT+" ;;
    4)
      echo "Enter comma-separated CIDT/culture methods (CIDT+,CX+,PARASITIC):"
      read -p "CIDT/culture methods: " cidt
      # Default if empty
      echo "${cidt:-CIDT+,CX+,PARASITIC}"
      ;;
    *)
      echo -e "${RED}Error: Invalid selection. Using default (All methods).${NC}"
      echo "$(date): Invalid cidt_choice: ${cidt_choice}" >> "$error_log"
      echo "CIDT+,CX+,PARASITIC"
      ;;
  esac
}

# Function to display analysis summary
display_summary() {
  local mode_string=$1
  local workflow_mode=$2
  local pathogens=$3
  local states=$4
  local travel=$5
  local cidt=$6
  local serotype_param=$7
  local mmwrFile=$8
  local censusFileB=$9
  local censusFileP=${10}
  local preprocessed_metadata=${11}
  local flag=${12}
  local chains=${13}
  local iterations=${14}
  local adapt_delta=${15}
  local max_treedepth=${16}
  local outDir=${17}
  local background=${18}
  local final_cmd=${19}
  
  echo ""
  echo -e "${BLUE}========= Analysis Summary ===========${NC}"
  echo -e "Mode: ${GREEN}${mode_string}${NC}"
  
  # Fixed conditional for preprocessed data:
  if [ "$workflow_mode" == "1" ] || [ "$workflow_mode" == "3" ]; then
    echo -e "Using preprocessed data: ${GREEN}Yes${NC}"
  else
    echo -e "Using preprocessed data: ${GREEN}No${NC}"
  fi
  
  echo -e "Pathogens: ${GREEN}$pathogens${NC}"
  echo -e "States: ${GREEN}$states${NC}"
  echo -e "Travel status: ${GREEN}$travel${NC}"
  echo -e "CIDT/culture methods: ${GREEN}$cidt${NC}"
  if [[ -n "$serotype_param" ]]; then
    # Extract just the serotype names from the parameter
    serotypes=$(echo $serotype_param | sed 's/--salmonella_serotypes //; s/"//g')
    echo -e "Salmonella serotypes: ${GREEN}$serotypes${NC}"
  fi
  echo -e "Input files:"
  echo -e "  Data file: ${GREEN}$mmwrFile${NC}"
  echo -e "  Census file (bacterial): ${GREEN}$censusFileB${NC}" 
  echo -e "  Census file (parasitic): ${GREEN}$censusFileP${NC}"
  if [[ -n "$preprocessed_metadata" ]]; then
    echo -e "  Metadata: ${GREEN}$preprocessed_metadata${NC}"
  fi
  
  if [[ "$flag" != "resume" ]]; then
    echo -e "Chains: ${GREEN}$chains${NC}"
    echo -e "Iterations: ${GREEN}$iterations${NC}"
    echo -e "Adapt delta: ${GREEN}$adapt_delta${NC}"
    echo -e "Max treedepth: ${GREEN}$max_treedepth${NC}"
  fi
  
  echo -e "Output directory: ${GREEN}$outDir${NC}"
  echo -e "Run in background: ${GREEN}$([ "$background" == true ] && echo "Yes" || echo "No")${NC}"
  echo ""
  echo -e "Command to run:"
  echo -e "${YELLOW}$final_cmd${NC}"
  echo ""
} 