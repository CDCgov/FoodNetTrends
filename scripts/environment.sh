#!/bin/bash

# Module for handling environment setup for FoodNet Trends pipeline
# This includes module loading, container checking, and path configuration
# Version: 1.0 (2025-05)

# Set up modules for HPC environment
setup_modules() {
  echo "DEBUG: Setting up modules"
  module purge
  module load nextflow/24.10.4
  module load singularity/4.1.4
  module load java/17.0.6
  echo "DEBUG: Modules loaded"
}

# Check if Singularity container exists
check_container() {
  if [ ! -f "foodnet.sif" ]; then
    echo "Warning: Singularity container (foodnet.sif) not found."
    echo "You may need to build the container first with:"
    echo "singularity build foodnet.sif foodnet.def"
    echo "$(date): Container not found: foodnet.sif" >> "$error_log"
    
    read -p "Continue anyway? This will likely fail. (y/n) [n]: " continue_choice
    continue_choice=${continue_choice:-n}
    if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
      echo "Exiting."
      exit 1
    fi
  fi
}

# Set up data directories based on environment
setup_data_paths() {
  echo "DEBUG: Setting up data paths"
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
  
  echo "DEBUG: Data directory: $DEFAULT_DATA_DIR"
  # Return the data directory for use in the main script
  echo "$DEFAULT_DATA_DIR"
}

# Setup environment variables
setup_environment() {
  echo "DEBUG: Starting setup_environment"

  # Initialize log file for error tracking
  timestamp=$(date +%Y%m%d_%H%M%S)
  error_log="foodnet_errors.log"
  echo "$(date): Starting FoodNet Trends Analysis Pipeline" > "$error_log"
  echo "DEBUG: Error log initialized at $error_log"
  
  # Make error_log globally available
  export error_log
  
  # Set up modules if we're in an HPC environment
  if command -v module &> /dev/null; then
    echo "DEBUG: Module command found, setting up modules"
    setup_modules
  else
    echo "DEBUG: Module command not found"
    echo "Warning: Module system not detected. Assuming dependencies are available in PATH."
    echo "$(date): Module system not detected" >> "$error_log"
  fi
  
  # Check for the Singularity container
  echo "DEBUG: Checking for Singularity container"
  check_container
  
  # Set up data paths
  echo "DEBUG: Setting up data paths"
  DEFAULT_DATA_DIR=$(setup_data_paths)
  export DEFAULT_DATA_DIR
  echo "DEBUG: DEFAULT_DATA_DIR set to $DEFAULT_DATA_DIR"
  
  # Create timestamp for automatic project ID
  timestamp=$(date +%Y%m%d_%H%M%S)
  export timestamp
  
  # Initialize output directory
  outDir="output"  # Default to "output" directory in current location
  export outDir
  
  echo "DEBUG: Finished setup_environment"
}

# Check if jq is available for JSON processing
check_jq() {
  echo "DEBUG: Checking for jq"
  have_jq=true
  if ! command -v jq &> /dev/null; then
    echo "Warning: jq is not installed. Basic functionality will work, but advanced serotype filtering will be limited."
    have_jq=false
  fi
  
  echo "DEBUG: have_jq = $have_jq"
  # Export have_jq so it's available to other scripts
  export have_jq
  
  echo "$have_jq"
} 