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

# Initialize log file for error tracking
error_log="foodnet_errors.log"
echo "$(date): Starting FoodNet Trends Analysis Pipeline" > "$error_log"

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

# Setup modules based on environment
if command -v module &> /dev/null; then
    module purge
    module load nextflow/24.10.4
    module load singularity/4.1.4
    module load java/17.0.6
else
    echo "Warning: Module system not detected. Dependencies must be in PATH."
    echo "$(date): Module system not detected" >> "$error_log"
fi

# Set up TMPDIR if not already defined
TMPDIR=${TMPDIR:-/scicomp/scratch/$(whoami)}
mkdir -p "$TMPDIR/nextflow" 2>/dev/null

# Create timestamp for ID
timestamp=$(date +%Y%m%d_%H%M%S)

# Display welcome banner
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

echo "Mode selected: $workflow_mode"

# Default pathogen values
ALL_PATHOGENS="CAMPYLOBACTER,CYCLOSPORA,SALMONELLA,SHIGELLA,STEC,VIBRIO,YERSINIA"
DEFAULT_PATHOGENS="CAMPYLOBACTER,CYCLOSPORA"
ALL_STATES="CA,CO,CT,GA,MD,MN,NM,NY,OR,TN"

# Initialize variables
preprocessed_data=""
preprocessed_metadata=""
mmwrFile=""
censusFileB=""
censusFileP=""

# Mode 1: Preprocessing
if [[ "$workflow_mode" == "1" ]]; then
    echo ""
    echo "======== Input Files ========"
    
    # Default data file
    defaultMmwrFile="${DEFAULT_DATA_DIR}/mmwr9624_May2025.sas7bdat"
    
    read -p "MMWR data file [${defaultMmwrFile}]: " mmwrFile
    mmwrFile=${mmwrFile:-$defaultMmwrFile}
    
    # Validate file exists
    if [ ! -f "${mmwrFile}" ]; then
        echo "Error: MMWR file does not exist: ${mmwrFile}"
        echo "$(date): Missing MMWR file: ${mmwrFile}" >> "$error_log"
        echo "Exiting."
        exit 1
    fi
    
    # Set output location
    echo ""
    echo "======== Output Settings ========"
    
    # Default output directory for preprocessed data
    defaultPreprocessedDir="preprocessed_$(date +%Y%m%d_%H%M%S)"
    read -p "Output directory for preprocessed data [${defaultPreprocessedDir}]: " preprocessedDir
    preprocessedDir=${preprocessedDir:-$defaultPreprocessedDir}
    
    # Default output base name (derived from input filename)
    fileBasename=$(basename "${mmwrFile}" | sed 's/\.[^.]*$//')
    defaultOutputBase="foodnet_data_${fileBasename}"
    read -p "Base name for output files [${defaultOutputBase}]: " outputBase
    outputBase=${outputBase:-$defaultOutputBase}
    
    # Ask about metadata generation
    echo ""
    read -p "Generate metadata JSON? (y/n) [y]: " generate_metadata
    generate_metadata=${generate_metadata:-y}
    if [[ "$generate_metadata" =~ ^[Yy]$ ]]; then
        metadata_param="--generateMetadata true"
    else
        metadata_param="--generateMetadata false"
    fi
    
    echo "Running preprocessing..."
    
    # Run the preprocessing workflow
    preprocess_cmd="nextflow run main.nf -profile singularity -entry PREPROCESS_WORKFLOW \
      --mmwrFile \"${mmwrFile}\" \
      --outdir \"${preprocessedDir}\" \
      --outputBase \"${outputBase}\" \
      ${metadata_param}"
    
    echo "$(date): Running preprocessing command: ${preprocess_cmd}" >> "$error_log"
    
    if ! eval $preprocess_cmd; then
        echo "Preprocessing failed."
        echo "Check .nextflow.log for details."
        echo "$(date): Preprocessing failed, check .nextflow.log" >> "$error_log"
        exit 1
    fi
    
    # Set paths to preprocessed data and metadata
    preprocessed_data="${preprocessedDir}/preprocessed/${outputBase}.csv"
    preprocessed_metadata="${preprocessedDir}/preprocessed/${outputBase}_metadata.json"
    
    # Check for metadata in alternate location (for backward compatibility)
    if [ ! -f "${preprocessed_metadata}" ] && [ -f "${preprocessedDir}/preprocessed/metadata/${outputBase}_metadata.json" ]; then
        echo "Found metadata file in alternate location, using it instead."
        preprocessed_metadata="${preprocessedDir}/preprocessed/metadata/${outputBase}_metadata.json"
    fi
    
    # Check if files were created
    if [ ! -f "${preprocessed_data}" ]; then
        echo "Preprocessing failed to create expected CSV file: ${preprocessed_data}"
        echo "$(date): Missing expected output CSV: ${preprocessed_data}" >> "$error_log"
        echo "Check logs for details."
        exit 1
    fi
    
    echo "Preprocessing complete!"
    echo "- Cleaned data: ${preprocessed_data}"
    
    # Continue to analysis?
    echo ""
    read -p "Continue to analysis? (y/n) [y]: " continue_analysis
    continue_analysis=${continue_analysis:-y}
    if [[ ! "$continue_analysis" =~ ^[Yy]$ ]]; then
        echo "Preprocessing complete. Exiting."
        exit 0
    fi
    
    # Set the MMWR file for analysis to the preprocessed data
    mmwrFile=$preprocessed_data
    
# Mode 2: Analysis with raw data
elif [[ "$workflow_mode" == "2" ]]; then
    echo ""
    echo "======== Input Files ========"
    
    # Default data file
    defaultMmwrFile="${DEFAULT_DATA_DIR}/mmwr9624_May2025.sas7bdat"
    
    read -p "MMWR data file [${defaultMmwrFile}]: " mmwrFile
    mmwrFile=${mmwrFile:-$defaultMmwrFile}
    
    # Validate file exists
    if [ ! -f "${mmwrFile}" ]; then
        echo "Error: MMWR file does not exist: ${mmwrFile}"
        echo "$(date): Missing MMWR file: ${mmwrFile}" >> "$error_log"
        echo "Exiting."
        exit 1
    fi

# Mode 3: Use existing preprocessed data
else
    echo ""
    echo "======== Input Files ========"

    # Search for preprocessed CSV files
    echo "Searching for preprocessed CSV files..."
    mapfile -t found_csv < <(find . -type f -path "*/preprocessed/*.csv" 2>/dev/null)
    if [[ ${#found_csv[@]} -gt 0 ]]; then
        echo "Found the following preprocessed CSV files:"
        for i in "${!found_csv[@]}"; do
            printf "%2d) %s\n" $((i+1)) "${found_csv[$i]}"
        done
        echo "$(( ${#found_csv[@]} + 1 ))) Enter a file path manually"
        read -p "Select a CSV file [1]: " csv_choice
        csv_choice=${csv_choice:-1}
        if [[ "$csv_choice" -ge 1 && "$csv_choice" -le ${#found_csv[@]} ]]; then
            preprocessed_data="${found_csv[$((csv_choice-1))]}"
        else
            read -p "Preprocessed data file: " preprocessed_data
        fi
    else
        echo "No preprocessed CSV files found. Please enter the path manually."
        read -p "Preprocessed data file: " preprocessed_data
    fi

    # Validate file exists
    if [ ! -f "${preprocessed_data}" ]; then
        echo "Error: Preprocessed data file does not exist: ${preprocessed_data}"
        echo "$(date): Missing preprocessed data file: ${preprocessed_data}" >> "$error_log"
        echo "Exiting."
        exit 1
    fi

    # Set MMWR file to preprocessed data
    mmwrFile=$preprocessed_data

    # Search for preprocessed JSON metadata files
    echo ""
    echo "Searching for preprocessed metadata JSON files..."
    mapfile -t found_json < <(find . -type f -path "*/preprocessed/*.json" 2>/dev/null)
    if [[ ${#found_json[@]} -gt 0 ]]; then
        echo "Found the following preprocessed metadata JSON files:"
        for i in "${!found_json[@]}"; do
            printf "%2d) %s\n" $((i+1)) "${found_json[$i]}"
        done
        echo "$(( ${#found_json[@]} + 1 ))) Enter a file path manually or leave blank for none"
        read -p "Select a metadata file [${#found_json[@]}+1]: " json_choice
        json_choice=${json_choice:-$(( ${#found_json[@]} + 1 ))}
        if [[ "$json_choice" -ge 1 && "$json_choice" -le ${#found_json[@]} ]]; then
            preprocessed_metadata="${found_json[$((json_choice-1))]}"
        else
            read -p "Enter path to metadata file (leave blank if none): " preprocessed_metadata
        fi
    else
        echo "No preprocessed metadata JSON files found. Enter path manually or leave blank for none."
        read -p "Enter path to metadata file (leave blank if none): " preprocessed_metadata
    fi

    # Validate metadata file exists if provided
    if [ -n "${preprocessed_metadata}" ] && [ ! -f "${preprocessed_metadata}" ]; then
        # Try to find metadata in alternate location (in metadata subdirectory)
        metadata_dir=$(dirname "${preprocessed_metadata}")
        metadata_basename=$(basename "${preprocessed_metadata}")
        alt_metadata_path="${metadata_dir}/metadata/${metadata_basename}"

        if [ -f "${alt_metadata_path}" ]; then
            echo "Found metadata file in alternate location: ${alt_metadata_path}"
            preprocessed_metadata="${alt_metadata_path}"
        else
            echo "Warning: Metadata file does not exist: ${preprocessed_metadata}"
            echo "$(date): Missing metadata file: ${preprocessed_metadata}" >> "$error_log"
            read -p "Continue anyway? (y/n) [n]: " continue_choice
            continue_choice=${continue_choice:-n}
            if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
                echo "Exiting."
                exit 1
            fi
            # Clear metadata file if continuing without it
            preprocessed_metadata=""
        fi
    fi
fi

# For all modes, we need census files
if [[ "$workflow_mode" == "1" || "$workflow_mode" == "2" || "$workflow_mode" == "3" ]]; then
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
    if [ ! -f "${censusFileB}" ]; then
        echo "Error: Census bacterial file does not exist: ${censusFileB}"
        echo "$(date): Missing census bacterial file: ${censusFileB}" >> "$error_log"
        echo "Exiting."
        exit 1
    fi
    
    if [ ! -f "${censusFileP}" ]; then
        echo "Error: Census parasitic file does not exist: ${censusFileP}"
        echo "$(date): Missing census parasitic file: ${censusFileP}" >> "$error_log"
        echo "Exiting."
        exit 1
    fi
fi

# Load metadata if available
has_metadata=false
has_serotypes=false

if [[ -n "${preprocessed_metadata}" && -f "${preprocessed_metadata}" ]]; then
    echo "Using metadata file: ${preprocessed_metadata}"
    has_metadata=true
    
    # Check if file contains serotypes
    if grep -q "salmonella_serotypes" "${preprocessed_metadata}"; then
        has_serotypes=true
        echo "Metadata contains Salmonella serotype information."
    fi
    
    # Extract pathogens from metadata if available
    if [[ "$have_jq" == true ]]; then
        metadata_pathogens=$(jq -r '.pathogens | join(",")' "${preprocessed_metadata}" 2>/dev/null)
        if [ -n "$metadata_pathogens" ]; then
            echo "Pathogens in dataset: $metadata_pathogens"
            ALL_PATHOGENS="$metadata_pathogens"
        fi
        
        metadata_states=$(jq -r '.states | join(",")' "${preprocessed_metadata}" 2>/dev/null)
        if [ -n "$metadata_states" ]; then
            echo "States in dataset: $metadata_states"
            ALL_STATES="$metadata_states"
        fi
    fi
fi

# Add CAMPYLOBACTER,CYCLOSPORA to ALL_PATHOGENS if not present
if [[ ! "$ALL_PATHOGENS" == *"CAMPYLOBACTER"* ]]; then
    ALL_PATHOGENS="$ALL_PATHOGENS,CAMPYLOBACTER"
fi
if [[ ! "$ALL_PATHOGENS" == *"CYCLOSPORA"* ]]; then
    ALL_PATHOGENS="$ALL_PATHOGENS,CYCLOSPORA"
fi

# Get travel status filter
echo ""
echo "======== Travel Status Filter ========"
echo "Select travel status to include:"
echo "1) All travel statuses (NO, UNKNOWN, YES)"
echo "2) Only non-travel related (NO only)"
echo "3) Non-travel and unknown (NO, UNKNOWN)"
echo "4) Travel-related only (YES only)"
echo "5) Custom selection"
read -p "Enter selection [1]: " travel_choice
travel_choice=${travel_choice:-1}

case $travel_choice in
    1) travel="NO,UNKNOWN,YES" ;;
    2) travel="NO" ;;
    3) travel="NO,UNKNOWN" ;;
    4) travel="YES" ;;
    5)
        echo "Enter comma-separated travel statuses (NO,UNKNOWN,YES):"
        read -p "Travel statuses: " travel
        # Default if empty
        travel=${travel:-NO,UNKNOWN,YES}
        ;;
    *)
        echo "Invalid selection. Using default (All travel statuses)."
        travel="NO,UNKNOWN,YES"
        ;;
esac

# Get CIDT/culture method filter
echo ""
echo "======== CIDT/Culture Method Filter ========"
echo "Select CIDT/culture methods to include:"
echo "1) All methods (CIDT+, CX+, PARASITIC)"
echo "2) Culture positive only (CX+)"
echo "3) CIDT positive only (CIDT+)"
echo "4) Custom selection"
read -p "Enter selection [1]: " cidt_choice
cidt_choice=${cidt_choice:-1}

case $cidt_choice in
    1) cidt="CIDT+,CX+,PARASITIC" ;;
    2) cidt="CX+" ;;
    3) cidt="CIDT+" ;;
    4)
        echo "Enter comma-separated CIDT/culture methods (CIDT+,CX+,PARASITIC):"
        read -p "CIDT/culture methods: " cidt
        # Default if empty
        cidt=${cidt:-CIDT+,CX+,PARASITIC}
        ;;
    *)
        echo "Invalid selection. Using default (All methods)."
        cidt="CIDT+,CX+,PARASITIC"
        ;;
esac

# Get pathogen selection
echo ""
echo "======== Pathogen Selection ========"
echo "Available pathogens in this dataset:"

# Display each pathogen
IFS=',' read -ra PATHOGEN_ARRAY <<< "$ALL_PATHOGENS"
for p in "${PATHOGEN_ARRAY[@]}"; do
    # Skip empty items or just "metadata"
    if [[ -n "$p" && "$p" != "metadata" ]]; then
        echo "- $p"
    fi
done
echo ""

echo "1) Run ALL available pathogens"
echo "2) Select specific pathogens"
read -p "Enter selection [2]: " pathogen_mode
pathogen_mode=${pathogen_mode:-2}

if [[ "$pathogen_mode" == "1" ]]; then
    # Use all pathogens
    echo "Selected: ALL pathogens"
    pathogens="$ALL_PATHOGENS"
else
    # Ask for specific pathogens
    echo ""
    echo "Enter pathogens to analyze (comma-separated with NO spaces)"
    read -p "Leave blank for default (${DEFAULT_PATHOGENS}): " pathogens
    pathogens=${pathogens:-"$DEFAULT_PATHOGENS"}
fi

# After pathogen selection, add STEC serotype selection

# Only prompt for STEC serotypes if STEC is among selected pathogens
if [[ "$pathogens" == *"STEC"* ]]; then
    echo ""
    echo "======== STEC Serotype Selection ========"
    stec_serotype_list=""
    # Try to get STEC serotypes from metadata if available and jq is present
    if [[ "$have_jq" == true && -n "$preprocessed_metadata" && -f "$preprocessed_metadata" ]]; then
        stec_serotype_list=$(jq -r '.stec_serotypes | join(",")' "$preprocessed_metadata" 2>/dev/null)
    fi
    # If not found in metadata, scan the preprocessed CSV for unique STEC serotypes
    if [[ -z "$stec_serotype_list" && -n "$mmwrFile" && -f "$mmwrFile" ]]; then
        # Try to find the serotype column (serotypesummary, sero2, or sero1)
        serotype_col=$(head -1 "$mmwrFile" | tr ',' '\n' | grep -i -m1 -E 'serotypesummary|sero2|sero1')
        if [[ -n "$serotype_col" ]]; then
            stec_serotype_list=$(awk -F',' -v col="$serotype_col" 'NR==1{for(i=1;i<=NF;i++)if(tolower($i)==tolower(col))c=i} NR>1 && toupper($1)=="STEC" && c{a[$c]++} END{for(k in a) printf "%s,", k}' "$mmwrFile" | sed 's/,
*$//')
        fi
    fi
    if [[ -n "$stec_serotype_list" ]]; then
        IFS=',' read -ra STEC_SEROTYPES_ARRAY <<< "$stec_serotype_list"
        echo "Detected STEC serotypes in data:"
        for i in "${!STEC_SEROTYPES_ARRAY[@]}"; do
            printf "%2d) %s\n" $((i+1)) "${STEC_SEROTYPES_ARRAY[$i]}"
        done
        echo "$(( ${#STEC_SEROTYPES_ARRAY[@]} + 1 ))) Enter a custom list manually"
        read -p "Select STEC serotypes (comma-separated indices, or leave blank for all): " stec_sero_choice
        if [[ -z "$stec_sero_choice" ]]; then
            stec_serotypes=""
        elif [[ "$stec_sero_choice" -eq $(( ${#STEC_SEROTYPES_ARRAY[@]} + 1 )) ]]; then
            read -p "Enter STEC serotypes (comma-separated): " stec_serotypes
        else
            # Convert indices to serotype names
            stec_serotypes=""
            IFS=',' read -ra IDX <<< "$stec_sero_choice"
            for idx in "${IDX[@]}"; do
                idx=$((idx-1))
                if [[ $idx -ge 0 && $idx -lt ${#STEC_SEROTYPES_ARRAY[@]} ]]; then
                    stec_serotypes+="${STEC_SEROTYPES_ARRAY[$idx]},"
                fi
            done
            stec_serotypes=$(echo "$stec_serotypes" | sed 's/,
*$//')
        fi
    else
        echo "No STEC serotype list detected. You may enter a custom list or leave blank for all."
        read -p "STEC serotypes to include: " stec_serotypes
    fi
else
    stec_serotypes=""
fi

# Get state selection
echo ""
echo "======== State Selection ========"
echo "Available states in this dataset:"

# Display each state
IFS=',' read -ra STATE_ARRAY <<< "$ALL_STATES"
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
    echo "Selected: ALL states"
    states="$ALL_STATES"
else
    # Ask for specific states
    echo ""
    echo "Enter states to analyze (comma-separated with NO spaces)"
    read -p "Leave blank for all states: " states
    states=${states:-"$ALL_STATES"}
fi

# Get run mode
echo ""
echo "Select run mode:"
echo "1) Test run (minimal settings - automatically sets: chains=1, iterations=100)"
echo "2) Full analysis (custom settings - you'll specify parameters)"
echo "3) Resume previous run (continues from last execution)"
read -p "Enter selection [1]: " run_mode
run_mode=${run_mode:-1}

# Convert run mode to flag for backward compatibility
case $run_mode in
    1) 
        flag="test"
        chains=1
        iterations=100
        adapt_delta=0.8
        max_treedepth=8
        echo ""
        echo "Test Mode: Using minimal settings"
        echo "Chains: $chains"
        echo "Iterations: $iterations"
        echo "Adapt delta: $adapt_delta"
        echo "Max treedepth: $max_treedepth"
        ;;
    2)
        flag="full"
        # Get MCMC parameters
        echo ""
        read -p "Number of chains [2]: " chains
        chains=${chains:-2}
        
        echo ""
        read -p "Number of iterations [500]: " iterations
        iterations=${iterations:-500}
        
        echo ""
        read -p "Adapt delta (0.0-1.0) [0.95]: " adapt_delta
        adapt_delta=${adapt_delta:-0.95}
        
        echo ""
        read -p "Max treedepth [10]: " max_treedepth
        max_treedepth=${max_treedepth:-10}
        ;;
    3)
        flag="resume"
        # Use default values for resumed runs
        chains=2
        iterations=500
        adapt_delta=0.95
        max_treedepth=10
        echo ""
        echo "Resume Mode: Using parameters from previous run"
        ;;
    *)
        echo "Invalid selection. Using default (Test run)."
        flag="test"
        chains=1
        iterations=100
        adapt_delta=0.8
        max_treedepth=8
        ;;
esac

# Get background execution preference
echo ""
read -p "Run in background? (y/n) [n]: " bg_choice
bg_choice=${bg_choice:-n}
if [[ "$bg_choice" =~ ^[Yy]$ ]]; then
    background=true
else
    background=false
fi

# Get output directory
echo ""
read -p "Output directory [${outDir}]: " user_outdir
outDir=${user_outdir:-$outDir}
# Trim any leading/trailing spaces
outDir=$(echo "$outDir" | xargs)

# Dashboard preferences - simplified and automatic
echo ""
echo "======== Dashboard Generation ========"
echo "An interactive HTML dashboard will be automatically generated with run details included."

# Automatically generate dashboard title with timestamp for identification
dashboard_title="FoodNet Trends Analysis - Run ${timestamp}"
dashboard_params="--enable_dashboard true --dashboard_title \"$dashboard_title\""

# Tell the user what's happening
echo "Dashboard will be created with title: \"$dashboard_title\""

# Build the command
cmd="nextflow run main.nf -profile singularity"
if [[ "$flag" == "resume" ]]; then
    cmd="$cmd -resume"
fi

cmd="$cmd --censusFileB \"${censusFileB}\""
cmd="$cmd --censusFileP \"${censusFileP}\""
cmd="$cmd --travel \"${travel}\""
cmd="$cmd --cidt \"${cidt}\""
cmd="$cmd --iterations ${iterations}"
cmd="$cmd --chains ${chains}"
cmd="$cmd --adapt_delta ${adapt_delta}"
cmd="$cmd --max_treedepth ${max_treedepth}"
cmd="$cmd --seed 123"
cmd="$cmd --outdir \"${outDir}\""
cmd="$cmd --pathogen \"${pathogens}\""
cmd="$cmd --states \"${states}\""

# Add MMWR file parameter based on workflow mode
if [[ "$workflow_mode" == "3" ]]; then
    # Using preprocessed data
    cmd="$cmd --mmwrFile \"$mmwrFile\" --preprocessed true --cleanFile \"$mmwrFile\""
    
    # Add metadata file if available
    if [[ -n "$preprocessed_metadata" ]]; then
        cmd="$cmd --metadata \"$preprocessed_metadata\""
    fi
else
    # Using raw data
    cmd="$cmd --mmwrFile \"$mmwrFile\""
    
    # If this is workflow_mode 1, we already preprocessed
    if [[ "$workflow_mode" == "1" ]]; then
        cmd="$cmd --preprocessed true --cleanFile \"$mmwrFile\""
        if [[ -n "$preprocessed_metadata" ]]; then
            cmd="$cmd --metadata \"$preprocessed_metadata\""
        fi
    fi
fi

# Add dashboard parameters
cmd="$cmd $dashboard_params"

# Add background option if needed
if [[ "$background" == true ]]; then
    # Create logs directory if it doesn't exist
    mkdir -p logs
    log_file="logs/foodnet_run_${timestamp}.log"
    cmd="nohup $cmd > \"${log_file}\" 2>&1 &"
    echo "Process will run in background with log: ${log_file}"
fi

# When building the Nextflow command, add:
if [[ -n "$stec_serotypes" ]]; then
    cmd="$cmd --stec_serotypes \"$stec_serotypes\""
fi

# Display summary
echo ""
echo "========= Analysis Summary ==========="
if [[ "$flag" == "test" ]]; then
    echo "Mode: Test run"
elif [[ "$flag" == "full" ]]; then
    echo "Mode: Full analysis"
else
    echo "Mode: Resume previous run"
fi

if [[ "$workflow_mode" == "1" || "$workflow_mode" == "3" ]]; then
    echo "Using preprocessed data: Yes"
else
    echo "Using preprocessed data: No"
fi

echo "Pathogens: $pathogens"
echo "States: $states"
echo "Travel status: $travel"
echo "CIDT/culture methods: $cidt"

echo "Input files:"
echo "  Data file: $mmwrFile"
echo "  Census file (bacterial): $censusFileB" 
echo "  Census file (parasitic): $censusFileP"
if [[ -n "$preprocessed_metadata" ]]; then
    echo "  Metadata: $preprocessed_metadata"
fi

if [[ "$flag" != "resume" ]]; then
    echo "Chains: $chains"
    echo "Iterations: $iterations"
    echo "Adapt delta: $adapt_delta"
    echo "Max treedepth: $max_treedepth"
fi

echo "Output directory: $outDir"
echo "Run in background: $([ "$background" == true ] && echo "Yes" || echo "No")"
echo ""
echo "Command to run:"
echo "$cmd"
echo ""

# Get confirmation from user
read -p "Execute command? (y/n) [y]: " execute
execute=${execute:-y}

if [[ "$execute" =~ ^[Yy]$ ]]; then
    echo "Starting analysis..."
    echo "$(date): Executing command: $cmd" >> "$error_log"
    
    if ! eval $cmd; then
        echo "Error running analysis command."
        echo "Check .nextflow.log for details."
        echo "$(date): Command execution failed" >> "$error_log"
        exit 1
    fi
else
    echo "Execution canceled."
    echo "$(date): User canceled execution" >> "$error_log"
fi 