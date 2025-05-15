#!/bin/bash

# set up paths, files
DEFAULT_DATA_DIR="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/"
outDir="output"  # Default to "output" directory in current location

# set up modules
module purge
module load nextflow/24.10.4
module load singularity/4.1.4
module load java/17.0.6

# Create timestamp for automatic project ID
timestamp=$(date +%Y%m%d_%H%M%S)

# Define color codes for better user experience
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "Your Nextflow temporary/cache files will be placed in $TMPDIR/nextflow/ by default"
echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}   FoodNet Trends Analysis Pipeline      ${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

# Ask for workflow mode
echo -e "Select mode:"
echo "1) Preprocess data (clean raw data files and generate metadata)"
echo "2) Run analysis (with complete pipeline)"
echo "3) Use existing preprocessed data"
read -p "Enter selection [1]: " workflow_mode
workflow_mode=${workflow_mode:-1}

# Handle preprocessed data
preprocessed_data=""
preprocessed_metadata=""

# Mode 1: Preprocessing
if [[ "$workflow_mode" == "1" ]]; then
    echo ""
    echo -e "${BLUE}======== Input Files ========${NC}"
    
    # Default data file
    defaultMmwrFile="${DEFAULT_DATA_DIR}/mmwr9624_May2025.sas7bdat"
    
    read -p "MMWR data file [${defaultMmwrFile}]: " mmwrFile
    mmwrFile=${mmwrFile:-$defaultMmwrFile}
    
    # Validate file exists
    if [ ! -f "$mmwrFile" ]; then
        echo -e "${RED}Error: MMWR file does not exist: $mmwrFile${NC}"
        echo -e "${RED}Exiting.${NC}"
        exit 1
    fi
    
    # Set output location
    echo ""
    echo -e "${BLUE}======== Output Settings ========${NC}"
    
    # Default output directory for preprocessed data
    defaultPreprocessedDir="preprocessed_$(date +%Y%m%d_%H%M%S)"
    read -p "Output directory for preprocessed data [${defaultPreprocessedDir}]: " preprocessedDir
    preprocessedDir=${preprocessedDir:-$defaultPreprocessedDir}
    
    # Default output base name (derived from input filename)
    fileBasename=$(basename "$mmwrFile" | sed 's/\.[^.]*$//')
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
    
    echo -e "${GREEN}Running preprocessing...${NC}"
    
    # Run the preprocessing workflow
    preprocess_cmd="nextflow run main.nf -profile singularity -entry PREPROCESS_WORKFLOW \
      --mmwrFile \"$mmwrFile\" \
      --outdir \"$preprocessedDir\" \
      --outputBase \"$outputBase\" \
      $metadata_param"
    
    if ! eval $preprocess_cmd; then
        echo -e "${RED}Preprocessing failed.${NC}"
        echo -e "${RED}Check .nextflow.log for details.${NC}"
        exit 1
    fi
    
    # Set paths to preprocessed data and metadata
    preprocessed_data="$preprocessedDir/preprocessed/${outputBase}.csv"
    preprocessed_metadata="$preprocessedDir/preprocessed/${outputBase}_metadata.json"
    
    # Check if files were created
    if [ ! -f "$preprocessed_data" ]; then
        echo -e "${RED}Preprocessing failed to create expected CSV file.${NC}"
        echo -e "${RED}Check logs for details.${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}Preprocessing complete!${NC}"
    echo -e "- Cleaned data: ${GREEN}$preprocessed_data${NC}"
    
    if [[ "$generate_metadata" =~ ^[Yy]$ ]] && [ -f "$preprocessed_metadata" ]; then
        echo -e "- Metadata: ${GREEN}$preprocessed_metadata${NC}"
    fi
    
    echo ""
    read -p "Proceed to analysis with this preprocessed data? (y/n) [y]: " proceed_to_analysis
    proceed_to_analysis=${proceed_to_analysis:-y}
    
    if [[ ! "$proceed_to_analysis" =~ ^[Yy]$ ]]; then
        echo -e "${GREEN}Preprocessing complete. You can run analysis later using mode 3.${NC}"
        exit 0
    fi
    
    # Continue to analysis using the preprocessed data
    echo -e "${GREEN}Proceeding to analysis...${NC}"

# Mode 3: Use existing preprocessed data
elif [[ "$workflow_mode" == "3" ]]; then
    echo ""
    echo -e "${BLUE}======== Preprocessed Data ========${NC}"
    
    read -p "Path to preprocessed CSV file: " preprocessed_data
    
    # Validate path exists
    if [ ! -f "$preprocessed_data" ]; then
        echo -e "${RED}Error: Preprocessed file does not exist: $preprocessed_data${NC}"
        echo -e "${RED}Exiting.${NC}"
        exit 1
    fi
    
    # Check for metadata file in same directory with _metadata.json suffix
    base_path=${preprocessed_data%.csv}
    auto_metadata="${base_path}_metadata.json"
    
    if [ -f "$auto_metadata" ]; then
        echo -e "${GREEN}Found metadata: $auto_metadata${NC}"
        preprocessed_metadata="$auto_metadata"
    else
        echo -e "${YELLOW}No metadata file found with naming pattern ${base_path}_metadata.json${NC}"
        read -p "Path to metadata JSON file (leave empty to skip): " user_metadata
        
        if [ -n "$user_metadata" ]; then
            if [ -f "$user_metadata" ]; then
                preprocessed_metadata="$user_metadata"
                echo -e "${GREEN}Using metadata: $preprocessed_metadata${NC}"
            else
                echo -e "${RED}Metadata file does not exist: $user_metadata${NC}"
                echo -e "${YELLOW}Continuing without metadata. Discovery will be limited.${NC}"
            fi
        else
            echo -e "${YELLOW}Continuing without metadata. Discovery will be limited.${NC}"
        fi
    fi
fi

# Load metadata if available
ALL_PATHOGENS="CAMPYLOBACTER,CYCLOSPORA,SALMONELLA,SHIGELLA,STEC,VIBRIO,YERSINIA"
ALL_STATES="CA,CO,CT,GA,MD,MN,NM,NY,OR,TN"
DEFAULT_PATHOGENS="CAMPYLOBACTER,CYCLOSPORA"
has_serotypes=false

if [[ -n "$preprocessed_metadata" && -f "$preprocessed_metadata" ]]; then
    echo -e "${BLUE}Reading metadata from: $preprocessed_metadata${NC}"
    
    if command -v jq &> /dev/null; then
        echo -e "${BLUE}======== Available Data ========${NC}"
        
        # Try to read pathogens with error handling
        if ! pathogens_list=$(jq -r '.pathogens[]' "$preprocessed_metadata" 2>/dev/null); then
            echo -e "${YELLOW}Error reading pathogens from metadata.${NC}"
        else
            echo -e "${BLUE}Pathogens in dataset:${NC}"
            echo "$pathogens_list" | sort | sed 's/^/- /'
            
            # Count pathogens
            if ! pathogen_count=$(jq -r '.pathogens | length' "$preprocessed_metadata" 2>/dev/null); then
                pathogen_count="Unknown"
            fi
            echo -e "${GREEN}Total: $pathogen_count pathogens${NC}"
            
            # Create comma-separated list
            ALL_PATHOGENS=$(echo "$pathogens_list" | tr '\n' ',' | sed 's/,$//')
        fi
        
        echo ""
        # Try to read states with error handling
        if ! states_list=$(jq -r '.states[]' "$preprocessed_metadata" 2>/dev/null); then
            echo -e "${YELLOW}Error reading states from metadata.${NC}"
        else
            echo -e "${BLUE}States in dataset:${NC}"
            echo "$states_list" | sort | sed 's/^/- /'
            
            # Count states
            if ! state_count=$(jq -r '.states | length' "$preprocessed_metadata" 2>/dev/null); then
                state_count="Unknown"
            fi
            echo -e "${GREEN}Total: $state_count states${NC}"
            
            # Create comma-separated list
            ALL_STATES=$(echo "$states_list" | tr '\n' ',' | sed 's/,$//')
        fi
        
        # Set intelligent defaults based on discovery
        # Pick the first 2 pathogens instead of hardcoding
        if [[ -n "$ALL_PATHOGENS" ]]; then
            DEFAULT_PATHOGENS=$(echo "$ALL_PATHOGENS" | cut -d',' -f1,2)
        fi
        
        # Check for Salmonella serotypes with robust error handling
        if jq -e '.salmonella_serotypes' "$preprocessed_metadata" > /dev/null 2>&1; then
            has_serotypes=true
            echo ""
            echo -e "${BLUE}Top Salmonella serotypes in dataset:${NC}"
            jq_cmd='.salmonella_serotypes | to_entries | sort_by(.value) | reverse | .[0:10] | .[] | "\(.key): \(.value) isolates"'
            if ! top_serotypes=$(jq -r "$jq_cmd" "$preprocessed_metadata" 2>/dev/null); then
                echo -e "${YELLOW}Could not process serotype information with jq. Displaying raw counts instead.${NC}"
                jq -r '.salmonella_serotypes | keys | .[0:10]' "$preprocessed_metadata" 2>/dev/null | sed 's/^/- /'
            else
                echo "$top_serotypes" | sed 's/^/- /'
            fi
        fi
    else
        echo -e "${YELLOW}jq not installed. Cannot parse JSON metadata.${NC}"
        echo -e "${YELLOW}Continuing with default values.${NC}"
    fi
fi

# For mode 2 (without preprocessing first), ask for input files
if [[ "$workflow_mode" == "2" ]]; then
    echo ""
    echo -e "${BLUE}======== Input Files ========${NC}"
    
    # Default data files
    defaultMmwrFile="${DEFAULT_DATA_DIR}/mmwr9624_May2025.sas7bdat"
    defaultCensusFileB="${DEFAULT_DATA_DIR}/cen9624.sas7bdat"
    defaultCensusFileP="${DEFAULT_DATA_DIR}/cen9624_para.sas7bdat"
    
    read -p "MMWR data file [${defaultMmwrFile}]: " mmwrFile
    mmwrFile=${mmwrFile:-$defaultMmwrFile}
    
    read -p "Census file (bacterial) [${defaultCensusFileB}]: " censusFileB
    censusFileB=${censusFileB:-$defaultCensusFileB}
    
    read -p "Census file (parasitic) [${defaultCensusFileP}]: " censusFileP
    censusFileP=${censusFileP:-$defaultCensusFileP}
    
    # Validate that files exist
    for file in "$mmwrFile" "$censusFileB" "$censusFileP"; do
        if [ ! -f "$file" ]; then
            echo -e "${RED}Warning: File does not exist: $file${NC}"
            read -p "Continue anyway? (y/n) [n]: " continue_choice
            continue_choice=${continue_choice:-n}
            if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
                echo -e "${RED}Exiting.${NC}"
                exit 1
            fi
        fi
    done
else
    # For modes 1 and 3, we need to ask for census files since they're not preprocessed
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
    for file in "$censusFileB" "$censusFileP"; do
        if [ ! -f "$file" ]; then
            echo -e "${RED}Warning: Census file does not exist: $file${NC}"
            read -p "Continue anyway? (y/n) [n]: " continue_choice
            continue_choice=${continue_choice:-n}
            if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
                echo -e "${RED}Exiting.${NC}"
                exit 1
            fi
        fi
    done
    
    # For preprocessed modes, set mmwrFile to the preprocessed CSV
    mmwrFile=$preprocessed_data
fi

# === Analysis Setup (All Modes) ===

# Ask for run mode
echo ""
echo -e "Select run mode:"
echo "1) Test run (minimal settings - automatically sets: chains=1, iterations=100)"
echo "2) Full analysis (custom settings - you'll specify parameters)"
echo "3) Resume previous run (continues from last execution)"
read -p "Enter selection [1]: " run_mode
run_mode=${run_mode:-1}

# Validate run mode
if [[ ! "$run_mode" =~ ^[1-3]$ ]]; then
    echo -e "${RED}Invalid selection. Using default (Test run).${NC}"
    run_mode=1
fi

# Convert run mode to flag for backward compatibility
case $run_mode in
    1) flag="test" ;;
    2) flag="full" ;;
    3) flag="resume" ;;
esac

# Ask if user wants to run in background
echo ""
read -p "Run in background? (y/n) [n]: " bg_choice
bg_choice=${bg_choice:-n}
if [[ "$bg_choice" =~ ^[Yy]$ ]]; then
    background=true
else
    background=false
fi

# Ask for output directory
echo ""
read -p "Output directory [${outDir}]: " user_outdir
outDir=${user_outdir:-$outDir}

# Ask about travel status
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
    1) travel="NO,UNKNOWN,YES" ;;
    2) travel="NO" ;;
    3) travel="NO,UNKNOWN" ;;
    4) travel="YES" ;;
    5)
        echo "Enter comma-separated travel statuses (NO,UNKNOWN,YES):"
        read -p "Travel statuses: " travel
        # Default if empty
        travel=${travel:-"NO,UNKNOWN,YES"}
        ;;
    *)
        echo -e "${RED}Invalid selection. Using default (All travel statuses).${NC}"
        travel="NO,UNKNOWN,YES"
        ;;
esac

# Ask about CIDT/culture method
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
    1) cidt="CIDT+,CX+,PARASITIC" ;;
    2) cidt="CX+" ;;
    3) cidt="CIDT+" ;;
    4)
        echo "Enter comma-separated CIDT/culture methods (CIDT+,CX+,PARASITIC):"
        read -p "CIDT/culture methods: " cidt
        # Default if empty
        cidt=${cidt:-"CIDT+,CX+,PARASITIC"}
        ;;
    *)
        echo -e "${RED}Invalid selection. Using default (All methods).${NC}"
        cidt="CIDT+,CX+,PARASITIC"
        ;;
esac

# Ask about pathogens
echo ""
echo -e "${BLUE}======== Pathogen Selection ========${NC}"
echo "Available pathogens in this dataset:"
# Parse the comma-separated list and display each pathogen
IFS=',' read -ra PATHOGEN_ARRAY <<< "$ALL_PATHOGENS"
for p in "${PATHOGEN_ARRAY[@]}"; do
    echo "- $p"
done
echo ""

echo "1) Run ALL available pathogens"
echo "2) Select specific pathogens"
read -p "Enter selection [2]: " pathogen_mode
pathogen_mode=${pathogen_mode:-2}

if [[ "$pathogen_mode" == "1" ]]; then
    # Use all pathogens
    pathogens=$ALL_PATHOGENS
    echo -e "${GREEN}Selected: ALL pathogens (${ALL_PATHOGENS})${NC}"
else
    # Ask for specific pathogens
    echo ""
    echo -e "Enter pathogens to analyze (comma-separated with NO spaces)"
    read -p "Leave blank for default (${DEFAULT_PATHOGENS}): " pathogens
    pathogens=${pathogens:-"$DEFAULT_PATHOGENS"}

    # Validate pathogens against the discovered list
    IFS=',' read -ra PATHOGEN_ARRAY <<< "$pathogens"
    IFS=',' read -ra VALID_PATHOGENS <<< "$ALL_PATHOGENS"
    invalid_found=false

    for p in "${PATHOGEN_ARRAY[@]}"; do
        valid=false
        for vp in "${VALID_PATHOGENS[@]}"; do
            if [[ "$p" == "$vp" ]]; then
                valid=true
                break
            fi
        done
        
        if [[ "$valid" == false ]]; then
            echo -e "${YELLOW}Warning: '$p' is not in the discovered pathogen list and may cause errors.${NC}"
            invalid_found=true
        fi
    done

    if [[ "$invalid_found" == true ]]; then
        echo ""
        read -p "Continue anyway? (y/n) [n]: " continue_choice
        continue_choice=${continue_choice:-n}
        if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
            echo -e "${RED}Exiting.${NC}"
            exit 1
        fi
    fi
fi

# Add state selection
echo ""
echo -e "${BLUE}======== State Selection ========${NC}"
echo "Available states in this dataset:"
# Parse the comma-separated list and display each state
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
    states=$ALL_STATES
    echo -e "${GREEN}Selected: ALL states (${ALL_STATES})${NC}"
else
    # Ask for specific states
    echo ""
    echo -e "Enter states to analyze (comma-separated with NO spaces)"
    read -p "Leave blank for all states: " states
    states=${states:-"$ALL_STATES"}

    # Validate states against the discovered list
    IFS=',' read -ra STATE_ARRAY <<< "$states"
    IFS=',' read -ra VALID_STATES <<< "$ALL_STATES"
    invalid_found=false

    for s in "${STATE_ARRAY[@]}"; do
        valid=false
        for vs in "${VALID_STATES[@]}"; do
            if [[ "$s" == "$vs" ]]; then
                valid=true
                break
            fi
        done
        
        if [[ "$valid" == false ]]; then
            echo -e "${YELLOW}Warning: '$s' is not in the discovered state list and may cause errors.${NC}"
            invalid_found=true
        fi
    done

    if [[ "$invalid_found" == true ]]; then
        echo ""
        read -p "Continue anyway? (y/n) [n]: " continue_choice
        continue_choice=${continue_choice:-n}
        if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
            echo -e "${RED}Exiting.${NC}"
            exit 1
        fi
    fi
fi

# Salmonella serotype selection - ONLY when SALMONELLA is in the pathogens list
serotype_param=""
if [[ ",$pathogens," == *",SALMONELLA,"* ]] && [[ "$has_serotypes" == true ]] && [[ -n "$preprocessed_metadata" ]]; then
    # Check if serotypes are available with error handling
    if jq -e '.salmonella_serotypes' "$preprocessed_metadata" > /dev/null 2>&1; then
        echo ""
        echo -e "${BLUE}======== Salmonella Serotype Analysis ========${NC}"
        echo "Salmonella was selected. Do you want to:"
        echo "1) Analyze ALL Salmonella serotypes together"
        echo "2) Analyze specific serotypes separately"
        echo "3) Focus on a single top serotype only"
        read -p "Enter selection [1]: " serotype_mode
        serotype_mode=${serotype_mode:-1}
        
        if [[ "$serotype_mode" == "2" ]]; then
            # Display available serotypes (top 20 to keep it manageable) with error handling
            echo ""
            echo -e "${BLUE}Top Salmonella serotypes in dataset:${NC}"
            jq_cmd='.salmonella_serotypes | to_entries | sort_by(.value) | reverse | .[0:20] | .[] | "\(.key): \(.value) isolates"'
            if ! top_serotypes=$(jq -r "$jq_cmd" "$preprocessed_metadata" 2>/dev/null); then
                echo -e "${YELLOW}Could not process serotype information. Showing serotype names only.${NC}"
                jq -r '.salmonella_serotypes | keys | .[0:20]' "$preprocessed_metadata" 2>/dev/null | sed 's/^/- /'
            else
                echo "$top_serotypes" | sed 's/^/- /'
            fi
            
            echo ""
            echo -e "Enter serotypes to analyze (comma-separated with NO spaces)"
            read -p "Serotypes: " serotypes
            
            # Validate serotypes
            if [[ -n "$serotypes" ]]; then
                # Add parameter for serotypes
                serotype_param="--salmonella_serotypes \"$serotypes\""
            else
                echo -e "${YELLOW}No serotypes specified. Analyzing all Salmonella together.${NC}"
                serotype_param=""
            fi
        elif [[ "$serotype_mode" == "3" ]]; then
            # Get top serotype automatically with error handling
            if ! top_serotype=$(jq -r '.salmonella_serotypes | to_entries | sort_by(.value) | reverse | .[0].key' "$preprocessed_metadata" 2>/dev/null); then
                echo -e "${YELLOW}Could not determine top serotype. Analyzing all Salmonella together.${NC}"
                serotype_param=""
            else
                echo -e "${GREEN}Will focus on top serotype: $top_serotype${NC}"
                serotype_param="--salmonella_serotypes \"$top_serotype\""
            fi
        else
            # Analyze all serotypes together (default)
            serotype_param=""
        fi
    else
        echo -e "${YELLOW}No Salmonella serotype information available.${NC}"
        serotype_param=""
    fi
fi

# Set MCMC parameters based on the run mode
if [[ "$flag" == "test" ]]; then
    # Test mode: use minimal settings
    chains=1
    iterations=100
    adapt_delta=0.8
    max_treedepth=8
    
    echo ""
    echo -e "${BLUE}Test Mode: Using minimal settings${NC}"
    echo -e "Chains: ${GREEN}$chains${NC}"
    echo -e "Iterations: ${GREEN}$iterations${NC}"
    echo -e "Adapt delta: ${GREEN}$adapt_delta${NC}"
    echo -e "Max treedepth: ${GREEN}$max_treedepth${NC}"
    
elif [[ "$flag" == "full" ]]; then
    # Full mode: ask for user input
    echo ""
    echo -e "${BLUE}======== MCMC Parameters ========${NC}"
    
    # Number of chains
    echo ""
    read -p "Number of chains [2]: " chains
    chains=${chains:-2}
    
    # Validate chains
    if ! [[ "$chains" =~ ^[0-9]+$ ]]; then
        echo -e "${RED}Invalid input. Using default value (2).${NC}"
        chains=2
    elif [ "$chains" -lt 2 ]; then
        echo -e "${YELLOW}Warning: At least 2 chains recommended for convergence diagnostics.${NC}"
        echo -e "${YELLOW}Continuing with $chains chain(s).${NC}"
    elif [ "$chains" -gt 8 ]; then
        echo -e "${YELLOW}Warning: Large number of chains may significantly increase runtime.${NC}"
        echo -e "${YELLOW}Continuing with $chains chains.${NC}"
    fi
    
    # Number of iterations
    echo ""
    read -p "Number of iterations [500]: " iterations
    iterations=${iterations:-500}
    
    # Validate iterations
    if ! [[ "$iterations" =~ ^[0-9]+$ ]]; then
        echo -e "${RED}Invalid input. Using default value (500).${NC}"
        iterations=500
    elif [ "$iterations" -lt 200 ]; then
        echo -e "${YELLOW}Warning: Low iteration count may lead to poor convergence.${NC}"
        echo -e "${YELLOW}Continuing with $iterations iterations.${NC}"
    elif [ "$iterations" -gt 2000 ]; then
        echo -e "${YELLOW}Warning: High iteration count will significantly increase runtime.${NC}"
        echo -e "${YELLOW}Continuing with $iterations iterations.${NC}"
    fi
    
    # Adapt delta
    echo ""
    read -p "Adapt delta (0.0-1.0) [0.95]: " adapt_delta
    adapt_delta=${adapt_delta:-0.95}
    
    # Validate adapt delta without bc
    if [[ ! "$adapt_delta" =~ ^0?\.[0-9]+$ ]]; then
        echo -e "${RED}Invalid input. Using default value (0.95).${NC}"
        adapt_delta=0.95
    elif [[ "$adapt_delta" == "0.7"* ]] || [[ "$adapt_delta" == "0.6"* ]] || [[ "$adapt_delta" == "0.5"* ]] || [[ "$adapt_delta" == "0.4"* ]] || [[ "$adapt_delta" == "0.3"* ]] || [[ "$adapt_delta" == "0.2"* ]] || [[ "$adapt_delta" == "0.1"* ]] || [[ "$adapt_delta" == "0.0"* ]]; then
        echo -e "${YELLOW}Warning: Low adapt_delta may cause algorithm issues.${NC}"
        echo -e "${YELLOW}Continuing with adapt_delta=$adapt_delta.${NC}"
    elif [[ "$adapt_delta" == "0.99"* ]] || [[ "$adapt_delta" == "1.0"* ]]; then
        echo -e "${YELLOW}Warning: Very high adapt_delta may significantly increase runtime.${NC}"
        echo -e "${YELLOW}Continuing with adapt_delta=$adapt_delta.${NC}"
    fi
    
    # Max treedepth
    echo ""
    read -p "Max treedepth [10]: " max_treedepth
    max_treedepth=${max_treedepth:-10}
    
    # Validate max treedepth
    if ! [[ "$max_treedepth" =~ ^[0-9]+$ ]]; then
        echo -e "${RED}Invalid input. Using default value (10).${NC}"
        max_treedepth=10
    elif [ "$max_treedepth" -lt 8 ]; then
        echo -e "${YELLOW}Warning: Low max_treedepth may cause truncated trajectories.${NC}"
        echo -e "${YELLOW}Continuing with max_treedepth=$max_treedepth.${NC}"
    elif [ "$max_treedepth" -gt 15 ]; then
        echo -e "${YELLOW}Warning: High max_treedepth will significantly increase runtime.${NC}"
        echo -e "${YELLOW}Continuing with max_treedepth=$max_treedepth.${NC}"
    fi
    
elif [[ "$flag" == "resume" ]]; then
    # Resume mode: use default values for command construction
    # (actual values will come from the resume cache)
    chains=2
    iterations=500
    adapt_delta=0.95
    max_treedepth=10
    
    echo ""
    echo -e "${BLUE}Resume Mode: Using parameters from previous run${NC}"
fi

# Build the command with appropriate parameters
if [[ "$flag" == "resume" ]]; then
    cmd="nextflow run main.nf -profile singularity -resume -entry SPLINE \
  --censusFileB \"$censusFileB\" \
  --censusFileP \"$censusFileP\" \
  --travel \"$travel\" \
  --cidt \"$cidt\" \
  --iterations $iterations \
  --chains $chains \
  --adapt_delta $adapt_delta \
  --max_treedepth $max_treedepth \
  --seed 123 \
  --outdir \"$outDir\" \
  --pathogen \"$pathogens\" \
  --states \"$states\" \
  $serotype_param"
else
    cmd="nextflow run main.nf -profile singularity -entry SPLINE \
  --censusFileB \"$censusFileB\" \
  --censusFileP \"$censusFileP\" \
  --travel \"$travel\" \
  --cidt \"$cidt\" \
  --iterations $iterations \
  --chains $chains \
  --adapt_delta $adapt_delta \
  --max_treedepth $max_treedepth \
  --seed 123 \
  --outdir \"$outDir\" \
  --pathogen \"$pathogens\" \
  --states \"$states\" \
  $serotype_param"
fi

# Add mmwrFile parameter, which could be either raw SAS or preprocessed CSV
if [[ -n "$mmwrFile" ]]; then
    cmd="$cmd --mmwrFile \"$mmwrFile\""
fi

# Add preprocessed flag if using preprocessed data
if [[ "$workflow_mode" == "1" || "$workflow_mode" == "3" ]]; then
    cmd="$cmd --preprocessed true"
fi

# Add metadata parameter if available
if [[ -n "$preprocessed_metadata" && -f "$preprocessed_metadata" ]]; then
    cmd="$cmd --metadata \"$preprocessed_metadata\""
fi

# Add background option if requested
if [[ $background == true ]]; then
    log_file="foodnet_run_${timestamp}.log"
    bg_cmd="nohup $cmd > $log_file 2>&1 &"
    final_cmd="$bg_cmd"
    echo -e "${YELLOW}Process will run in background with log: $log_file${NC}"
else
    final_cmd="$cmd"
fi

# Review and confirm
echo ""
echo -e "${BLUE}========= Analysis Summary ===========${NC}"
# Fixed conditional for mode:
if [ "$flag" == "test" ]; then
  mode_string="Test run"
elif [ "$flag" == "full" ]; then
  mode_string="Full analysis"
else
  mode_string="Resume previous run"
fi
echo -e "Mode: ${GREEN}${mode_string}${NC}"

# Fixed conditional for preprocessed data:
if [ "$workflow_mode" == "1" ] || [ "$workflow_mode" == "3" ]; then
  echo -e "Using preprocessed data: ${GREEN}Yes${NC}"
else
  echo -e "Using preprocessed data: ${GREEN}No${NC}"
fi

echo -e "Pathogens: ${GREEN}$pathogens${NC}"
echo -e "States: ${GREEN}$states${NC}"echo -e "Pathogens: ${GREEN}$pathogens${NC}"
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

read -p "Proceed with analysis? (y/n) [y]: " proceed
proceed=${proceed:-y}

if [[ "$proceed" =~ ^[Yy]$ ]]; then
    echo -e "${GREEN}Starting analysis...${NC}"
    
    # Create the output directory if it doesn't exist
    mkdir -p "$outDir"
    
    # Run the command
    if eval $final_cmd; then
        if [[ $background == true ]]; then
            echo -e "${GREEN}Process started in background. Check status with:${NC}"
            echo -e "${YELLOW}tail -f $log_file${NC}"
        else
            echo -e "${GREEN}Analysis completed successfully.${NC}"
            echo -e "${GREEN}Results are available in: $outDir${NC}"
        fi
    else
        echo -e "${RED}Error running analysis command.${NC}"
        echo -e "${RED}Check .nextflow.log for details.${NC}"
        exit 1
    fi
else
    echo -e "${RED}Analysis cancelled.${NC}"
fi
