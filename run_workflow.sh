#!/bin/bash

# set up paths, files
dataDir="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/"
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

# Define all available pathogens
ALL_PATHOGENS="CAMPYLOBACTER,CYCLOSPORA,SALMONELLA,SHIGELLA,STEC,VIBRIO,YERSINIA"

# Initialize optional variables to prevent undefined errors
serotype_config=""
catchment_config=""
matching_sensitivity="MEDIUM"

# Function to handle pathogen grouping decisions
handle_pathogen_grouping() {
    local pathogen="$1"
    local preprocessed_file="$2"
    local grouping=""
    
    case "$pathogen" in
        STEC)
            echo "" >&2
            echo -e "${BLUE}STEC Grouping Options:${NC}" >&2
            echo "STEC can be analyzed as:" >&2
            echo "1) Combined - All STEC together" >&2
            echo "2) O157 only - Just STEC O157" >&2
            echo "3) Non-O157 only - Just non-O157 STEC" >&2
            echo "4) Both separately - O157 and non-O157 as separate analyses" >&2
            echo "" >&2
            read -p "Select STEC grouping option [1]: " stec_choice
            stec_choice=${stec_choice:-1}
            
            case "$stec_choice" in
                1) grouping="STEC:combined" ;;
                2) grouping="STEC:O157" ;;
                3) grouping="STEC:nonO157" ;;
                4) grouping="STEC:O157|STEC:nonO157" ;;
                *) 
                    echo -e "${YELLOW}Invalid choice. Using combined.${NC}" >&2
                    grouping="STEC:combined" 
                    ;;
            esac
            ;;
        
        SALMONELLA)
            echo "" >&2
            echo -e "${BLUE}Salmonella Grouping Options:${NC}" >&2
            echo "1) Combined - All serotypes together" >&2
            echo "2) Custom selection - Choose specific serotypes" >&2
            echo "" >&2
            read -p "Select Salmonella grouping option [1]: " sal_choice
            sal_choice=${sal_choice:-1}
            
            if [[ "$sal_choice" == "2" ]]; then
                # Extract and rank serotypes
                echo "" >&2
                echo "Analyzing serotypes in data..." >&2
                
                # Extract serotypes with counts from preprocessed data
                # This command extracts the serotypesummary column and counts occurrences
                serotype_data=$(awk -F',' -v pathogen="SALMONELLA" '
                    BEGIN { OFS="\t" }
                    NR==1 { 
                        for(i=1; i<=NF; i++) {
                            gsub(/^"|"$/, "", $i)
                            if($i == "pathogen") p_col=i
                            if($i == "serotypesummary") s_col=i
                        }
                    }
                    NR>1 && p_col && s_col {
                        gsub(/^"|"$/, "", $p_col)
                        gsub(/^"|"$/, "", $s_col)
                        if($p_col == pathogen && $s_col != "") {
                            serotypes[$s_col]++
                        }
                    }
                    END {
                        for(s in serotypes) {
                            print serotypes[s], s
                        }
                    }
                ' "$preprocessed_file" | sort -rn)
                
                if [[ -z "$serotype_data" ]]; then
                    echo -e "${YELLOW}No serotype data found. Using combined analysis.${NC}" >&2
                    grouping="SALMONELLA:combined"
                else
                    # Display ranked serotypes
                    echo -e "${GREEN}Top serotypes found:${NC}" >&2
                    echo "$serotype_data" | head -20 | nl -nln -w3 | while read num count serotype; do
                        printf "[%s] %-30s (n=%s)\n" "$num" "$serotype" "$count" >&2
                    done
                    
                    # Check if there are more
                    total_serotypes=$(echo "$serotype_data" | wc -l)
                    if [[ $total_serotypes -gt 20 ]]; then
                        echo "" >&2
                        echo "... and $((total_serotypes - 20)) more serotypes" >&2
                        read -p "Show all serotypes? (y/n) [n]: " show_all
                        show_all=${show_all:-n}
                        if [[ "$show_all" =~ ^[Yy]$ ]]; then
                            echo "$serotype_data" | tail -n +21 | nl -nln -w3 -v 21 | while read num count serotype; do
                                printf "[%s] %-30s (n=%s)\n" "$num" "$serotype" "$count" >&2
                            done
                        fi
                    fi
                    
                    echo "" >&2
                    echo "Enter serotype numbers to analyze (comma-separated, e.g., 1,2,5)" >&2
                    echo "Or press Enter to analyze all serotypes combined" >&2
                    read -p "Selection: " serotype_selection
                    
                    if [[ -z "$serotype_selection" ]]; then
                        grouping="SALMONELLA:combined"
                    else
                        # Convert numbers to serotype names
                        selected_serotypes=""
                        IFS=',' read -ra selections <<< "$serotype_selection"
                        for sel in "${selections[@]}"; do
                            serotype_name=$(echo "$serotype_data" | sed -n "${sel}p" | cut -f2-)
                            if [[ -n "$serotype_name" ]]; then
                                if [[ -n "$selected_serotypes" ]]; then
                                    selected_serotypes="${selected_serotypes}|SALMONELLA:${serotype_name}"
                                else
                                    selected_serotypes="SALMONELLA:${serotype_name}"
                                fi
                            fi
                        done
                        # Check if any valid serotypes were selected
                        if [[ -z "$selected_serotypes" ]]; then
                            echo -e "${YELLOW}No valid serotypes selected. Using combined analysis.${NC}" >&2
                            grouping="SALMONELLA:combined"
                        else
                            grouping="$selected_serotypes"
                        fi
                    fi
                fi
            else
                grouping="SALMONELLA:combined"
            fi
            ;;
            
        *)
            # For other pathogens, analyze as single group
            grouping="${pathogen}:combined"
            ;;
    esac
    
    echo "$grouping"
}

echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}   FoodNet Trends Analysis Pipeline      ${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

# Ask for run mode
echo -e "Select run mode:"
echo "1) Test"
echo "2) Publication" 
echo "3) Max"
echo "4) Custom"
echo "5) Resume previous run"
read -p "Enter selection [1]: " run_mode
run_mode=${run_mode:-1}

# Validate run mode
if [[ ! "$run_mode" =~ ^[1-5]$ ]]; then
    echo -e "${RED}Invalid selection. Using default (Test).${NC}"
    run_mode=1
fi

# Convert run mode to flag
case $run_mode in
    1) flag="test" ;;
    2) flag="publication" ;;
    3) flag="max" ;;
    4) flag="custom" ;;
    5) flag="resume" ;;
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

# Check for existing preprocessed data FIRST (before pathogen selection)
# This feature allows reusing previously cleaned data to save processing time
# Searches for clean_mmwr.csv files in output directories from previous runs
echo ""
echo -e "${BLUE}======== Preprocessing Options ========${NC}"
echo "Checking for existing preprocessed data..."

# Search for preprocessed files in output directories
# Uses find with -print0 for safe handling of filenames with spaces
# Sorted by path for consistent ordering
preprocessed_files=()
if [[ -d "output" ]]; then
    while IFS= read -r -d '' file; do
        preprocessed_files+=("$file")
    done < <(find "output" -name "clean_mmwr.csv" -type f -print0 2>/dev/null | sort -z)
fi

use_preprocessed=false
preprocessed_file=""

if [[ ${#preprocessed_files[@]} -gt 0 ]]; then
    echo ""
    echo -e "${GREEN}Found existing preprocessed data:${NC}"
    echo ""
    
    # Display options with metadata
    for i in "${!preprocessed_files[@]}"; do
        file="${preprocessed_files[$i]}"
        file_date=$(stat -c %y "$file" 2>/dev/null | cut -d' ' -f1,2 | cut -d'.' -f1)
        file_size=$(du -h "$file" 2>/dev/null | cut -f1)
        echo "$((i+1))) $file"
        echo "    Created: $file_date, Size: $file_size"
        
        # Check if preprocessing report exists
        report_file="${file%.csv}_preprocessing_report.csv"
        if [[ -f "$report_file" ]]; then
            echo "    ✓ Preprocessing report available"
        fi
        echo ""
    done
    
    echo "0) Skip - run preprocessing again"
    echo ""
    read -p "Select preprocessed file to use [0]: " selection
    selection=${selection:-0}
    
    if [[ "$selection" =~ ^[1-9][0-9]*$ ]] && [[ $selection -le ${#preprocessed_files[@]} ]]; then
        use_preprocessed=true
        preprocessed_file="${preprocessed_files[$((selection-1))]}"
        echo -e "${GREEN}Using preprocessed data: $preprocessed_file${NC}"
    else
        echo -e "${YELLOW}Will run preprocessing step${NC}"
    fi
else
    echo -e "${YELLOW}No preprocessed data found. Will run preprocessing step.${NC}"
fi

# NOW ask about pathogens (after we know about preprocessed data)
echo ""
echo -e "Pathogen selection:"
echo "1) Run ALL available pathogens"
echo "2) Select specific pathogens"
read -p "Enter selection [2]: " pathogen_mode
pathogen_mode=${pathogen_mode:-2}

if [[ "$pathogen_mode" == "1" ]]; then
    # Use all pathogens
    if [[ "$use_preprocessed" == true ]] && [[ -f "$preprocessed_file" ]]; then
        # Extract all pathogens from preprocessed data
        echo "Extracting all pathogens from preprocessed data..."
        all_from_data=$(cut -d',' -f1 "$preprocessed_file" | tail -n +2 | \
            sed 's/^"//;s/"$//' | \
            grep '^[A-Z][A-Z]*$' | \
            sort -u | \
            tr '\n' ',' | sed 's/,$//')
        if [[ -n "$all_from_data" ]]; then
            pathogens=$all_from_data
            echo -e "${GREEN}Selected: ALL pathogens found in data (${pathogens})${NC}"
        else
            pathogens=$ALL_PATHOGENS
            echo -e "${YELLOW}Could not extract pathogens from data. Using default list.${NC}"
            echo -e "${GREEN}Selected: ALL pathogens (${ALL_PATHOGENS})${NC}"
        fi
    else
        pathogens=$ALL_PATHOGENS
        echo -e "${GREEN}Selected: ALL pathogens (${ALL_PATHOGENS})${NC}"
    fi
else
    # Ask for specific pathogens
    echo ""
    
    # If using preprocessed data, try to extract available pathogens
    if [[ "$use_preprocessed" == true ]] && [[ -f "$preprocessed_file" ]]; then
        echo "Analyzing preprocessed data for available pathogens..."
        # More robust pathogen extraction that handles quoted CSV fields
        # Use cut to get first column, then clean and filter
        available_pathogens=$(cut -d',' -f1 "$preprocessed_file" | tail -n +2 | \
            sed 's/^"//;s/"$//' | \
            grep '^[A-Z][A-Z]*$' | \
            sort -u | \
            tr '\n' ',' | sed 's/,$//')
        
        if [[ -n "$available_pathogens" ]]; then
            echo -e "${GREEN}Found pathogens in preprocessed data:${NC}"
            IFS=',' read -ra pathogen_array <<< "$available_pathogens"
            for p in "${pathogen_array[@]}"; do
                echo "- $p"
            done
            echo ""
            echo "Enter pathogens to analyze (comma-separated with NO spaces)"
            echo -e "${YELLOW}Available: $available_pathogens${NC}"
            read -p "Leave blank to analyze all found pathogens: " pathogens
            pathogens=${pathogens:-"$available_pathogens"}
            # Ensure pathogens is not empty
            if [[ -z "$pathogens" ]]; then
                echo -e "${YELLOW}No pathogens specified. Using default.${NC}"
                pathogens="CAMPYLOBACTER,CYCLOSPORA"
            fi
        else
            echo -e "${YELLOW}Could not extract pathogen list from preprocessed data.${NC}"
            echo "Using standard pathogen list..."
            echo "- CAMPYLOBACTER"
            echo "- CYCLOSPORA"
            echo "- SALMONELLA"
            echo "- SHIGELLA"
            echo "- STEC"
            echo "- VIBRIO"
            echo "- YERSINIA"
            echo ""
            echo "Enter pathogens to analyze (comma-separated with NO spaces)"
            read -p "Leave blank for default (CAMPYLOBACTER,CYCLOSPORA): " pathogens
            pathogens=${pathogens:-"CAMPYLOBACTER,CYCLOSPORA"}
        fi
    else
        echo -e "Available pathogens:"
        echo "- CAMPYLOBACTER"
        echo "- CYCLOSPORA"
        echo "- SALMONELLA"
        echo "- SHIGELLA"
        echo "- STEC"
        echo "- VIBRIO"
        echo "- YERSINIA"
        echo ""
        echo "Enter pathogens to analyze (comma-separated with NO spaces)"
        read -p "Leave blank for default (CAMPYLOBACTER,CYCLOSPORA): " pathogens
        pathogens=${pathogens:-"CAMPYLOBACTER,CYCLOSPORA"}
    fi

    # Validate pathogens
    if [[ "$use_preprocessed" == false ]]; then
        # Strict validation when not using preprocessed data
        valid_pathogens=("CAMPYLOBACTER" "CYCLOSPORA" "SALMONELLA" "SHIGELLA" "STEC" "VIBRIO" "YERSINIA")
        IFS=',' read -ra pathogen_array <<< "$pathogens"
        invalid_found=false

        for p in "${pathogen_array[@]}"; do
            valid=false
            for vp in "${valid_pathogens[@]}"; do
                if [[ "$p" == "$vp" ]]; then
                    valid=true
                    break
                fi
            done
            
            if [[ "$valid" == false ]]; then
                echo -e "${YELLOW}Warning: '$p' is not a recognized pathogen and may cause errors.${NC}"
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
    else
        # More lenient validation when using preprocessed data
        # Just check that something was entered
        if [[ -z "$pathogens" ]]; then
            echo -e "${YELLOW}No pathogens selected. Using default.${NC}"
            pathogens="CAMPYLOBACTER,CYCLOSPORA"
        fi
        echo -e "${GREEN}Using pathogens from preprocessed data: $pathogens${NC}"
    fi
fi

# The preprocessing selection has been moved earlier in the script

# Ask about data configuration files
echo ""
echo -e "${BLUE}======== Data Configuration Options ========${NC}"
echo "These are optional CSV files for custom analysis settings."
echo ""

# Serotype configuration
read -p "Use custom serotype configuration? (y/n) [n]: " use_serotype_config
use_serotype_config=${use_serotype_config:-n}
serotype_config=""
if [[ "$use_serotype_config" =~ ^[Yy]$ ]]; then
    echo "Enter path to serotype configuration CSV file"
    echo "(Example: analysis_configs/examples/serotype_config.csv)"
    read -p "Path: " serotype_config
    if [[ -n "$serotype_config" ]] && [[ ! -f "$serotype_config" ]]; then
        echo -e "${YELLOW}Warning: File not found: $serotype_config${NC}"
        echo -e "${YELLOW}Pipeline will fail if file doesn't exist at runtime.${NC}"
    fi
fi

# Catchment configuration
echo ""
read -p "Use custom catchment configuration? (y/n) [n]: " use_catchment_config
use_catchment_config=${use_catchment_config:-n}
catchment_config=""
if [[ "$use_catchment_config" =~ ^[Yy]$ ]]; then
    echo "Enter path to catchment configuration CSV file"
    echo "(Example: analysis_configs/examples/catchment_config.csv)"
    read -p "Path: " catchment_config
    if [[ -n "$catchment_config" ]] && [[ ! -f "$catchment_config" ]]; then
        echo -e "${YELLOW}Warning: File not found: $catchment_config${NC}"
        echo -e "${YELLOW}Pipeline will fail if file doesn't exist at runtime.${NC}"
    fi
fi

# Pathogen matching sensitivity (only if not using preprocessed data)
# matching_sensitivity already initialized at top with default "MEDIUM"
if [[ "$use_preprocessed" == false ]]; then
    echo ""
    echo -e "${BLUE}======== Pathogen Name Standardization ========${NC}"
    echo "Choose sensitivity level for pathogen name matching:"
    echo "  STRICT  - Exact matches only"
    echo "  MEDIUM  - Exact + prefix matching + fuzzy (1 char diff)"
    echo "  RELAXED - All matching methods + fuzzy (2 char diff)"
    echo ""
    read -p "Matching sensitivity [MEDIUM]: " matching_sensitivity
    matching_sensitivity=${matching_sensitivity:-MEDIUM}
    # Validate input
    if [[ ! "$matching_sensitivity" =~ ^(STRICT|MEDIUM|RELAXED)$ ]]; then
        echo -e "${YELLOW}Invalid input. Using default (MEDIUM).${NC}"
        matching_sensitivity="MEDIUM"
    fi
fi

# Set MCMC parameters based on the run mode
if [[ "$flag" == "test" ]]; then
    chains=1
    iterations=100
    adapt_delta=0.8
    max_treedepth=8
    
    echo ""
    echo -e "${BLUE}Test Profile${NC}"
    echo -e "Chains: ${GREEN}$chains${NC}"
    echo -e "Iterations: ${GREEN}$iterations${NC}"
    echo -e "Adapt delta: ${GREEN}$adapt_delta${NC}"
    echo -e "Max treedepth: ${GREEN}$max_treedepth${NC}"
    
elif [[ "$flag" == "publication" ]]; then
    chains=6
    iterations=10001
    adapt_delta=0.99
    max_treedepth=15
    
    echo ""
    echo -e "${BLUE}Publication Profile${NC}"
    echo -e "Chains: ${GREEN}$chains${NC}"
    echo -e "Iterations: ${GREEN}$iterations${NC}"
    echo -e "Adapt delta: ${GREEN}$adapt_delta${NC}"
    echo -e "Max treedepth: ${GREEN}$max_treedepth${NC}"
    
elif [[ "$flag" == "max" ]]; then
    chains=8
    iterations=20000
    adapt_delta=0.99
    max_treedepth=15
    
    echo ""
    echo -e "${BLUE}Max Profile${NC}"
    echo -e "Chains: ${GREEN}$chains${NC}"
    echo -e "Iterations: ${GREEN}$iterations${NC}"
    echo -e "Adapt delta: ${GREEN}$adapt_delta${NC}"
    echo -e "Max treedepth: ${GREEN}$max_treedepth${NC}"
    
elif [[ "$flag" == "custom" ]]; then
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
    echo -e "${YELLOW}Note: Resume mode uses settings from the previous run.${NC}"
    echo -e "${YELLOW}Configuration files and pathogen groupings are preserved.${NC}"
    
    # Need to ensure pathogens is set for resume mode
    if [[ -z "$pathogens" ]]; then
        echo -e "${YELLOW}Warning: No pathogens specified for resume. Using default.${NC}"
        pathogens="CAMPYLOBACTER,CYCLOSPORA"
    fi
fi

# Handle pathogen grouping for STEC and Salmonella
# This section processes user preferences for how to analyze pathogen subgroups
pathogen_grouping=""
if [[ "$flag" != "resume" ]]; then
    # Check if STEC or Salmonella are in the selected pathogens
    IFS=',' read -ra selected_pathogens <<< "$pathogens"
    grouped_pathogens=""
    
    for pathogen in "${selected_pathogens[@]}"; do
        if [[ "$pathogen" == "STEC" ]] || [[ "$pathogen" == "SALMONELLA" ]]; then
            # Get the preprocessed file to use for analysis
            if [[ "$use_preprocessed" == true ]]; then
                data_file="$preprocessed_file"
            else
                # If not using preprocessed, still ask for grouping preferences
                # For STEC, we can offer options without seeing the data
                # For Salmonella, we'll default to combined unless using preprocessed
                data_file=""
            fi
            
            if [[ -n "$data_file" ]] && [[ -f "$data_file" ]]; then
                # Can analyze the file for serotypes
                grouping=$(handle_pathogen_grouping "$pathogen" "$data_file")
            else
                # Handle without data file
                if [[ "$pathogen" == "STEC" ]]; then
                    # STEC grouping can be decided without seeing data
                    grouping=$(handle_pathogen_grouping "$pathogen" "")
                else
                    # For Salmonella without data, default to combined
                    echo "" >&2
                    echo -e "${YELLOW}Note: Salmonella serotype selection requires preprocessed data.${NC}" >&2
                    echo -e "${YELLOW}Using combined analysis for all Salmonella serotypes.${NC}" >&2
                    grouping="SALMONELLA:combined"
                fi
            fi
            
            # Add to grouped pathogens
            if [[ -n "$grouped_pathogens" ]]; then
                grouped_pathogens="${grouped_pathogens}|$grouping"
            else
                grouped_pathogens="$grouping"
            fi
        else
            # Non-grouped pathogens stay as-is
            if [[ -n "$grouped_pathogens" ]]; then
                grouped_pathogens="${grouped_pathogens}|${pathogen}:combined"
            else
                grouped_pathogens="${pathogen}:combined"
            fi
        fi
    done
    
    # Store the grouping configuration
    pathogen_grouping="$grouped_pathogens"
fi

# Build the base command
cmd="nextflow run main.nf -profile singularity -entry SPLINE \
  --mmwrFile \"$dataDir/mmwr9624_May2025.sas7bdat\" \
  --censusFileB \"$dataDir/cen9624.sas7bdat\" \
  --censusFileP \"$dataDir/cen9624_para.sas7bdat\" \
  --iterations $iterations \
  --chains $chains \
  --adapt_delta $adapt_delta \
  --max_treedepth $max_treedepth \
  --seed 123 \
  --outdir \"$outDir\" \
  --pathogen \"$pathogens\" \
  --projID \"$timestamp\""

# Add resume flag if in resume mode
if [[ "$flag" == "resume" ]]; then
    cmd="nextflow run main.nf -profile singularity -resume -entry SPLINE \
  --mmwrFile \"$dataDir/mmwr9624_May2025.sas7bdat\" \
  --censusFileB \"$dataDir/cen9624.sas7bdat\" \
  --censusFileP \"$dataDir/cen9624_para.sas7bdat\" \
  --iterations $iterations \
  --chains $chains \
  --adapt_delta $adapt_delta \
  --max_treedepth $max_treedepth \
  --seed 123 \
  --outdir \"$outDir\" \
  --pathogen \"$pathogens\" \
  --projID \"$timestamp\""
fi

# Add optional parameters (same for both resume and non-resume)
# Add configuration files if provided
if [[ -n "$serotype_config" ]]; then
    cmd="$cmd --serotype_config \"$serotype_config\""
fi
if [[ -n "$catchment_config" ]]; then
    cmd="$cmd --catchment_config \"$catchment_config\""
fi
# Add matching sensitivity
cmd="$cmd --matching_sensitivity \"$matching_sensitivity\""
# Add preprocessed data flags if using
if [[ "$use_preprocessed" == true ]]; then
    cmd="$cmd --preprocessed true --cleanFile \"$preprocessed_file\""
fi
# Add pathogen grouping if specified
if [[ -n "$pathogen_grouping" ]]; then
    cmd="$cmd --pathogen_grouping \"$pathogen_grouping\""
fi

# Add background option if requested
if [[ $background == true ]]; then
  bg_cmd="nohup $cmd > foodnet_run_${timestamp}.log 2>&1 &"
  final_cmd="$bg_cmd"
  echo -e "${YELLOW}Process will run in background with log: foodnet_run_${timestamp}.log${NC}"
else
  final_cmd="$cmd"
fi

# Review and confirm
echo ""
echo -e "${BLUE}========= Analysis Summary ==========${NC}"
echo -e "Mode: ${GREEN}$([ "$flag" == "test" ] && echo "Test" || [ "$flag" == "publication" ] && echo "Publication" || [ "$flag" == "max" ] && echo "Max" || [ "$flag" == "custom" ] && echo "Custom" || echo "Resume previous run")${NC}"
echo -e "Pathogens: ${GREEN}$pathogens${NC}"

if [[ "$flag" != "resume" ]]; then
    echo -e "Chains: ${GREEN}$chains${NC}"
    echo -e "Iterations: ${GREEN}$iterations${NC}"
    echo -e "Adapt delta: ${GREEN}$adapt_delta${NC}"
    echo -e "Max treedepth: ${GREEN}$max_treedepth${NC}"
fi

echo -e "Output directory: ${GREEN}$outDir${NC}"
echo -e "Run in background: ${GREEN}$([ "$background" == true ] && echo "Yes" || echo "No")${NC}"

if [[ "$use_preprocessed" == true ]]; then
    echo -e "Preprocessed data: ${GREEN}Yes - $preprocessed_file${NC}"
else
    echo -e "Preprocessed data: ${GREEN}No - will run preprocessing${NC}"
    echo -e "Pathogen matching: ${GREEN}$matching_sensitivity${NC}"
fi

if [[ -n "$serotype_config" ]] || [[ -n "$catchment_config" ]]; then
    echo ""
    echo -e "${BLUE}Data Configuration Files:${NC}"
    if [[ -n "$serotype_config" ]]; then
        echo -e "Serotype config: ${GREEN}$serotype_config${NC}"
    fi
    if [[ -n "$catchment_config" ]]; then
        echo -e "Catchment config: ${GREEN}$catchment_config${NC}"
    fi
fi
echo ""
echo -e "Command to run:"
echo -e "${YELLOW}$final_cmd${NC}"
echo ""

read -p "Proceed with analysis? (y/n) [y]: " proceed
proceed=${proceed:-y}

if [[ "$proceed" =~ ^[Yy]$ ]]; then
    echo -e "${GREEN}Starting analysis...${NC}"
    eval $final_cmd
    
    if [[ $background == true ]]; then
        echo -e "${GREEN}Process started in background. Check status with:${NC}"
        echo -e "${YELLOW}tail -f foodnet_run_${timestamp}.log${NC}"
    fi
else
    echo -e "${RED}Analysis cancelled.${NC}"
fi
