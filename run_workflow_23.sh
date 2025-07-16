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
pathogens=""  # Initialize to prevent undefined variable errors
use_preprocessed=false
preprocessed_file=""

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
            
            while true; do
                read -p "Select STEC grouping option [1]: " stec_choice
                stec_choice=${stec_choice:-1}
                
                case "$stec_choice" in
                    1) grouping="STEC:combined"; break ;;
                    2) grouping="STEC:O157"; break ;;
                    3) grouping="STEC:nonO157"; break ;;
                    4) grouping="STEC:O157|STEC:nonO157"; break ;;
                    *) 
                        echo -e "${RED}Invalid choice: '$stec_choice'. Please enter 1-4.${NC}" >&2
                        ;;
                esac
            done
            ;;
        
        SALMONELLA)
            echo "" >&2
            echo -e "${BLUE}Salmonella Grouping Options:${NC}" >&2
            echo "1) Combined - All serotypes together" >&2
            echo "2) Custom selection - Choose specific serotypes" >&2
            echo "" >&2
            
            while true; do
                read -p "Select Salmonella grouping option [1]: " sal_choice
                sal_choice=${sal_choice:-1}
                
                if [[ "$sal_choice" =~ ^[12]$ ]]; then
                    break
                else
                    echo -e "${RED}Invalid choice: '$sal_choice'. Please enter 1 or 2.${NC}" >&2
                fi
            done
            
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
                        
                        while true; do
                            read -p "Show all serotypes? (y/n) [n]: " show_all
                            show_all=${show_all:-n}
                            
                            if [[ "$show_all" =~ ^[YyNn]$ ]]; then
                                break
                            else
                                echo -e "${RED}Invalid input: '$show_all'. Please enter y or n.${NC}" >&2
                            fi
                        done
                        
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
    
    # Validate that grouping is not empty
    if [[ -z "$grouping" ]] || [[ -z "${grouping// }" ]]; then
        echo -e "${YELLOW}Warning: Empty grouping detected. Using default.${NC}" >&2
        grouping="${pathogen}:combined"
    fi
    
    echo "$grouping"
}

# Function to extract states from metadata file
extract_states_from_data() {
    local metadata_dir=$1
    local states_file="${metadata_dir}/metadata_states.csv"
    
    if [[ -f "$states_file" ]]; then
        # Read and format state data
        awk -F',' 'NR>1 {
            printf "%s,%d,%d,%d\n", $1, $2, $3, $4
        }' "$states_file"
        return 0
    else
        return 1
    fi
}

# Function to extract CIDT methods from metadata file
extract_cidt_from_data() {
    local metadata_dir=$1
    local cidt_file="${metadata_dir}/metadata_cidt.csv"
    
    if [[ -f "$cidt_file" ]]; then
        # Read and format CIDT data
        awk -F',' 'NR>1 {
            printf "%s,%d,%.1f\n", $1, $2, $5
        }' "$cidt_file"
        return 0
    else
        return 1
    fi
}

# Function to extract travel statuses from metadata file
extract_travel_from_data() {
    local metadata_dir=$1
    local travel_file="${metadata_dir}/metadata_travel.csv"
    
    if [[ -f "$travel_file" ]]; then
        # Read and format travel data
        awk -F',' 'NR>1 {
            printf "%s,%d,%.1f\n", $1, $2, $3
        }' "$travel_file"
        return 0
    else
        return 1
    fi
}

# Function to display state names
get_state_name() {
    case $1 in
        CA) echo "California" ;;
        CO) echo "Colorado" ;;
        CT) echo "Connecticut" ;;
        GA) echo "Georgia" ;;
        MD) echo "Maryland" ;;
        MN) echo "Minnesota" ;;
        NM) echo "New Mexico" ;;
        NY) echo "New York" ;;
        OR) echo "Oregon" ;;
        TN) echo "Tennessee" ;;
        *) echo "$1" ;;
    esac
}

echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}   FoodNet Trends Analysis Pipeline      ${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

# Ask for run mode
while true; do
    echo -e "Select run mode:"
    echo "1) Test"
    echo "2) Publication" 
    echo "3) Max"
    echo "4) Custom"
    echo "5) Resume previous run"
    read -p "Enter selection [1]: " run_mode
    run_mode=${run_mode:-1}
    
    # Validate run mode
    if [[ "$run_mode" =~ ^[1-5]$ ]]; then
        break
    else
        echo -e "${RED}Invalid selection: '$run_mode'. Please enter 1-5.${NC}"
        echo ""
    fi
done

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
while true; do
    read -p "Run in background? (y/n) [n]: " bg_choice
    bg_choice=${bg_choice:-n}
    
    if [[ "$bg_choice" =~ ^[YyNn]$ ]]; then
        break
    else
        echo -e "${RED}Invalid input: '$bg_choice'. Please enter y or n.${NC}"
    fi
done

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
        
        # Check if resource profile exists
        resource_profile="$(dirname "$file")/resource_profile.csv"
        if [[ -f "$resource_profile" ]]; then
            echo "    ✓ Resource profile available"
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
        # When using preprocessed data with "ALL pathogens", use AUTO_DISCOVER
        echo -e "${GREEN}Using AUTO_DISCOVER mode to extract all pathogens from preprocessed data${NC}"
        pathogens="AUTO_DISCOVER"
        pathogen_grouping=""  # Clear any grouping - will be handled by workflow
    elif [[ "$use_preprocessed" == true ]]; then
        # This shouldn't happen but handle it
        echo -e "${YELLOW}Preprocessed file not found. Using default pathogen list.${NC}"
        pathogens=$ALL_PATHOGENS
        echo -e "${GREEN}Selected: ALL pathogens (${ALL_PATHOGENS})${NC}"
    else
        # Not using preprocessed data - use AUTO_DISCOVER
        echo -e "${GREEN}Using AUTO_DISCOVER mode for fresh data processing${NC}"
        pathogens="AUTO_DISCOVER"
        pathogen_grouping=""  # Clear any grouping - will be handled by workflow
    fi
else
    # Ask for specific pathogens
    echo ""
    
    # If using preprocessed data, try to extract available pathogens
    if [[ "$use_preprocessed" == true ]] && [[ -f "$preprocessed_file" ]]; then
        echo "Analyzing preprocessed data for available pathogens..."
        # More robust pathogen extraction that handles quoted CSV fields
        # Use awk for proper CSV parsing, handling quoted fields
        available_pathogens=$(awk -F',' '
            function unquote(s) {
                gsub(/^"/, "", s)
                gsub(/"$/, "", s)
                gsub(/""/, "\"", s)  # Handle escaped quotes
                return s
            }
            NR>1 {
                # Get first field, handling quoted values
                field1 = $1
                # If the field starts with a quote, we need to handle embedded commas
                if (substr(field1, 1, 1) == "\"") {
                    # Find the closing quote
                    full_field = field1
                    for (i = 2; i <= NF; i++) {
                        full_field = full_field "," $i
                        if (substr($i, length($i), 1) == "\"" && substr($i, length($i)-1, 1) != "\"") {
                            break
                        }
                    }
                    field1 = full_field
                }
                pathogen = unquote(field1)
                # Only include valid pathogen names (starts with letter, allows letters, numbers, dots, hyphens, spaces)
                # Must be primarily uppercase but allow some flexibility
                if (pathogen ~ /^[A-Za-z][A-Za-z0-9. -]*$/ && length(pathogen) > 0) {
                    # Convert to uppercase for consistency
                    pathogen = toupper(pathogen)
                    print pathogen
                }
            }' "$preprocessed_file" | sort -u | tr '\n' ',' | sed 's/,$//')
        
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
            # Ensure pathogens is not empty after defaulting
            if [[ -z "$pathogens" ]] || [[ -z "${pathogens// }" ]]; then
                echo -e "${YELLOW}No pathogens could be determined. Using default.${NC}"
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

    # Validate and normalize pathogen names
    if [[ "$use_preprocessed" == false ]]; then
        # Strict validation when not using preprocessed data
        valid_pathogens=("CAMPYLOBACTER" "CYCLOSPORA" "SALMONELLA" "SHIGELLA" "STEC" "VIBRIO" "YERSINIA")
        IFS=',' read -ra pathogen_array <<< "$pathogens"
        normalized_pathogens=""
        invalid_found=false

        for p in "${pathogen_array[@]}"; do
            # Convert to uppercase and trim whitespace
            p_upper=$(echo "$p" | tr '[:lower:]' '[:upper:]' | xargs)
            valid=false
            
            for vp in "${valid_pathogens[@]}"; do
                if [[ "$p_upper" == "$vp" ]]; then
                    valid=true
                    # Add to normalized list
                    if [[ -n "$normalized_pathogens" ]]; then
                        normalized_pathogens="${normalized_pathogens},${vp}"
                    else
                        normalized_pathogens="${vp}"
                    fi
                    break
                fi
            done
            
            if [[ "$valid" == false ]]; then
                echo -e "${RED}Error: '$p' is not a recognized pathogen.${NC}"
                echo -e "${YELLOW}Valid pathogens are: ${valid_pathogens[*]}${NC}"
                invalid_found=true
            fi
        done

        if [[ "$invalid_found" == true ]]; then
            echo ""
            echo "Please enter valid pathogen names."
            # Loop back to pathogen entry
            continue
        else
            # Update pathogens with normalized values
            pathogens="$normalized_pathogens"
        fi
    else
        # More lenient validation when using preprocessed data
        # Just check that something was entered
        if [[ -z "$pathogens" ]] || [[ -z "${pathogens// }" ]]; then
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
while true; do
    read -p "Use custom serotype configuration? (y/n) [n]: " use_serotype_config
    use_serotype_config=${use_serotype_config:-n}
    
    if [[ "$use_serotype_config" =~ ^[YyNn]$ ]]; then
        break
    else
        echo -e "${RED}Invalid input: '$use_serotype_config'. Please enter y or n.${NC}"
    fi
done

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
while true; do
    read -p "Use custom catchment configuration? (y/n) [n]: " use_catchment_config
    use_catchment_config=${use_catchment_config:-n}
    
    if [[ "$use_catchment_config" =~ ^[YyNn]$ ]]; then
        break
    else
        echo -e "${RED}Invalid input: '$use_catchment_config'. Please enter y or n.${NC}"
    fi
done

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
    while true; do
        read -p "Matching sensitivity [MEDIUM]: " matching_sensitivity
        matching_sensitivity=${matching_sensitivity:-MEDIUM}
        # Convert to uppercase for comparison
        matching_sensitivity_upper="${matching_sensitivity^^}"
        
        if [[ "$matching_sensitivity_upper" =~ ^(STRICT|MEDIUM|RELAXED)$ ]]; then
            matching_sensitivity="$matching_sensitivity_upper"
            break
        else
            echo -e "${RED}Invalid input: '$matching_sensitivity'. Please enter STRICT, MEDIUM, or RELAXED.${NC}"
        fi
    done
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
if [[ "$flag" != "resume" ]] && [[ "$pathogens" != "AUTO_DISCOVER" ]]; then
    # Only ask about grouping if we have specific pathogens selected
    # Skip if using AUTO_DISCOVER since we don't know what pathogens exist yet
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

# State selection
selected_states=""
if [[ "$flag" != "resume" ]]; then
    echo ""
    echo -e "${BLUE}======== State Selection ========${NC}"
    
    # Try to get metadata if using preprocessed data
    metadata_available=false
    if [[ "$use_preprocessed" == true ]] && [[ -f "$preprocessed_file" ]]; then
        metadata_dir=$(dirname "$preprocessed_file")
        if [[ -f "${metadata_dir}/metadata_states.csv" ]]; then
            metadata_available=true
        fi
    fi
    
    echo "Select states to analyze:"
    echo "1) ALL states (default)"
    echo "2) Select specific states"
    read -p "Enter selection [1]: " state_mode
    state_mode=${state_mode:-1}
    
    if [[ "$state_mode" == "2" ]]; then
        if [[ "$metadata_available" == true ]]; then
            echo ""
            echo "Available states in data:"
            echo ""
            
            # Read state data and display
            i=1
            declare -a state_array
            while IFS=',' read -r state first_year last_year total_cases; do
                state_name=$(get_state_name "$state")
                printf "%2d. %-2s - %-15s (%d-%d, %'d cases)\n" $i "$state" "$state_name" "$first_year" "$last_year" "$total_cases"
                state_array[$i]=$state
                ((i++))
            done < <(extract_states_from_data "$metadata_dir")
            
            echo ""
            echo "Enter state numbers separated by spaces (e.g., 1 3 5), or 'all' for all states:"
            read -p "Selection: " state_selection
            
            if [[ "$state_selection" == "all" ]]; then
                selected_states=""
            else
                # Convert numbers to state codes
                state_list=""
                for num in $state_selection; do
                    if [[ -n "${state_array[$num]}" ]]; then
                        if [[ -n "$state_list" ]]; then
                            state_list="${state_list},${state_array[$num]}"
                        else
                            state_list="${state_array[$num]}"
                        fi
                    fi
                done
                selected_states=$state_list
                echo -e "${GREEN}Selected states: $selected_states${NC}"
            fi
        else
            # Manual entry without metadata
            echo ""
            echo "Enter state codes separated by commas (e.g., CA,NY,GA):"
            echo "Available states: CA, CO, CT, GA, MD, MN, NM, NY, OR, TN"
            read -p "States: " selected_states
            # Remove spaces
            selected_states=$(echo "$selected_states" | tr -d ' ')
        fi
    fi
fi

# CIDT selection
selected_cidt=""
if [[ "$flag" != "resume" ]]; then
    echo ""
    echo -e "${BLUE}======== Diagnostic Method Selection ========${NC}"
    
    # Try to get metadata if using preprocessed data
    metadata_available=false
    if [[ "$use_preprocessed" == true ]] && [[ -f "$preprocessed_file" ]]; then
        metadata_dir=$(dirname "$preprocessed_file")
        if [[ -f "${metadata_dir}/metadata_cidt.csv" ]]; then
            metadata_available=true
        fi
    fi
    
    echo "Select diagnostic methods:"
    echo "1) ALL methods (default)"
    
    if [[ "$metadata_available" == true ]]; then
        # Read CIDT data to show options
        cidt_data=$(extract_cidt_from_data "$metadata_dir")
        if echo "$cidt_data" | grep -q "CX+"; then
            cx_info=$(echo "$cidt_data" | grep "CX+" | awk -F',' '{printf "%d cases (%.1f%%)", $2, $3}')
            echo "2) Culture only (CX+) - $cx_info"
        fi
        if echo "$cidt_data" | grep -q "CIDT+"; then
            cidt_info=$(echo "$cidt_data" | grep "CIDT+" | awk -F',' '{printf "%d cases (%.1f%%)", $2, $3}')
            echo "3) CIDT only (CIDT+) - $cidt_info"
        fi
        if echo "$cidt_data" | grep -q "PARASITIC"; then
            para_info=$(echo "$cidt_data" | grep "PARASITIC" | awk -F',' '{printf "%d cases (%.1f%%)", $2, $3}')
            echo "4) Parasitic only - $para_info"
        fi
        echo "5) Culture + CIDT (CX+,CIDT+)"
        echo "6) Custom selection"
    else
        echo "2) Culture only (CX+)"
        echo "3) CIDT only (CIDT+)"
        echo "4) Parasitic only (PARASITIC)"
        echo "5) Culture + CIDT (CX+,CIDT+)"
        echo "6) Custom selection"
    fi
    
    read -p "Enter selection [1]: " cidt_mode
    cidt_mode=${cidt_mode:-1}
    
    case $cidt_mode in
        1) selected_cidt="" ;;  # Use default (all)
        2) selected_cidt="CX+" ;;
        3) selected_cidt="CIDT+" ;;
        4) selected_cidt="PARASITIC" ;;
        5) selected_cidt="CX+,CIDT+" ;;
        6) 
            echo "Enter diagnostic methods separated by commas (e.g., CX+,CIDT+):"
            read -p "Methods: " selected_cidt
            selected_cidt=$(echo "$selected_cidt" | tr -d ' ')
            ;;
    esac
    
    if [[ -n "$selected_cidt" ]]; then
        echo -e "${GREEN}Selected diagnostic methods: $selected_cidt${NC}"
    fi
fi

# Travel selection
selected_travel=""
if [[ "$flag" != "resume" ]]; then
    echo ""
    echo -e "${BLUE}======== Travel Status Selection ========${NC}"
    
    # Try to get metadata if using preprocessed data
    metadata_available=false
    if [[ "$use_preprocessed" == true ]] && [[ -f "$preprocessed_file" ]]; then
        metadata_dir=$(dirname "$preprocessed_file")
        if [[ -f "${metadata_dir}/metadata_travel.csv" ]]; then
            metadata_available=true
        fi
    fi
    
    echo "Select travel statuses:"
    echo "1) ALL statuses (default)"
    
    if [[ "$metadata_available" == true ]]; then
        # Read travel data to show percentages
        travel_data=$(extract_travel_from_data "$metadata_dir")
        if echo "$travel_data" | grep -q "NO"; then
            no_info=$(echo "$travel_data" | grep "NO" | awk -F',' '{printf "%d cases (%.1f%%)", $2, $3}')
            echo "2) Domestic only (NO) - $no_info"
        fi
        if echo "$travel_data" | grep -q "YES"; then
            yes_info=$(echo "$travel_data" | grep "YES" | awk -F',' '{printf "%d cases (%.1f%%)", $2, $3}')
            echo "3) Travel-related only (YES) - $yes_info"
        fi
        echo "4) Domestic + Unknown (NO,UNKNOWN)"
        echo "5) Travel + Unknown (YES,UNKNOWN)"
        echo "6) Custom selection"
    else
        echo "2) Domestic only (NO)"
        echo "3) Travel-related only (YES)"
        echo "4) Domestic + Unknown (NO,UNKNOWN)"
        echo "5) Travel + Unknown (YES,UNKNOWN)"
        echo "6) Custom selection"
    fi
    
    read -p "Enter selection [1]: " travel_mode
    travel_mode=${travel_mode:-1}
    
    case $travel_mode in
        1) selected_travel="" ;;  # Use default (all)
        2) selected_travel="NO" ;;
        3) selected_travel="YES" ;;
        4) selected_travel="NO,UNKNOWN" ;;
        5) selected_travel="YES,UNKNOWN" ;;
        6) 
            echo "Enter travel statuses separated by commas (NO,YES,UNKNOWN):"
            read -p "Statuses: " selected_travel
            selected_travel=$(echo "$selected_travel" | tr -d ' ')
            ;;
    esac
    
    if [[ -n "$selected_travel" ]]; then
        echo -e "${GREEN}Selected travel statuses: $selected_travel${NC}"
    fi
fi

# Build the base command
cmd="nextflow run main.nf -profile singularity -entry SPLINE \
  --mmwrFile \"$dataDir/mmwr9623.sas7bdat\" \
  --censusFileB \"$dataDir/cen9623.sas7bdat\" \
  --censusFileP \"$dataDir/cen9623_para.sas7bdat\" \
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
  --mmwrFile \"$dataDir/mmwr9623.sas7bdat\" \
  --censusFileB \"$dataDir/cen9623.sas7bdat\" \
  --censusFileP \"$dataDir/cen9623_para.sas7bdat\" \
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
# Add pathogen grouping if specified and not empty
if [[ -n "$pathogen_grouping" ]] && [[ -n "${pathogen_grouping// }" ]]; then
    cmd="$cmd --pathogen_grouping \"$pathogen_grouping\""
fi
# Add state filter if specified
if [[ -n "$selected_states" ]]; then
    cmd="$cmd --states \"$selected_states\""
fi
# Add CIDT filter if specified
if [[ -n "$selected_cidt" ]]; then
    cmd="$cmd --cidt \"$selected_cidt\""
fi
# Add travel filter if specified
if [[ -n "$selected_travel" ]]; then
    cmd="$cmd --travel \"$selected_travel\""
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
if [[ "$pathogens" == "AUTO_DISCOVER" ]]; then
    echo -e "Pathogens: ${GREEN}All pathogens found in data (auto-discovery)${NC}"
else
    echo -e "Pathogens: ${GREEN}$pathogens${NC}"
fi
# Show selected filters
if [[ -n "$selected_states" ]]; then
    echo -e "States: ${GREEN}$selected_states${NC}"
else
    echo -e "States: ${GREEN}ALL states${NC}"
fi
if [[ -n "$selected_cidt" ]]; then
    echo -e "Diagnostic methods: ${GREEN}$selected_cidt${NC}"
else
    echo -e "Diagnostic methods: ${GREEN}ALL methods (CIDT+,CX+,PARASITIC)${NC}"
fi
if [[ -n "$selected_travel" ]]; then
    echo -e "Travel status: ${GREEN}$selected_travel${NC}"
else
    echo -e "Travel status: ${GREEN}ALL statuses (NO,UNKNOWN,YES)${NC}"
fi

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

while true; do
    read -p "Proceed with analysis? (y/n) [y]: " proceed
    proceed=${proceed:-y}
    
    if [[ "$proceed" =~ ^[YyNn]$ ]]; then
        break
    else
        echo -e "${RED}Invalid input: '$proceed'. Please enter y or n.${NC}"
    fi
done

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
