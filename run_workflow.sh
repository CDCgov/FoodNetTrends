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

echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}   FoodNet Trends Analysis Pipeline      ${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

# Ask for run mode
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

# Ask about pathogens for all modes
echo ""
echo -e "Pathogen selection:"
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

    # Validate pathogens
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

# Build the base command - USING EXACT SAME FORMAT AS ORIGINAL
if [[ "$flag" == "resume" ]]; then
    cmd="nextflow run main.nf -profile singularity -resume -entry SPLINE \
  --mmwrFile \"$dataDir/mmwr9623_Jan2024.sas7bdat\" \
  --censusFileB \"$dataDir/cen9623.sas7bdat\" \
  --censusFileP \"$dataDir/cen9623_para.sas7bdat\" \
  --iterations $iterations \
  --chains $chains \
  --adapt_delta $adapt_delta \
  --max_treedepth $max_treedepth \
  --seed 123 \
  --outdir \"$outDir\" \
  --pathogen \"$pathogens\""
else
    cmd="nextflow run main.nf -profile singularity -entry SPLINE \
  --mmwrFile \"$dataDir/mmwr9623_Jan2024.sas7bdat\" \
  --censusFileB \"$dataDir/cen9623.sas7bdat\" \
  --censusFileP \"$dataDir/cen9623_para.sas7bdat\" \
  --iterations $iterations \
  --chains $chains \
  --adapt_delta $adapt_delta \
  --max_treedepth $max_treedepth \
  --seed 123 \
  --outdir \"$outDir\" \
  --pathogen \"$pathogens\""
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
echo -e "Mode: ${GREEN}$([ "$flag" == "test" ] && echo "Test run" || [ "$flag" == "full" ] && echo "Full analysis" || echo "Resume previous run")${NC}"
echo -e "Pathogens: ${GREEN}$pathogens${NC}"

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
    eval $final_cmd
    
    if [[ $background == true ]]; then
        echo -e "${GREEN}Process started in background. Check status with:${NC}"
        echo -e "${YELLOW}tail -f foodnet_run_${timestamp}.log${NC}"
    fi
else
    echo -e "${RED}Analysis cancelled.${NC}"
fi
