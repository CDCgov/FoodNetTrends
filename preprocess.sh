#!/bin/bash

# FoodNetTrends Preprocessing Script - Simplified Output Version
# This script performs data cleaning and discovers available pathogens and states

# set up paths, files
DEFAULT_DATA_DIR="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/"
DEFAULT_MMWR_FILE="${DEFAULT_DATA_DIR}/mmwr9623_Jan2024.sas7bdat"

# set up modules
module purge
module load nextflow/24.10.4
module load singularity/4.1.4
module load java/17.0.6

# Define color codes for better user experience
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=========================================${NC}"
echo -e "${BLUE}   FoodNet Trends Data Preprocessing     ${NC}"
echo -e "${BLUE}=========================================${NC}"
echo ""

# Ask for input file first
echo -e "${BLUE}======== Input File ========${NC}"
read -p "MMWR data file [${DEFAULT_MMWR_FILE}]: " mmwrFile
mmwrFile=${mmwrFile:-$DEFAULT_MMWR_FILE}

# Validate that file exists
if [ ! -f "$mmwrFile" ]; then
  echo -e "${RED}Error: MMWR file does not exist: $mmwrFile${NC}"
  echo -e "${RED}Exiting.${NC}"
  exit 1
fi

# Extract the base name of the file for the default output directory
FILENAME=$(basename "$mmwrFile")
BASENAME="${FILENAME%.*}"  # Remove extension
# Default output directory based on the file name
outDir="preprocessed_${BASENAME}"

# Ask for output directory
echo ""
echo -e "Where should preprocessed files be saved?"
read -p "Output directory [${outDir}]: " user_outdir
outDir=${user_outdir:-$outDir}

# Run Nextflow preprocessing command
echo -e "${BLUE}Running preprocessing through Nextflow...${NC}"

# Create output directory
mkdir -p "$outDir"

# Build the Nextflow command using the existing PREPROCESS process
preprocess_cmd="nextflow run main.nf -profile singularity -entry PREPROCESS_ONLY \
  --mmwrFile \"$mmwrFile\" \
  --outdir \"$outDir\""

# Display the command
echo -e "Command to run:"
echo -e "${YELLOW}$preprocess_cmd${NC}"
echo ""

# Run the command
read -p "Proceed with preprocessing? (y/n) [y]: " proceed
proceed=${proceed:-y}

if [[ "$proceed" =~ ^[Yy]$ ]]; then
  echo -e "${GREEN}Starting preprocessing...${NC}"
  
  # Run the command and capture success/failure
  if eval $preprocess_cmd; then
    # Success - Nextflow completed without errors
    echo -e "${GREEN}Nextflow preprocessing completed successfully!${NC}"
    
    # Wait a moment for files to be fully written
    sleep 2
    
    # Check for files in the root output directory
    CLEAN_CSV="$outDir/clean_mmwr.csv"
    METADATA_FILE="$outDir/clean_mmwr_metadata.json"
    
    # If files aren't in the root, run the cleanup script
    if [ ! -f "$CLEAN_CSV" ]; then
      echo -e "${YELLOW}Files may be in nested directories. Looking for them...${NC}"
      
      # Find the clean_mmwr.csv file in the nested directory structure
      FOUND_CSV=$(find "$outDir" -name "clean_mmwr.csv" -type f | head -n 1)
      
      if [ -n "$FOUND_CSV" ]; then
        echo -e "${GREEN}Found clean_mmwr.csv at: $FOUND_CSV${NC}"
        cp "$FOUND_CSV" "$CLEAN_CSV"
        echo -e "${GREEN}Copied to: $CLEAN_CSV${NC}"
        
        # Find metadata file if it exists
        FOUND_METADATA=$(find "$(dirname "$FOUND_CSV")" -name "*metadata.json" -type f | head -n 1)
        
        if [ -n "$FOUND_METADATA" ]; then
          echo -e "${GREEN}Found metadata file at: $FOUND_METADATA${NC}"
          cp "$FOUND_METADATA" "$METADATA_FILE"
          echo -e "${GREEN}Copied to: $METADATA_FILE${NC}"
        fi
      else
        echo -e "${RED}Could not find clean_mmwr.csv in $outDir or subdirectories.${NC}"
        echo -e "${RED}Preprocessing may have failed.${NC}"
        exit 1
      fi
    fi
    
    # At this point, we should have the clean_mmwr.csv in the root directory
    if [ -f "$CLEAN_CSV" ]; then
      echo -e "${GREEN}Preprocessing completed successfully!${NC}"
      echo -e "${GREEN}Cleaned data: $CLEAN_CSV${NC}"
      
      # If we have metadata, display it
      if [ -f "$METADATA_FILE" ]; then
        echo -e "${GREEN}Metadata: $METADATA_FILE${NC}"
        
        # Display available pathogens and states
        echo ""
        echo -e "${BLUE}======== Discovered Data ========${NC}"
        
        if command -v jq &> /dev/null; then
          # Use jq to parse JSON
          echo -e "${BLUE}Available pathogens:${NC}"
          jq -r '.pathogens | keys[]' "$METADATA_FILE" | sort | sed 's/^/- /'
          
          pathogen_count=$(jq -r '.pathogens | keys | length' "$METADATA_FILE")
          echo -e "${GREEN}Total: $pathogen_count pathogens${NC}"
          
          echo ""
          echo -e "${BLUE}Available states:${NC}"
          jq -r '.states | keys[]' "$METADATA_FILE" | sort | sed 's/^/- /'
          
          state_count=$(jq -r '.states | keys | length' "$METADATA_FILE")
          echo -e "${GREEN}Total: $state_count states${NC}"
          
          # Check for Salmonella serotypes
          if jq -e '.salmonella_serotypes' "$METADATA_FILE" > /dev/null; then
            serotype_count=$(jq -r '.salmonella_serotypes | keys | length' "$METADATA_FILE")
            
            if [ "$serotype_count" -gt 0 ]; then
              echo ""
              echo -e "${BLUE}Found $serotype_count Salmonella serotypes${NC}"
              echo -e "${BLUE}Top 10 Salmonella serotypes:${NC}"
              # Fixed command - just sort alphabetically
              jq -r '.salmonella_serotypes | keys | sort | .[0:10]' "$METADATA_FILE" | sed 's/^/- /'
            fi
          fi
        else
          echo -e "${YELLOW}jq not found. Install jq for better JSON parsing.${NC}"
          echo -e "${YELLOW}Metadata file saved to: $METADATA_FILE${NC}"
        fi
      else
        # No metadata file - generate one
        echo -e "${YELLOW}No metadata file found. Generating metadata...${NC}"
        
        # Use Rscript directly to generate metadata
        Rscript -e "
        library(jsonlite)
        library(dplyr)
        
        # Read CSV file
        data <- read.csv('$CLEAN_CSV')
        
        # Extract pathogens
        pathogens <- sort(unique(data\$pathogen))
        pathogen_counts <- table(data\$pathogen)
        pathogen_list <- as.list(pathogen_counts)
        
        # Extract states
        states <- sort(unique(data\$state))
        state_counts <- table(data\$state)
        state_list <- as.list(state_counts)
        
        # Extract Salmonella serotypes if possible
        salmonella_serotypes <- list()
        if ('serotypesummary' %in% names(data)) {
          sal_data <- data[data\$pathogen == 'SALMONELLA', ]
          sero_counts <- table(sal_data\$serotypesummary)
          salmonella_serotypes <- as.list(sero_counts)
        }
        
        # Create metadata object
        metadata <- list(
          pathogens = pathogen_list,
          states = state_list,
          salmonella_serotypes = salmonella_serotypes,
          generated_date = as.character(Sys.time()),
          source_file = basename('$CLEAN_CSV')
        )
        
        # Write to JSON file
        write_json(metadata, '$METADATA_FILE', pretty = TRUE)
        
        cat('Metadata saved to: $METADATA_FILE\n')
        " 2>/dev/null
        
        # Check if metadata was generated
        if [ -f "$METADATA_FILE" ]; then
          echo -e "${GREEN}Generated metadata: $METADATA_FILE${NC}"
          
          # Display metadata summary
          echo ""
          echo -e "${BLUE}======== Discovered Data ========${NC}"
          
          if command -v jq &> /dev/null; then
            echo -e "${BLUE}Available pathogens:${NC}"
            jq -r '.pathogens | keys[]' "$METADATA_FILE" | sort | sed 's/^/- /'
            
            pathogen_count=$(jq -r '.pathogens | keys | length' "$METADATA_FILE")
            echo -e "${GREEN}Total: $pathogen_count pathogens${NC}"
            
            echo ""
            echo -e "${BLUE}Available states:${NC}"
            jq -r '.states | keys[]' "$METADATA_FILE" | sort | sed 's/^/- /'
            
            state_count=$(jq -r '.states | keys | length' "$METADATA_FILE")
            echo -e "${GREEN}Total: $state_count states${NC}"
          else
            echo -e "${YELLOW}jq not found. Metadata saved but cannot display summary.${NC}"
          fi
        else
          echo -e "${YELLOW}Could not generate metadata, but preprocessing was successful.${NC}"
        fi
      fi
      
      # Cleanup nested directories (optional)
      echo ""
      read -p "Clean up nested directories? (y/n) [y]: " cleanup
      cleanup=${cleanup:-y}
      
      if [[ "$cleanup" =~ ^[Yy]$ ]]; then
        echo -e "${BLUE}Cleaning up nested directories...${NC}"
        
        # Move pipeline_info if it exists but in a nested location
        PIPELINE_INFO=$(find "$outDir" -name "pipeline_info" -type d | head -n 1)
        if [ -n "$PIPELINE_INFO" ] && [ "$PIPELINE_INFO" != "$outDir/pipeline_info" ]; then
          mkdir -p "$outDir/pipeline_info"
          cp -r "$PIPELINE_INFO"/* "$outDir/pipeline_info/"
          echo -e "${GREEN}Preserved pipeline_info in root directory.${NC}"
        fi
        
        # Remove all subdirectories except pipeline_info in the root
        for dir in $(find "$outDir" -mindepth 1 -type d | grep -v "^$outDir/pipeline_info$"); do
          echo -e "${YELLOW}Removing: $dir${NC}"
          rm -rf "$dir"
        done
        
        echo -e "${GREEN}Directory cleanup complete.${NC}"
        echo -e "${GREEN}Final directory structure:${NC}"
        ls -la "$outDir"
      fi
      
      echo ""
      echo -e "${GREEN}You can now run the analysis with the preprocessed data:${NC}"
      echo -e "${YELLOW}./run_workflow.sh${NC}"
      echo -e "${GREEN}When prompted, specify this preprocessing directory:${NC}"
      echo -e "${YELLOW}$outDir${NC}"
    else
      echo -e "${RED}Preprocessing failed or did not generate expected files.${NC}"
      echo -e "${RED}Check .nextflow.log for details.${NC}"
    fi
  else
    # Command failed with non-zero exit code
    echo -e "${RED}Preprocessing failed with errors.${NC}"
    echo -e "${RED}Check .nextflow.log for details.${NC}"
  fi
else
  echo -e "${RED}Preprocessing cancelled.${NC}"
fi
