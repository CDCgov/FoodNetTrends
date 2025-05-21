#!/bin/bash
# =========================================================================
# FINAL FIX RUN - Simplified Direct Execution
# =========================================================================
#
# This script runs the FoodNet Trends pipeline using the simplified workflow
# to avoid all variable scope and dashboard issues.
#
# Usage:
#   ./FINAL_FIX_RUN.sh
#
# =========================================================================

# Clear screen and show header
clear
echo "============================================================"
echo "FoodNet Trends FINAL FIX - Direct Run with Dashboard"
echo "============================================================"
echo

# Determine paths - use same as previous runs
MMWR_FILE=${1:-"./preprocessed_20250521_123534/preprocessed/foodnet_data_mmwr9624_May2025.csv"}
CENSUS_B=${2:-"/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9624.sas7bdat"}
CENSUS_P=${3:-"/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9624_para.sas7bdat"}

# Verify paths
echo "Verifying paths..."
for path in "$MMWR_FILE" "$CENSUS_B" "$CENSUS_P"; do
  if [ ! -f "$path" ]; then
    echo "ERROR: File not found: $path"
    exit 1
  else
    echo "✓ Found: $path"
  fi
done

# Create timestamp for output
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
PROJ_ID="run_${TIMESTAMP}"
OUTPUT_DIR="output"

# Display configuration
echo
echo "============================================================"
echo "Configuration:"
echo "============================================================"
echo "MMWR File:              $MMWR_FILE"
echo "Bacterial Census File:  $CENSUS_B"
echo "Parasitic Census File:  $CENSUS_P"
echo "Project ID:             $PROJ_ID"
echo "Output Directory:       $OUTPUT_DIR/$PROJ_ID"
echo "Dashboard:              ENABLED"
echo "Dashboard Title:        FoodNet Trends Analysis - $PROJ_ID"
echo 

echo "Ready to run? Press ENTER to continue or CTRL+C to cancel..."
read -r

# Run the workflow with the simplified workflow
nextflow run main.nf \
  -profile singularity,production \
  -process.memory 64.GB \
  -process.cpus 32 \
  -executor.queueSize 100 \
  -executor.submitRateLimit 10/1min \
  --travel NO,UNKNOWN,YES \
  --cidt "CIDT+,CX+,PARASITIC" \
  --iterations 5000 \
  --chains 16 \
  --adapt_delta 0.99 \
  --max_treedepth 15 \
  --cores 32 \
  --seed 123 \
  --outdir "$OUTPUT_DIR" \
  --pathogen CAMPYLOBACTER,CYCLOSPORA,SALMONELLA,SHIGELLA,STEC,VIBRIO,YERSINIA \
  --states CA,CO,CT,GA,MD,MN,NM,NY,OR,TN \
  --mmwrFile "$MMWR_FILE" \
  --censusFileB "$CENSUS_B" \
  --censusFileP "$CENSUS_P" \
  --preprocessed true \
  --cleanFile "$MMWR_FILE" \
  --enable_dashboard true \
  --dashboard_title "FoodNet Trends Analysis - $PROJ_ID" \
  --projID "$PROJ_ID"

EXIT_CODE=$?

# Check if successful
if [ $EXIT_CODE -eq 0 ]; then
  echo
  echo "============================================================"
  echo "SUCCESS: Pipeline completed successfully"
  echo "============================================================"
  echo "You can find the results in: $OUTPUT_DIR/$PROJ_ID"
  echo "Dashboard should be at: $OUTPUT_DIR/$PROJ_ID/${PROJ_ID}_dashboard.html"
  
  # Make a note about the emergency dashboard option
  echo
  echo "If the dashboard is missing or incomplete, you can create an"
  echo "emergency dashboard with:"
  echo "./emergency_dashboard.sh $OUTPUT_DIR $PROJ_ID"
else
  echo
  echo "============================================================"
  echo "ERROR: Pipeline failed with exit code $EXIT_CODE"
  echo "============================================================"
  echo "Creating emergency dashboard..."
  
  # Create emergency dashboard
  ./emergency_dashboard.sh "$OUTPUT_DIR" "$PROJ_ID"
  
  echo "Emergency dashboard created at: $OUTPUT_DIR/$PROJ_ID/${PROJ_ID}_dashboard.html"
fi

echo
echo "FINISHED: $(date)"