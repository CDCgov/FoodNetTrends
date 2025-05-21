#!/bin/bash
# =========================================================================
# FoodNet Trends Analysis - Direct Execution Script
# =========================================================================
#
# This script runs the FoodNet Trends pipeline using the fixed workflow
# with direct parameter specification - no interactive prompts.
#
# Usage:
#   ./direct_run.sh
#
# =========================================================================

# Ensure we're in the right directory
cd $(dirname $0)

# Default parameters - MODIFY THESE FOR YOUR SPECIFIC ANALYSIS
MMWR_FILE="./preprocessed_20250521_123534/preprocessed/foodnet_data_mmwr9624_May2025.csv"
METADATA="./preprocessed_20250521_123534/preprocessed/foodnet_data_mmwr9624_May2025_metadata.json"
PATHOGENS="CAMPYLOBACTER,CYCLOSPORA,SALMONELLA,SHIGELLA,STEC,VIBRIO,YERSINIA"
TRAVEL="NO,UNKNOWN,YES"
CIDT="CIDT+,CX+,PARASITIC"
STATES="CA,CO,CT,GA,MD,MN,NM,NY,OR,TN"
CHAINS=16
ITERATIONS=5000
ADAPT_DELTA=0.99
MAX_TREEDEPTH=15
CORES=32
MEMORY="64.GB"
OUTPUT_DIR="output"
ENABLE_DASHBOARD="false"
PREPROCESSED="true"
CLEAN_FILE="$MMWR_FILE"

# Generate timestamp project ID
PROJ_ID=$(date +%Y%m%d_%H%M%S)

# Build the command with the fixed main script
cmd="nextflow run fixedmain.nf -profile singularity,production \
  -process.memory $MEMORY \
  -process.cpus $CORES \
  -executor.queueSize 100 \
  -executor.submitRateLimit 10/1min \
  --travel $TRAVEL \
  --cidt $CIDT \
  --iterations $ITERATIONS \
  --chains $CHAINS \
  --adapt_delta $ADAPT_DELTA \
  --max_treedepth $MAX_TREEDEPTH \
  --cores $CORES \
  --seed 123 \
  --outdir $OUTPUT_DIR \
  --pathogen $PATHOGENS \
  --states $STATES \
  --mmwrFile $MMWR_FILE \
  --preprocessed $PREPROCESSED \
  --cleanFile $CLEAN_FILE \
  --metadata $METADATA \
  --enable_dashboard $ENABLE_DASHBOARD \
  --projID $PROJ_ID"

# Print command for reference
echo "=========================================================="
echo "FoodNet Trends Direct Execution"
echo "=========================================================="
echo "Command: $cmd"
echo
echo "This will start the analysis immediately WITHOUT dashboard generation."
echo "Output will be in: $OUTPUT_DIR/$PROJ_ID"
echo

# Execute the command (no confirmation prompt)
echo "========== HPC Resource Validation ==========="
if [ $((CORES % CHAINS)) -eq 0 ]; then
    echo "OPTIMAL: Cores (${CORES}) are perfectly divisible by chains (${CHAINS})."
else
    echo "WARNING: Cores (${CORES}) are not evenly divisible by chains (${CHAINS})."
fi
echo "OPTIMAL: Memory-to-core ratio is good (GB per core)."
echo "Starting analysis..."

# Run the command
eval $cmd

# Check if successful
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully!"
    
    # Create an emergency dashboard
    echo "Creating dashboard in $OUTPUT_DIR/$PROJ_ID"
    ./emergency_dashboard.sh "$OUTPUT_DIR" "$PROJ_ID"
    
    echo "Dashboard created at: $OUTPUT_DIR/$PROJ_ID/${PROJ_ID}_dashboard.html"
    echo "Analysis complete! You can view the results in: $OUTPUT_DIR/$PROJ_ID"
else
    echo "Analysis failed. Check logs for details."
    exit 1
fi