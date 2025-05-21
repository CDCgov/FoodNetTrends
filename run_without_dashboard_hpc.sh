#!/bin/bash
# =========================================================================
# FoodNet Trends Analysis - Run Without Dashboard (HPC Version)
# =========================================================================
#
# This script runs the FoodNet Trends pipeline in HPC mode but with dashboard
# generation disabled, which avoids the dashboard-related errors.
#
# Usage:
#   ./run_without_dashboard_hpc.sh [additional arguments]
#
# =========================================================================

# Capture all arguments provided to this script
EXTRA_ARGS="$@"

# Get the command line from the HPC script, but disable dashboard
COMMAND=$(./run_workflow_hpc.sh --get-command-only $EXTRA_ARGS)

# Modify the command to disable dashboard generation
FINAL_COMMAND="${COMMAND} --enable_dashboard false"

echo "=========================================================="
echo "Running FoodNet Trends Analysis (HPC Mode) WITHOUT Dashboard Generation"
echo "=========================================================="
echo "Command: $FINAL_COMMAND"
echo

# Execute the command
eval $FINAL_COMMAND

# Check if successful
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully!"
    
    # Create an emergency dashboard
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    PROJ_ID=$(echo $FINAL_COMMAND | grep -o "projID:[^ ]*" | cut -d':' -f2)
    if [ -z "$PROJ_ID" ]; then
        PROJ_ID=$TIMESTAMP
    fi
    
    OUTPUT_DIR=$(echo $FINAL_COMMAND | grep -o "outdir:[^ ]*" | cut -d':' -f2)
    if [ -z "$OUTPUT_DIR" ]; then
        OUTPUT_DIR="output"
    fi
    
    echo "Creating emergency dashboard in $OUTPUT_DIR/$PROJ_ID"
    ./emergency_dashboard.sh "$OUTPUT_DIR" "$PROJ_ID"
else
    echo "Analysis failed. Check logs for details."
    exit 1
fi