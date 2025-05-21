#!/bin/bash
# =========================================================================
# FoodNet Trends Analysis - Run Without Dashboard
# =========================================================================
#
# This script runs the FoodNet Trends pipeline with dashboard generation
# disabled, which avoids the dashboard-related errors that can occur.
#
# Usage:
#   ./run_without_dashboard.sh [additional arguments]
#
# =========================================================================

# Capture all arguments provided to this script
EXTRA_ARGS="$@"

# Get the command line from the original script, but disable dashboard
COMMAND=$(./run_workflow.sh --get-command-only $EXTRA_ARGS)

# Modify the command to disable dashboard generation
FINAL_COMMAND="${COMMAND} --enable_dashboard false"

echo "=========================================================="
echo "Running FoodNet Trends Analysis WITHOUT Dashboard Generation"
echo "=========================================================="
echo "Command: $FINAL_COMMAND"
echo

# Execute the command
eval $FINAL_COMMAND

# Check if successful
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully!"
else
    echo "Analysis failed. Check logs for details."
    exit 1
fi