#!/bin/bash
# =========================================================================
# FoodNet Trends Dashboard Generator Script - STANDALONE OPTION
# =========================================================================
#
# Purpose:
#   Run just the dashboard generation workflow with Nextflow and container.
#   This is a STANDALONE option that can be used to regenerate the dashboard
#   after the main analysis has completed.
#
# Usage:
#   ./dashboard.sh [project_id]
#
# =========================================================================

# Get project ID from argument or use current timestamp
PROJECT_ID=${1:-$(date +%Y%m%d_%H%M%S)}
OUTPUT_DIR="output"

echo "============================================================"
echo "FoodNet Trends Dashboard Generator"
echo "============================================================"
echo "Project ID:     $PROJECT_ID"
echo "Output Directory: $OUTPUT_DIR/$PROJECT_ID"
echo "Generated at:   $(date)"
echo "============================================================"
echo "This script will generate an ENHANCED dashboard with interactive"
echo "visualizations, trend charts, and data tables for the analysis"
echo "results located in $OUTPUT_DIR/$PROJECT_ID."
echo "============================================================"
echo

# Check if output directory exists and has result files
if [ ! -d "$OUTPUT_DIR/$PROJECT_ID" ]; then
  echo "ERROR: Output directory does not exist: $OUTPUT_DIR/$PROJECT_ID"
  echo "This suggests that the analysis has not been run yet."
  echo "Please run the main analysis pipeline first, then generate the dashboard."
  exit 1
fi

# Check for result files before proceeding
RESULT_FILES=$(find "$OUTPUT_DIR/$PROJECT_ID" -name "*_IRCatch.csv" | wc -l)
if [ "$RESULT_FILES" -eq 0 ]; then
  echo "ERROR: No result files found in $OUTPUT_DIR/$PROJECT_ID"
  echo "This suggests that the analysis has not completed yet or failed."
  echo "Please ensure the main pipeline has finished successfully before running this script."
  exit 1
else
  echo "Found $RESULT_FILES result files in output directory."
  echo "Proceeding with dashboard generation..."
fi

# Load necessary modules if available
if command -v module &> /dev/null; then
    module purge
    module load nextflow/24.10.4
    module load singularity/4.1.4
    module load java/17.0.6
    echo "Loaded required modules for HPC environment"
else
    echo "Module system not detected, using environment PATH"
fi

# Run the dashboard-only workflow
nextflow run main.nf -entry DASHBOARD_GEN \
  -profile singularity,production \
  -process.memory 16.GB \
  -process.cpus 4 \
  -executor.queueSize 100 \
  -executor.submitRateLimit '10/1min' \
  --projID "$PROJECT_ID" \
  --outdir "$OUTPUT_DIR" \
  --enable_dashboard true \
  --dashboard_title "FoodNet Trends Analysis: $PROJECT_ID" \
  --dashboard_theme "modern" \
  --publish_dir_mode copy

# Check if dashboard was created
DASHBOARD_FILE="$OUTPUT_DIR/$PROJECT_ID/${PROJECT_ID}_dashboard.html"
if [ -f "$DASHBOARD_FILE" ]; then
  echo "Dashboard created successfully at: $DASHBOARD_FILE"
  echo "File size: $(du -h "$DASHBOARD_FILE" | cut -f1)"
else
  echo "Warning: Dashboard file not found at expected location: $DASHBOARD_FILE"
  echo "Checking for alternate naming patterns..."
  
  # Look for any dashboard files
  FOUND_DASHBOARDS=$(find "$OUTPUT_DIR/$PROJECT_ID" -name "*dashboard*.html" | wc -l)
  if [ "$FOUND_DASHBOARDS" -gt 0 ]; then
    echo "Found $FOUND_DASHBOARDS dashboard files:"
    find "$OUTPUT_DIR/$PROJECT_ID" -name "*dashboard*.html" -exec ls -lh {} \;
  else
    echo "ERROR: No dashboard files found."
    exit 1
  fi
fi

echo
echo "Dashboard generation completed at: $(date)"
echo "============================================================"