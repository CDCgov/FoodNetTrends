#!/bin/bash
# =========================================================================
# FoodNet Trends Dashboard Generator - Standalone Script
# =========================================================================
#
# Purpose:
#   Generate the production dashboard outside of the main workflow
#   as a separate and reliable step.
#
# Usage:
#   ./generate_dashboard_only.sh [project_id]
#
# =========================================================================

# Set strict error handling
set -euo pipefail

# Get project ID from command line or use the latest directory in output/
PROJECT_ID=${1:-$(ls -t output/ | grep -v "pipeline_info" | head -n 1)}

# Check if project ID is valid
if [ ! -d "output/$PROJECT_ID" ]; then
  echo "ERROR: Project directory not found: output/$PROJECT_ID"
  echo "Please provide a valid project ID or ensure output directory exists."
  exit 1
fi

echo "==========================================================="
echo "  FoodNet Trends Dashboard Generator"
echo "  Project ID: $PROJECT_ID"
echo "  Output Dir: output/$PROJECT_ID"
echo "==========================================================="

# Ensure output directory exists
mkdir -p "output/$PROJECT_ID"

# Set up environment for R script
export R_MAX_VSIZE=12G
export R_GC_MEM_GROW=0

# List of possible R installation paths
R_PATHS=(
  "Rscript"
  "/usr/bin/Rscript"
  "/usr/local/bin/Rscript"
  "/opt/R/bin/Rscript"
)

# Find a working R
R_CMD=""
for path in "${R_PATHS[@]}"; do
  if command -v "$path" &> /dev/null; then
    R_CMD="$path"
    break
  fi
done

if [ -z "$R_CMD" ]; then
  echo "ERROR: Could not find R installation."
  exit 1
fi

echo "Using R at: $R_CMD"

# Verify we have result files
echo "Checking for result files..."
RESULT_FILES=$(find "output/$PROJECT_ID" -name "*_IRCatch.csv" | wc -l)
echo "Found $RESULT_FILES incidence rate files."

if [ "$RESULT_FILES" -eq 0 ]; then
  echo "WARNING: No result files found - dashboard may be empty."
else
  echo "Result files available for dashboard generation."
fi

# Check for the original generate_dashboard.R script
DASHBOARD_SCRIPT="bin/generate_dashboard.R"
if [ ! -f "$DASHBOARD_SCRIPT" ]; then
  echo "ERROR: Dashboard generation script not found: $DASHBOARD_SCRIPT"
  exit 1
fi

# Create a dashboard directory
DASHBOARD_DIR="output/$PROJECT_ID"
echo "Dashboard will be generated in: $DASHBOARD_DIR"

# Create a dashboard file path
DASHBOARD_FILE="$DASHBOARD_DIR/${PROJECT_ID}_dashboard.html"
echo "Dashboard file will be: $DASHBOARD_FILE"

# Find the dashboard template
TEMPLATE_FILE="assets/dashboard_template.html"
if [ ! -f "$TEMPLATE_FILE" ]; then
  echo "WARNING: Dashboard template not found: $TEMPLATE_FILE"
  echo "Will use default template."
  TEMPLATE_FILE=""
fi

# Execute the dashboard generation with increased resources
echo "Generating production dashboard..."

# Run with detailed error logging
"$R_CMD" \
  "$DASHBOARD_SCRIPT" \
  --outDir="$DASHBOARD_DIR" \
  --resultDir="$DASHBOARD_DIR" \
  --outputFile="$DASHBOARD_FILE" \
  --title="FoodNet Trends Analysis: $PROJECT_ID" \
  --templateFile="$TEMPLATE_FILE"

# Check if dashboard was created
if [ -f "$DASHBOARD_FILE" ]; then
  echo "SUCCESS: Dashboard generated successfully at:"
  echo "$DASHBOARD_FILE"
  echo ""
  echo "Dashboard size: $(du -h "$DASHBOARD_FILE" | cut -f1)"
  echo "Dashboard timestamp: $(date -r "$DASHBOARD_FILE")"
else
  echo "ERROR: Dashboard generation failed - file was not created."
  exit 1
fi

echo ""
echo "Dashboard generation completed at: $(date)"
echo "==========================================================="