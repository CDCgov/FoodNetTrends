#!/bin/bash
# =========================================================================
# FoodNet Trends Pipeline - Path Verification Script
# =========================================================================
#
# This script verifies that all required files and directories exist
# for the FoodNet Trends pipeline, specifically for dashboard generation.
#
# Usage:
#   ./verify_paths.sh
#
# =========================================================================

# Strict error handling
set -e

# Text formatting
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RESET='\033[0m'

echo "=========================================================="
echo "FoodNet Trends Dashboard Path Verification"
echo "=========================================================="

# Get top-level project directory
PROJECT_DIR=$(pwd)
echo "Project directory: $PROJECT_DIR"

# Check if we're in the right directory
if [ ! -f "main.nf" ]; then
    echo -e "${RED}ERROR: Not in the FoodNet Trends project directory (main.nf not found)${RESET}"
    exit 1
fi

# Check workflow files
echo -e "\nChecking workflow files..."

# Check main workflow files
check_file() {
    if [ -f "$1" ]; then
        echo -e "  ${GREEN}✓${RESET} $1 exists"
        return 0
    else
        echo -e "  ${RED}✗${RESET} $1 MISSING"
        return 1
    fi
}

check_dir() {
    if [ -d "$1" ]; then
        echo -e "  ${GREEN}✓${RESET} $1 directory exists"
        return 0
    else
        echo -e "  ${RED}✗${RESET} $1 directory MISSING"
        return 1
    fi
}

ERRORS=0

# Check main workflow files
check_file "main.nf" || ERRORS=$((ERRORS+1))
check_file "workflows/spline_fixed.nf" || ERRORS=$((ERRORS+1))
check_file "workflows/preprocess.nf" || ERRORS=$((ERRORS+1))
check_file "modules/local/generate_dashboard.nf" || ERRORS=$((ERRORS+1))
check_file "modules/local/trendy.nf" || ERRORS=$((ERRORS+1))

# Check dashboard files
echo -e "\nChecking dashboard files..."
check_file "assets/dashboard_template.html" || ERRORS=$((ERRORS+1))
check_file "bin/generate_dashboard.R" || ERRORS=$((ERRORS+1))

# Check execution scripts
echo -e "\nChecking execution scripts..."
check_file "run_workflow.sh" || ERRORS=$((ERRORS+1))
check_file "run_workflow_hpc.sh" || ERRORS=$((ERRORS+1))

# Check alternative scripts
echo -e "\nChecking alternative scripts..."
check_file "run_without_dashboard.sh" || ERRORS=$((ERRORS+1))
check_file "run_without_dashboard_hpc.sh" || ERRORS=$((ERRORS+1))
check_file "emergency_dashboard.sh" || ERRORS=$((ERRORS+1))
check_file "direct_run.sh" || ERRORS=$((ERRORS+1))

# Check for any existing work directories
echo -e "\nChecking for existing work directories..."
check_dir "output" && ls -la output | head -n 5
find . -name "work" -type d | sort | while read -r workdir; do
    if [ -d "$workdir" ]; then
        echo -e "  ${YELLOW}!${RESET} Found work directory: $workdir"
    fi
done

# Check for preprocessed data
echo -e "\nChecking for preprocessed data..."
if [ -d "preprocessed_20250521_123534" ]; then
    echo -e "  ${GREEN}✓${RESET} Preprocessed data directory exists"
    check_file "preprocessed_20250521_123534/preprocessed/foodnet_data_mmwr9624_May2025.csv" || ERRORS=$((ERRORS+1))
    check_file "preprocessed_20250521_123534/preprocessed/foodnet_data_mmwr9624_May2025_metadata.json" || ERRORS=$((ERRORS+1))
else
    echo -e "  ${YELLOW}!${RESET} Preprocessed data directory not found"
    ERRORS=$((ERRORS+1))
fi

# Check NextFlow configuration
echo -e "\nChecking NextFlow configuration..."
check_file "nextflow.config" || ERRORS=$((ERRORS+1))
if grep -q "enable_dashboard.*true" nextflow.config; then
    echo -e "  ${GREEN}✓${RESET} Dashboard is enabled by default in nextflow.config"
else
    echo -e "  ${YELLOW}!${RESET} Dashboard might not be enabled by default in nextflow.config"
fi

# Summary
echo -e "\n=========================================================="
if [ $ERRORS -eq 0 ]; then
    echo -e "${GREEN}All path verifications passed!${RESET}"
    echo -e "The FoodNet Trends pipeline appears to be correctly set up."
    echo -e "\nYou can safely run the pipeline with dashboard generation using:"
    echo -e "  ${GREEN}./run_workflow_hpc.sh${RESET}"
else
    echo -e "${RED}Found $ERRORS issues with paths or files!${RESET}"
    echo -e "Please fix these issues before running the pipeline."
    echo -e "\nYou can try running without dashboard generation using:"
    echo -e "  ${YELLOW}./run_without_dashboard_hpc.sh${RESET}"
fi
echo -e "=========================================================="