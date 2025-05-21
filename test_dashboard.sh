#!/usr/bin/env bash
# =========================================================================
# FoodNet Trends Dashboard Generation Test
# =========================================================================
#
# Purpose:
#   Verify dashboard generation process by creating minimal test data
#   and running the dashboard generation process manually.
#
# Usage:
#   ./test_dashboard.sh [output_directory]
#

# Set error handling
set -e

# Set defaults
OUTPUT_DIR=${1:-./test_output}
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
DASHBOARD_NAME="${TIMESTAMP}_dashboard.html"

echo "========================================================="
echo "Testing FoodNet Trends Dashboard Generation"
echo "========================================================="
echo "Output directory: $OUTPUT_DIR"
echo "Dashboard file: $DASHBOARD_NAME"
echo ""

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
  echo "Creating output directory..."
  mkdir -p "$OUTPUT_DIR"
fi

# Create a test directory for our files
TEST_DIR="${OUTPUT_DIR}/test_data"
mkdir -p "$TEST_DIR"

# Create test files - minimal incidence rate data
echo "Creating test incidence rate data..."
cat > "${TEST_DIR}/CAMPYLOBACTER_IRCatch.csv" << EOF
state,year,median_incidence,lower_hdi,upper_hdi,pathogen,travel,culture
CA,2022,12.5,10.2,14.8,CAMPYLOBACTER,NO,CX+
CA,2023,13.1,10.9,15.3,CAMPYLOBACTER,NO,CX+
GA,2022,9.8,8.1,11.5,CAMPYLOBACTER,NO,CX+
GA,2023,10.2,8.5,11.9,CAMPYLOBACTER,NO,CX+
EOF

cat > "${TEST_DIR}/SALMONELLA_IRCatch.csv" << EOF
state,year,median_incidence,lower_hdi,upper_hdi,pathogen,travel,culture
CA,2022,8.3,7.1,9.5,SALMONELLA,NO,CX+
CA,2023,7.9,6.8,9.1,SALMONELLA,NO,CX+
GA,2022,7.5,6.4,8.6,SALMONELLA,NO,CX+
GA,2023,7.2,6.1,8.3,SALMONELLA,NO,CX+
EOF

# Create test summary files
echo "Creating test summary files..."
cat > "${TEST_DIR}/CAMPYLOBACTER_summary.txt" << EOF
===========================================================
FoodNet Trends Analysis: CAMPYLOBACTER
===========================================================
Analysis complete at: $(date)
Iterations: 500
Chains: 2
Travel: NO
CIDT: CX+
===========================================================
EOF

cat > "${TEST_DIR}/SALMONELLA_summary.txt" << EOF
===========================================================
FoodNet Trends Analysis: SALMONELLA
===========================================================
Analysis complete at: $(date)
Iterations: 500
Chains: 2
Travel: NO
CIDT: CX+
===========================================================
EOF

# Create test relative risk files
echo "Creating test relative risk data..."
cat > "${TEST_DIR}/CAMPYLOBACTER_EstIRRCatch_2022-2023.csv" << EOF
state,year,comparison_period,current_incidence,period_incidence,relative_risk,percent_change,pathogen
CA,2023,2022,13.1,12.5,1.048,4.8,CAMPYLOBACTER
GA,2023,2022,10.2,9.8,1.041,4.1,CAMPYLOBACTER
EOF

cat > "${TEST_DIR}/SALMONELLA_EstIRRCatch_2022-2023.csv" << EOF
state,year,comparison_period,current_incidence,period_incidence,relative_risk,percent_change,pathogen
CA,2023,2022,7.9,8.3,0.952,-4.8,SALMONELLA
GA,2023,2022,7.2,7.5,0.96,-4.0,SALMONELLA
EOF

# Create data quality JSON for testing
echo "Creating test data quality metadata..."
cat > "${TEST_DIR}/test_data_quality.json" << EOF
{
  "projectId": "test_dashboard",
  "generationDate": "$(date -Iseconds)",
  "dataQuality": {
    "usesPlaceholderData": false,
    "affectedPathogens": [],
    "warnings": [],
    "dataConsistency": {
      "incidenceRateFiles": 4,
      "dataQualityWarnings": 0
    }
  }
}
EOF

# Copy template file for testing
echo "Setting up dashboard template..."
cp -v "assets/dashboard_template.html" "${TEST_DIR}/"

# Run the dashboard generation script directly
echo "Running dashboard generation script..."
cd "$TEST_DIR"
Rscript "${BIN_DIR:-../bin}/generate_dashboard.R" \
  --outDir="$OUTPUT_DIR" \
  --resultDir="." \
  --outputFile="$DASHBOARD_NAME" \
  --title="Test Dashboard: $TIMESTAMP" \
  --templateFile="dashboard_template.html" \
  --qualityDataPath="test_data_quality.json"

# Check if dashboard was successfully created
if [ -f "$DASHBOARD_NAME" ]; then
  echo "Success! Dashboard generated at: ${TEST_DIR}/$DASHBOARD_NAME"
  echo "File size: $(du -h "$DASHBOARD_NAME" | cut -f1)"
  mv "$DASHBOARD_NAME" "$OUTPUT_DIR/"
  echo "Dashboard moved to: ${OUTPUT_DIR}/$DASHBOARD_NAME"
else
  echo "ERROR: Dashboard generation failed!"
  echo "Check error messages above."
  exit 1
fi

echo ""
echo "========================================================="
echo "Dashboard testing complete"
echo "View your dashboard at: file://${OUTPUT_DIR}/$DASHBOARD_NAME"
echo "========================================================="