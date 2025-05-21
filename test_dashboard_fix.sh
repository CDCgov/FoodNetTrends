#!/bin/bash
# =========================================================================
# FoodNet Trends Dashboard Generation Fix Test
# =========================================================================
#
# Purpose: 
#   Test that dashboard generation produces files with the correct patterns
#   that Nextflow expects, including both projID and timestamp formats.
#
# Usage:
#   ./test_dashboard_fix.sh
#
# =========================================================================

set -e  # Exit on error

echo "==============================================================="
echo "Testing Dashboard Generation Fix"
echo "==============================================================="

# Create test directory
TEST_DIR="./test_dashboard_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"
echo "Created test directory: $TEST_DIR"

# Copy required files
cp bin/generate_dashboard.R "$TEST_DIR/"
cp assets/dashboard_template.html "$TEST_DIR/"

# Create minimal test data
echo "Creating test data..."
mkdir -p "$TEST_DIR/data"

# Create test data quality JSON
cat > "$TEST_DIR/data/test_data_quality.json" << EOF
{
  "projectId": "test_project",
  "generationDate": "$(date -Iseconds)",
  "dataQuality": {
    "usesPlaceholderData": false,
    "affectedPathogens": [],
    "warnings": [],
    "dataConsistency": {
      "incidenceRateFiles": 2,
      "dataQualityWarnings": 0
    }
  }
}
EOF

# Create test incidence rate data
cat > "$TEST_DIR/data/CAMPYLOBACTER_IRCatch.csv" << EOF
state,year,median_incidence,lower_hdi,upper_hdi,pathogen,travel,culture
CA,2022,12.5,10.2,14.8,CAMPYLOBACTER,NO,CX+
CA,2023,13.1,10.9,15.3,CAMPYLOBACTER,NO,CX+
EOF

# Test running the script directly
cd "$TEST_DIR"

echo "==============================================================="
echo "Test 1: Standard Output Filename"
echo "==============================================================="
PROJECT_ID="test_project"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

# Run the R script for standard output
Rscript generate_dashboard.R \
  --outDir="." \
  --resultDir="./data" \
  --outputFile="${PROJECT_ID}_dashboard.html" \
  --title="Test Dashboard" \
  --templateFile="dashboard_template.html" \
  --qualityDataPath="data/test_data_quality.json"

# Check if both patterns exist
if [ -f "${PROJECT_ID}_dashboard.html" ]; then
  echo "✅ SUCCESS: Project ID dashboard created: ${PROJECT_ID}_dashboard.html"
else
  echo "❌ FAILURE: Project ID dashboard not created"
  exit 1
fi

# Check if timestamp format was also created
TIMESTAMP_FILES=(*_dashboard.html)
TIMESTAMP_COUNT="${#TIMESTAMP_FILES[@]}"

if [ "$TIMESTAMP_COUNT" -gt 1 ]; then
  echo "✅ SUCCESS: Multiple dashboard files created (${TIMESTAMP_COUNT})"
  ls -la *_dashboard.html
else
  echo "⚠️ WARNING: Only one dashboard file created"
  ls -la *_dashboard.html 
fi

echo "==============================================================="
echo "Test 2: Timestamp Output Filename"
echo "==============================================================="
rm -f *_dashboard.html
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

# Run the R script with timestamp format
Rscript generate_dashboard.R \
  --outDir="." \
  --resultDir="./data" \
  --outputFile="${TIMESTAMP}_dashboard.html" \
  --title="Test Dashboard (Timestamp)" \
  --templateFile="dashboard_template.html" \
  --qualityDataPath="data/test_data_quality.json"

# Check if timestamp format exists
if [ -f "${TIMESTAMP}_dashboard.html" ]; then
  echo "✅ SUCCESS: Timestamp dashboard created: ${TIMESTAMP}_dashboard.html"
else
  echo "❌ FAILURE: Timestamp dashboard not created"
  exit 1
fi

# Check if multiple files were created
DASHBOARD_FILES=(*_dashboard.html)
DASHBOARD_COUNT="${#DASHBOARD_FILES[@]}"

if [ "$DASHBOARD_COUNT" -gt 1 ]; then
  echo "✅ SUCCESS: Multiple dashboard files created (${DASHBOARD_COUNT})"
  ls -la *_dashboard.html
else
  echo "⚠️ WARNING: Only one dashboard file created"
  ls -la *_dashboard.html
fi

echo "==============================================================="
echo "Test 3: Emergency Fallback Generation"
echo "==============================================================="
rm -f *_dashboard.html
mkdir -p empty_data
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

# Run the R script with minimal data to trigger emergency fallback
Rscript generate_dashboard.R \
  --outDir="." \
  --resultDir="./empty_data" \
  --outputFile="${PROJECT_ID}_dashboard.html" \
  --title="Emergency Dashboard" \
  --templateFile="nonexistent_template.html" \
  --qualityDataPath="empty_data/nonexistent.json"

# Check if emergency dashboard exists
if [ -f "${PROJECT_ID}_dashboard.html" ]; then
  echo "✅ SUCCESS: Emergency dashboard created: ${PROJECT_ID}_dashboard.html"
else
  echo "❌ FAILURE: Emergency dashboard not created"
  exit 1
fi

# Check if timestamp format was created
TIMESTAMP_FILES=(*_dashboard.html)
TIMESTAMP_COUNT="${#TIMESTAMP_FILES[@]}"

if [ "$TIMESTAMP_COUNT" -gt 1 ]; then
  echo "✅ SUCCESS: Multiple emergency dashboards created (${TIMESTAMP_COUNT})"
  ls -la *_dashboard.html
else
  echo "⚠️ WARNING: Only one emergency dashboard created"
  ls -la *_dashboard.html
fi

# Final results
echo "==============================================================="
echo "All dashboard generation tests completed"
echo "Test files are in: $TEST_DIR"
echo "==============================================================="