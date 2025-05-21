#!/bin/bash
# Test script for census file path metadata storage and retrieval

# Set up variables
MMWR_FILE="/path/to/mmwr/data.sas7bdat"  # Replace with actual path
CENSUS_B="/path/to/census/bacterial.csv"  # Replace with actual path
CENSUS_P="/path/to/census/parasitic.csv"  # Replace with actual path
PROJ_ID="census_metadata_test_$(date +%Y%m%d_%H%M%S)"
OUTPUT_DIR="results"

echo "====================================================="
echo "Testing census file path metadata storage and retrieval"
echo "====================================================="
echo "MMWR File: $MMWR_FILE"
echo "Census Bacterial: $CENSUS_B"
echo "Census Parasitic: $CENSUS_P"
echo "Project ID: $PROJ_ID"
echo "Output Directory: $OUTPUT_DIR"
echo "====================================================="

# First run preprocessing to generate metadata with census file paths
echo "Step 1: Running preprocessing to generate metadata..."
nextflow run workflows/preprocess.nf \
  --mmwrFile "$MMWR_FILE" \
  --censusFileB "$CENSUS_B" \
  --censusFileP "$CENSUS_P" \
  --projID "$PROJ_ID" \
  --outdir "$OUTPUT_DIR" \
  --generate_metadata true

# Check if preprocessing succeeded
if [ $? -ne 0 ]; then
  echo "Preprocessing failed! Check logs for errors."
  exit 1
fi

# Check if metadata file was created
METADATA_FILE="$OUTPUT_DIR/preprocessed/${PROJ_ID}_metadata.json"
if [ ! -f "$METADATA_FILE" ]; then
  echo "Metadata file not created: $METADATA_FILE"
  exit 1
fi

# Check if metadata contains census file paths
echo "Checking metadata file for census file paths..."
grep -q "census_file_bacterial" "$METADATA_FILE"
if [ $? -ne 0 ]; then
  echo "Metadata does not contain bacterial census file path!"
  exit 1
fi

grep -q "census_file_parasitic" "$METADATA_FILE"
if [ $? -ne 0 ]; then
  echo "Metadata does not contain parasitic census file path!"
  exit 1
fi

echo "Metadata file contains both census file paths. Good!"

# Now run the main workflow without specifying census files
echo "Step 2: Running main workflow without specifying census files..."
echo "This should automatically get census files from metadata"

nextflow run main.nf \
  --mmwrFile "$OUTPUT_DIR/preprocessed/${PROJ_ID}.csv" \
  --preprocessed true \
  --projID "$PROJ_ID" \
  --outdir "$OUTPUT_DIR" \
  --pathogen "CAMPYLOBACTER,SALMONELLA" \
  --skip_dashboard true

# Check if main workflow succeeded
if [ $? -ne 0 ]; then
  echo "Main workflow failed! Check logs for errors."
  exit 1
fi

echo "====================================================="
echo "Test completed successfully!"
echo "Metadata file: $METADATA_FILE"
echo "====================================================="

# Optionally display the metadata file
echo "Metadata file contents:"
cat "$METADATA_FILE"