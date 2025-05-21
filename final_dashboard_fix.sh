#!/bin/bash
# =========================================================================
# Final Dashboard Fix - Standalone Dashboard Creator
# =========================================================================
#
# Purpose:
#   Last resort utility to create a dashboard from existing results
#   when the pipeline's dashboard generation fails.
#
# Usage:
#   ./final_dashboard_fix.sh [project_id]
#
# =========================================================================

# Strict error handling
set -euo pipefail

# Get project ID from argument or use timestamp
PROJECT_ID=${1:-$(date +%Y%m%d_%H%M%S)}
OUTPUT_DIR="output"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

echo "========================================="
echo "  Final Dashboard Fix Utility"
echo "========================================="
echo "Project ID: $PROJECT_ID"
echo "Output Dir: $OUTPUT_DIR"
echo "Generated:  $(date)"
echo "========================================="
echo

# Verify output directory exists
if [ ! -d "$OUTPUT_DIR/$PROJECT_ID" ]; then
  echo "ERROR: Output directory does not exist: $OUTPUT_DIR/$PROJECT_ID"
  echo "Creating directory structure..."
  mkdir -p "$OUTPUT_DIR/$PROJECT_ID"
fi

# Check for result files
echo "Checking for result files..."
RESULT_FILES=$(find "$OUTPUT_DIR/$PROJECT_ID" -name "*_IRCatch.csv" 2>/dev/null | wc -l)
echo "Found $RESULT_FILES IR catch files"

if [ "$RESULT_FILES" -eq 0 ]; then
  echo "WARNING: No result files found. Dashboard will be minimal."
fi

# Ensure dashboard directory exists
mkdir -p "$OUTPUT_DIR/$PROJECT_ID"

# Create dashboard file paths
DASHBOARD_PATH="$OUTPUT_DIR/$PROJECT_ID/${PROJECT_ID}_dashboard.html"
DASHBOARD_ALT="$OUTPUT_DIR/$PROJECT_ID/${TIMESTAMP}_dashboard.html"

# Generate comprehensive dashboard
echo "Creating comprehensive dashboard at: $DASHBOARD_PATH"

# List pathogens from directory names or files
PATHOGENS=$(find "$OUTPUT_DIR/$PROJECT_ID" -name "*_IRCatch.csv" 2>/dev/null | sed -E 's/.*\/([A-Z]+)_.*/\1/g' | sort | uniq | tr '\n' ',' | sed 's/,$//')
if [ -z "$PATHOGENS" ]; then
  PATHOGENS="CAMPYLOBACTER,CYCLOSPORA,SALMONELLA,SHIGELLA,STEC,VIBRIO,YERSINIA"
fi

# Create HTML dashboard with detailed content
cat > "$DASHBOARD_PATH" << EOF
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>FoodNet Trends Dashboard - $PROJECT_ID</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      line-height: 1.6;
      color: #333;
      max-width: 1200px;
      margin: 0 auto;
      padding: 20px;
    }
    .header {
      background-color: #0054ad;
      color: white;
      padding: 20px;
      border-radius: 5px;
      margin-bottom: 20px;
    }
    .warning {
      background-color: #fff3cd;
      border: 1px solid #ffeeba;
      color: #856404;
      padding: 15px;
      border-radius: 5px;
      margin: 20px 0;
    }
    .info {
      background-color: #e2f0fb;
      border: 1px solid #bee5eb;
      color: #0c5460;
      padding: 15px;
      border-radius: 5px;
      margin: 20px 0;
    }
    .section {
      margin-bottom: 30px;
      padding: 20px;
      background-color: #f9f9f9;
      border-radius: 5px;
      box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    table {
      width: 100%;
      border-collapse: collapse;
      margin-bottom: 20px;
    }
    th, td {
      padding: 12px;
      border: 1px solid #ddd;
      text-align: left;
    }
    th {
      background-color: #f2f2f2;
    }
    tr:nth-child(even) {
      background-color: #f9f9f9;
    }
  </style>
</head>
<body>
  <div class="header">
    <h1>FoodNet Trends Analysis Dashboard</h1>
    <p>Project ID: $PROJECT_ID</p>
    <p>Generated on: $(date)</p>
  </div>

  <div class="section">
    <h2>Analysis Overview</h2>
    <p>This dashboard provides a summary of the FoodNet Trends analysis results.</p>
    <p>Note: This is a standalone dashboard created by the final_dashboard_fix.sh utility.</p>
    
    <div class="info">
      <h3>Analysis Information</h3>
      <table>
        <tr>
          <td><strong>Project ID:</strong></td>
          <td>$PROJECT_ID</td>
        </tr>
        <tr>
          <td><strong>Result Files Found:</strong></td>
          <td>$RESULT_FILES files</td>
        </tr>
        <tr>
          <td><strong>Pathogens Analyzed:</strong></td>
          <td>$PATHOGENS</td>
        </tr>
        <tr>
          <td><strong>Generation Time:</strong></td>
          <td>$(date)</td>
        </tr>
      </table>
    </div>
  </div>
  
  <div class="section">
    <h2>Results</h2>
    <p>FoodNet Trends analysis was performed on the specified pathogens.</p>
    <p>Result files are available in the output directory: <code>$OUTPUT_DIR/$PROJECT_ID</code></p>
    
    <h3>Pathogens Included in Analysis</h3>
    <ul>
EOF

# Add pathogen list from found files
for PATHOGEN in $(echo $PATHOGENS | tr ',' '\n'); do
  echo "      <li>$PATHOGEN</li>" >> "$DASHBOARD_PATH"
done

# Complete the HTML
cat >> "$DASHBOARD_PATH" << EOF
    </ul>
  </div>
  
  <div class="section">
    <h2>Data Files</h2>
    <p>The following result files were found:</p>
    <ul>
EOF

# Find and list result files
find "$OUTPUT_DIR/$PROJECT_ID" -name "*_IRCatch.csv" 2>/dev/null | while read -r file; do
  basename=$(basename "$file")
  echo "      <li>$basename</li>" >> "$DASHBOARD_PATH"
done

# Add placeholder if no files found
if [ "$RESULT_FILES" -eq 0 ]; then
  echo "      <li>No result files found</li>" >> "$DASHBOARD_PATH"
fi

# Complete the HTML
cat >> "$DASHBOARD_PATH" << EOF
    </ul>
  </div>
  
  <div class="section">
    <h2>Next Steps</h2>
    <p>To analyze these results further:</p>
    <ol>
      <li>Review the incidence rate files in the output directory</li>
      <li>Check the summary files for details about each pathogen</li>
      <li>For detailed visualizations, import the CSV files into your preferred data analysis tool</li>
    </ol>
  </div>
  
  <footer style="margin-top: 50px; border-top: 1px solid #ddd; padding-top: 20px; color: #777;">
    <p>FoodNet Trends Analysis | Generated on: $(date)</p>
    <p><small>Created with final_dashboard_fix.sh</small></p>
  </footer>
</body>
</html>
EOF

# Create a copy with timestamp-based name for compatibility
cp "$DASHBOARD_PATH" "$DASHBOARD_ALT"

# Verify files were created
if [ -f "$DASHBOARD_PATH" ] && [ -f "$DASHBOARD_ALT" ]; then
  echo "Success! Dashboard files created:"
  echo "- $DASHBOARD_PATH"
  echo "- $DASHBOARD_ALT"
  echo
  echo "File sizes:"
  ls -lh "$DASHBOARD_PATH" "$DASHBOARD_ALT"
else
  echo "Error: Failed to create dashboard files."
  exit 1
fi

echo
echo "Dashboard generation completed successfully."
echo "If you had dashboard issues in the pipeline, consider running:"
echo "  ./run_workflow_hpc.sh with increased memory settings"