process GENERATE_DASHBOARD {
    tag "Generate dashboard"
    label 'process_medium'
    shell "/bin/bash"
    container 'foodnet.sif'

    input:
    path ir_outputs
    val resultDir
    val projID
    path dashboardTemplate
    val dashboardScript

    output:
    path "${projID}_dashboard.html", emit: dashboard
    path "${projID}_data_quality.log", optional: true, emit: quality_log
    path "${projID}_data_quality.json", optional: true, emit: quality_json
    publishDir "${params.outdir}/${params.projID}", mode: params.publish_dir_mode

    script:
    """
    # Create a data quality assessment log
    echo "=========================================" > ${projID}_data_quality.log
    echo "  FoodNet Trends Data Quality Assessment" >> ${projID}_data_quality.log
    echo "  Project ID: ${projID}" >> ${projID}_data_quality.log
    echo "  Generated: `date`" >> ${projID}_data_quality.log
    echo "=========================================" >> ${projID}_data_quality.log
    echo "" >> ${projID}_data_quality.log
    
    # Initialize JSON structure for data quality metadata
    ISO_DATE=`date -Iseconds`
    cat > ${projID}_data_quality.json << EOF
    {
      "projectId": "${projID}",
      "generationDate": "\$ISO_DATE",
      "dataQuality": {
        "usesPlaceholderData": false,
        "affectedPathogens": [],
        "warnings": [],
        "dataConsistency": {}
      }
    }
    EOF
    
    # Check for placeholder data warnings in any result files
    echo "Checking for data quality issues..." >> ${projID}_data_quality.log
    
    # Create placeholder and warnings flags
    PLACEHOLDER_DATA_FOUND=0
    PATHOGEN_WARNINGS=""
    AFFECTED_PATHOGENS=""
    DATA_QUALITY_WARNINGS=0
    
    # Look for placeholder warnings in summary files
    find . -name "*_data_summary.txt" > summary_files.txt
    if [ -s summary_files.txt ]; then
        while read -r SUMMARY_FILE; do
            if [ -f "\$SUMMARY_FILE" ]; then
                # Extract pathogen name
                BASE_NAME=`basename "\$SUMMARY_FILE"`
                PATHOGEN=`echo "\$BASE_NAME" | sed 's/_data_summary.txt//'`
                
                # Check for placeholders
                if grep -q "PLACEHOLDER DATA\\|placeholder data\\|SYNTHETIC DATA" "\$SUMMARY_FILE"; then
                    # Update pathogen list
                    if [ -z "\$AFFECTED_PATHOGENS" ]; then
                        AFFECTED_PATHOGENS="\\\"\$PATHOGEN\\\""
                    else
                        AFFECTED_PATHOGENS="\$AFFECTED_PATHOGENS, \\\"\$PATHOGEN\\\""
                    fi
                    
                    # Add warning
                    echo "WARNING: Placeholder data detected for pathogen: \$PATHOGEN" >> ${projID}_data_quality.log
                    PATHOGEN_WARNINGS="\$PATHOGEN_WARNINGS\n- \$PATHOGEN: Uses placeholder data (results NOT suitable for production use)"
                    DATA_QUALITY_WARNINGS=\$((DATA_QUALITY_WARNINGS + 1))
                    PLACEHOLDER_DATA_FOUND=1
                fi
            fi
        done < summary_files.txt
    fi
    
    # Update JSON if placeholder data found
    if [ \$PLACEHOLDER_DATA_FOUND -eq 1 ]; then
        # Update JSON to indicate placeholder data is being used
        sed -i 's/"usesPlaceholderData": false/"usesPlaceholderData": true/' ${projID}_data_quality.json
        # Add the affected pathogens to the JSON
        sed -i "s/\\\"affectedPathogens\\\": \\[\\]/\\\"affectedPathogens\\\": [\$AFFECTED_PATHOGENS]/" ${projID}_data_quality.json
        
        # Add warning to log
        echo "WARNING: Placeholder data detected" >> ${projID}_data_quality.log
        echo "*** RESULTS MAY NOT BE SUITABLE FOR PRODUCTION USE ***" >> ${projID}_data_quality.log
        echo "" >> ${projID}_data_quality.log
        
        # Create a warning banner for HTML
        WARNING_BANNER="<div style='background-color: #f8d7da; color: #721c24; padding: 15px; margin: 20px 0; border: 1px solid #f5c6cb; border-radius: 5px;'><h3 style='margin-top:0'>⚠️ WARNING: Placeholder Data Detected</h3><p>This analysis contains placeholder data which may not represent real-world conditions.</p><p><strong>Results are NOT suitable for production or public health decision-making!</strong></p></div>"
        
        # Insert warning into template
        if [ -f "${dashboardTemplate}" ]; then
            cp "${dashboardTemplate}" modified_template.html
            sed -i "s|<body>|<body>\\n\$WARNING_BANNER|" modified_template.html
            TEMPLATE_TO_USE="modified_template.html"
        else
            # Fallback if template missing
            cat > fallback_template.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Fallback Template</title></head>
<body>
\$WARNING_BANNER
<h1>FoodNet Trends Analysis</h1>
<p>This is a fallback template due to missing template file.</p>
EOF
            TEMPLATE_TO_USE="fallback_template.html"
        fi
    else
        # No placeholder warnings
        echo "No placeholder data warnings detected." >> ${projID}_data_quality.log
        TEMPLATE_TO_USE="${dashboardTemplate}"
    fi
    
    # Count IR files for data consistency
    echo "" >> ${projID}_data_quality.log
    echo "Checking data consistency..." >> ${projID}_data_quality.log
    IR_FILE_COUNT=`ls -1 *_IRCatch.csv 2>/dev/null | wc -l`
    echo "Found \$IR_FILE_COUNT incidence rate files." >> ${projID}_data_quality.log
    
    # Add to JSON
    sed -i "s/\"dataConsistency\": {}/\"dataConsistency\": {\\n      \\\"incidenceRateFiles\\\": \$IR_FILE_COUNT,\\n      \\\"dataQualityWarnings\\\": \$DATA_QUALITY_WARNINGS\\n    }/" ${projID}_data_quality.json
    
    # Generate the dashboard with improved error handling
    echo "Generating dashboard with template: \$TEMPLATE_TO_USE" >> ${projID}_data_quality.log
    
    # Verify dashboard script exists and run it
    if [ ! -f "${dashboardScript}" ]; then
        echo "ERROR: Dashboard script not found: ${dashboardScript}" >> ${projID}_data_quality.log
        # Create fallback dashboard
        cat > ${projID}_dashboard.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Error</title></head>
<body>
<h1>FoodNet Trends Analysis Dashboard</h1>
<h2>Error: Script Not Found</h2>
<p>The dashboard generation script was not found. This is a simplified fallback dashboard.</p>
<p>Analysis completed at: `date`</p>
<p>Project ID: ${projID}</p>
</body>
</html>
EOF
    else
        # Run script with error handling
        set +e
        Rscript \\
          "${dashboardScript}" \\
          --outDir="${resultDir}" \\
          --resultDir="${resultDir}" \\
          --outputFile="${projID}_dashboard.html" \\
          --title="FoodNet Trends Analysis: ${projID}" \\
          --templateFile="\$TEMPLATE_TO_USE" \\
          --qualityDataPath="${projID}_data_quality.json"
          
        SCRIPT_EXIT_CODE=\$?
        set -e
        
        # Handle script errors
        if [ \$SCRIPT_EXIT_CODE -ne 0 ]; then
            echo "ERROR: Dashboard generation script failed with exit code \$SCRIPT_EXIT_CODE" >> ${projID}_data_quality.log
            # Create fallback dashboard
            cat > ${projID}_dashboard.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Error</title></head>
<body>
<h1>FoodNet Trends Analysis Dashboard</h1>
<h2>Error: Script Failed</h2>
<p>The dashboard generation script failed with exit code \$SCRIPT_EXIT_CODE.</p>
<p>Analysis completed at: `date`</p>
<p>Project ID: ${projID}</p>
</body>
</html>
EOF
        fi
    fi
    
    # Verify dashboard was created
    if [ ! -f "${projID}_dashboard.html" ]; then
        echo "ERROR: Dashboard file was not created" >> ${projID}_data_quality.log
        # Create a simple fallback dashboard
        cat > ${projID}_dashboard.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Error</title></head>
<body>
<h1>FoodNet Trends Analysis Dashboard</h1>
<h2>Error: Missing Output</h2>
<p>The dashboard output file was not created properly.</p>
<p>Analysis completed at: `date`</p>
<p>Project ID: ${projID}</p>
</body>
</html>
EOF
    fi
    
    # Final note
    echo "Dashboard generation completed at `date`" >> ${projID}_data_quality.log
    """
}