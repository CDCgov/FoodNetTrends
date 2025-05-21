process GENERATE_DASHBOARD {
    // Add error retry strategy to handle transient failures
    errorStrategy 'retry'
    maxRetries 2
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
    publishDir "${params.outdir}/${projID}", mode: params.publish_dir_mode

    script:
    """
    # Create a data quality assessment log
    echo "=========================================" > ${projID}_data_quality.log
    echo "  FoodNet Trends Data Quality Assessment" >> ${projID}_data_quality.log
    echo "  Project ID: ${projID}" >> ${projID}_data_quality.log
    # Get date safely with shell function
    DATE_NOW=\$(date)
    echo "  Generated: \$DATE_NOW" >> ${projID}_data_quality.log
    echo "=========================================" >> ${projID}_data_quality.log
    echo "" >> ${projID}_data_quality.log
    
    # Initialize JSON structure for data quality metadata
    # Get ISO date safely with direct assignment
    ISO_DATE=\$(date -Iseconds)
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
                BASE_NAME=\$(basename "\$SUMMARY_FILE")
                PATHOGEN=\$(echo "\$BASE_NAME" | sed 's/_data_summary.txt//')
                
                # Check for placeholders
                if grep -q "PLACEHOLDER DATA\\|placeholder data\\|SYNTHETIC DATA" "\$SUMMARY_FILE"; then
                    # Update pathogen list
                    if [ -z "\$AFFECTED_PATHOGENS" ]; then
                        AFFECTED_PATHOGENS="\\\"\$PATHOGEN\\\""
                    else
                        AFFECTED_PATHOGENS="\$AFFECTED_PATHOGENS, \\\"\$PATHOGEN\\\""
                    fi
                    
                    # Add warning - safely handle variable
                    if [ -n "\$PATHOGEN" ]; then
                        echo "WARNING: Placeholder data detected for pathogen: \$PATHOGEN" >> ${projID}_data_quality.log
                        PATHOGEN_WARNINGS="\$PATHOGEN_WARNINGS\n- \$PATHOGEN: Uses placeholder data (results NOT suitable for production use)"
                    else
                        echo "WARNING: Placeholder data detected in unidentified file" >> ${projID}_data_quality.log
                        PATHOGEN_WARNINGS="\$PATHOGEN_WARNINGS\n- UNKNOWN: Uses placeholder data (results NOT suitable for production use)"
                    fi
                    DATA_QUALITY_WARNINGS=\$((DATA_QUALITY_WARNINGS + 1))
                    PLACEHOLDER_DATA_FOUND=1
                fi
            fi
        done < summary_files.txt
    fi
    
    # First set default template
    TEMPLATE_TO_USE="${dashboardTemplate}"
    
    # Handle placeholder data - complete rewrite to avoid nested if statements
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
        
        # Template handling - avoid nesting if statements
        if [ -f "${dashboardTemplate}" ]; then
            cp "${dashboardTemplate}" modified_template.html
            sed -i "s|<body>|<body>\\n\$WARNING_BANNER|" modified_template.html
            TEMPLATE_TO_USE="modified_template.html"
        else
            # Create a simple fallback template
            echo "<!DOCTYPE html>" > fallback_template.html
            echo "<html>" >> fallback_template.html
            echo "<head><title>FoodNet Trends Dashboard - Fallback Template</title></head>" >> fallback_template.html
            echo "<body>" >> fallback_template.html
            echo "\$WARNING_BANNER" >> fallback_template.html
            echo "<h1>FoodNet Trends Analysis</h1>" >> fallback_template.html
            echo "<p>This is a fallback template due to missing template file.</p>" >> fallback_template.html
            echo "</body>" >> fallback_template.html
            echo "</html>" >> fallback_template.html
            TEMPLATE_TO_USE="fallback_template.html"
        fi
    else
        # No placeholder warnings
        echo "No placeholder data warnings detected." >> ${projID}_data_quality.log
    fi
    
    # Count IR files for data consistency
    echo "" >> ${projID}_data_quality.log
    echo "Checking data consistency..." >> ${projID}_data_quality.log
    IR_FILE_COUNT=\$(ls -1 *_IRCatch.csv 2>/dev/null | wc -l)
    echo "Found \$IR_FILE_COUNT incidence rate files." >> ${projID}_data_quality.log
    
    # Add to JSON
    sed -i "s/\"dataConsistency\": {}/\"dataConsistency\": {\\n      \\\"incidenceRateFiles\\\": \$IR_FILE_COUNT,\\n      \\\"dataQualityWarnings\\\": \$DATA_QUALITY_WARNINGS\\n    }/" ${projID}_data_quality.json
    
    # Create temporary debug info for troubleshooting
    echo "Creating debug information for dashboard generation" >> ${projID}_data_quality.log
    mkdir -p debuginfo
    ls -la > debuginfo/directory_listing.txt
    ls -la ${resultDir} > debuginfo/resultdir_listing.txt 2>/dev/null || echo "Could not list result directory" > debuginfo/resultdir_listing.txt
    
    # Collect info about input and result files
    find . -name "*_IRCatch.csv" -o -name "*_summary.txt" -o -name "*EstIRRCatch*.csv" > debuginfo/result_files.txt
    RESULT_COUNT=\$(wc -l < debuginfo/result_files.txt)
    echo "Found \$RESULT_COUNT result files" >> ${projID}_data_quality.log
    
    # Check for common result files that should exist without using shell variables
    # Use a simpler approach to check each pathogen individually
    ls CAMPYLOBACTER_*_IRCatch.csv 1>/dev/null 2>&1 && echo "Found IR files for CAMPYLOBACTER" >> debuginfo/file_checks.txt || echo "Missing IR files for CAMPYLOBACTER" >> debuginfo/file_checks.txt
    ls CYCLOSPORA_*_IRCatch.csv 1>/dev/null 2>&1 && echo "Found IR files for CYCLOSPORA" >> debuginfo/file_checks.txt || echo "Missing IR files for CYCLOSPORA" >> debuginfo/file_checks.txt
    ls SALMONELLA_*_IRCatch.csv 1>/dev/null 2>&1 && echo "Found IR files for SALMONELLA" >> debuginfo/file_checks.txt || echo "Missing IR files for SALMONELLA" >> debuginfo/file_checks.txt
    ls SHIGELLA_*_IRCatch.csv 1>/dev/null 2>&1 && echo "Found IR files for SHIGELLA" >> debuginfo/file_checks.txt || echo "Missing IR files for SHIGELLA" >> debuginfo/file_checks.txt
    
    # Generate the dashboard with improved error handling
    echo "Generating dashboard with template: \$TEMPLATE_TO_USE" >> ${projID}_data_quality.log
    
    # Verify dashboard script exists and run it
    if [ ! -f "${dashboardScript}" ]; then
        echo "ERROR: Dashboard script not found: ${dashboardScript}" >> ${projID}_data_quality.log
        # Create fallback dashboard with proper closing tags
        cat > ${projID}_dashboard.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Error</title></head>
<body>
<h1>FoodNet Trends Analysis Dashboard</h1>
<h2>Error: Script Not Found</h2>
<p>The dashboard generation script was not found. This is a simplified fallback dashboard.</p>
<p>Analysis completed at: \$(date)</p>
<p>Project ID: ${projID}</p>
</body>
</html>
EOF
    else
        # Run script with debugging environment and error handling
        echo "Preparing dashboard data..." >> ${projID}_data_quality.log
        
        # Set environment variables to help R script with debugging
        export R_DEBUG_LEVEL=1
        export R_LIBS_USER="${params.outdir}/Rlibs"
        export TRENDY_DEBUG=1
        
        # Create directory for R libs if it doesn't exist (helps with permissions)
        mkdir -p "${params.outdir}/Rlibs" 2>/dev/null
        
        # Log which result files we have before running R
        echo "Incidence rate files available:" >> ${projID}_data_quality.log
        find . -name "*_IRCatch.csv" -ls >> ${projID}_data_quality.log 2>/dev/null || echo "  None found" >> ${projID}_data_quality.log
        
        # Run script with error handling and detailed logging
        set +e
        Rscript \\
          "${dashboardScript}" \\
          --outDir="${resultDir}" \\
          --resultDir="${resultDir}" \\
          --outputFile="${projID}_dashboard.html" \\
          --title="FoodNet Trends Analysis: ${projID}" \\
          --templateFile="\$TEMPLATE_TO_USE" \\
          --qualityDataPath="${projID}_data_quality.json" 2>&1 | tee -a dashboard_generation.log
          
        SCRIPT_EXIT_CODE=\$?
        set -e
        
        # Handle script errors
        if [ \$SCRIPT_EXIT_CODE -ne 0 ]; then
            echo "ERROR: Dashboard generation script failed with exit code \$SCRIPT_EXIT_CODE" >> ${projID}_data_quality.log
            # Save the log file for debugging
            cp dashboard_generation.log debuginfo/
            
            # Create a more informative fallback dashboard - using simpler commands
            cat > ${projID}_dashboard.html << EOF
<!DOCTYPE html>
<html>
<head>
  <title>FoodNet Trends Dashboard - Error</title>
  <style>
    body { font-family: Arial, sans-serif; line-height: 1.6; padding: 20px; max-width: 1000px; margin: 0 auto; }
    .error-banner { background-color: #f8d7da; color: #721c24; padding: 15px; border-radius: 5px; }
    .debug-info { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-top: 20px; }
    h1 { color: #0066cc; }
    pre { background-color: #f1f1f1; padding: 10px; overflow-x: auto; }
  </style>
</head>
<body>
  <h1>FoodNet Trends Analysis Dashboard</h1>
  <div class="error-banner">
    <h2>Error: Dashboard Generation Failed</h2>
    <p>The dashboard generation script failed with exit code \$SCRIPT_EXIT_CODE.</p>
  </div>
  
  <div class="debug-info">
    <h3>Debug Information</h3>
    <p><strong>Project ID:</strong> ${projID}</p>
    
    <h4>Possible solutions:</h4>
    <ul>
      <li>Check that all required result files exist</li>
      <li>Verify that the census files are properly formatted</li>
      <li>Run the workflow with proper input files instead of placeholders</li>
    </ul>
  </div>
</body>
</html>
EOF

            # Add runtime information separately - properly escaped for Nextflow
            CURRENT_DATE=\$(date)
            echo "<script>document.getElementsByClassName('debug-info')[0].insertAdjacentHTML('afterbegin', '<p><strong>Analysis time:</strong> \$CURRENT_DATE</p>');</script>" >> ${projID}_dashboard.html
            
            # Add error log information - properly escaped for Nextflow
            ERROR_LOG=\$(tail -n 10 dashboard_generation.log 2>/dev/null || echo "No log file available")
            echo "<script>document.getElementsByClassName('debug-info')[0].insertAdjacentHTML('beforeend', '<h4>Script Errors:</h4><pre>\$ERROR_LOG</pre>');</script>" >> ${projID}_dashboard.html
        fi
    fi
    
    # Verify dashboard was created and collect debug info
    if [ ! -f "${projID}_dashboard.html" ]; then
        echo "ERROR: Dashboard file was not created" >> ${projID}_data_quality.log
        
        # Collect diagnostic information
        echo "Collecting diagnostic information for troubleshooting" >> ${projID}_data_quality.log
        R --version > debuginfo/r_version.txt 2>&1 || echo "R not available" > debuginfo/r_version.txt
        env > debuginfo/environment.txt
        df -h > debuginfo/disk_space.txt
        
        # Save R package info if possible
        Rscript -e "installed.packages()[,c('Package', 'Version')]" > debuginfo/r_packages.txt 2>/dev/null || echo "Cannot list R packages" > debuginfo/r_packages.txt
        
        # Create a simple fallback dashboard with diagnostic info
        cat > ${projID}_dashboard.html << EOF
<!DOCTYPE html>
<html>
<head>
  <title>FoodNet Trends Dashboard - Error</title>
  <style>
    body { font-family: Arial, sans-serif; line-height: 1.6; padding: 20px; max-width: 1000px; margin: 0 auto; }
    .error-banner { background-color: #f8d7da; color: #721c24; padding: 15px; border-radius: 5px; }
    .debug-info { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-top: 20px; }
    h1 { color: #0066cc; }
  </style>
</head>
<body>
  <h1>FoodNet Trends Analysis Dashboard</h1>
  <div class="error-banner">
    <h2>Error: Missing Dashboard Output</h2>
    <p>The dashboard HTML file was not created properly.</p>
  </div>
  
  <div class="debug-info">
    <h3>Analysis Information</h3>
    <p><strong>Analysis completed at:</strong> \$(date)</p>
    <p><strong>Project ID:</strong> ${projID}</p>
    
    <h4>Diagnostic Information</h4>
    <p>Diagnostic information has been saved to the 'debuginfo' directory.</p>
    
    <h4>Possible solutions:</h4>
    <ul>
      <li>Check that all required input files are available</li>
      <li>Verify that the analysis produced valid result files</li>
      <li>Run the workflow with proper input files instead of placeholders if applicable</li>
      <li>Check the disk space and permissions in the output directory</li>
    </ul>
  </div>
</body>
</html>
EOF
    else
        # Dashboard was created, copy debug info if available
        if [ -d "debuginfo" ]; then
            mkdir -p "${params.outdir}/${projID}/debuginfo"
            cp -r debuginfo/* "${params.outdir}/${projID}/debuginfo/" 2>/dev/null
        fi
    fi
    
    # Final note
    echo "Dashboard generation completed at \$(date)" >> ${projID}_data_quality.log
    """
}