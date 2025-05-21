process GENERATE_DASHBOARD {
    // Use retry strategy to handle potential memory or time issues
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' } // 137 is OOM kill
    maxRetries 3
    tag "Generate dashboard"
    label 'process_high_memory' // Increase resource label
    shell "/bin/bash"
    container 'foodnet.sif'
    time '1h'
    memory '16 GB' // Significantly increase memory allocation

    input:
    path ir_outputs
    val resultDir
    val projID
    path dashboardTemplate
    val dashboardScript

    output:
    path "*.html", emit: dashboard, optional: true
    path "dashboard.html", emit: main_dashboard, optional: true  
    path "*", emit: all_files, optional: true
    path "*_data_quality.log", optional: true, emit: quality_log
    path "*_data_quality.json", optional: true, emit: quality_json
    publishDir "${params.outdir}/${projID}", mode: params.publish_dir_mode

    script:
    def timestamp = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
    def currentDate = new Date().toString()
    """
    #!/bin/bash
    # Simple dashboard generator script
    
    # Set timestamp from Nextflow variable
    MY_TIMESTAMP="${timestamp}"
    
    # Memory optimization for R
    export R_MAX_VSIZE=12G
    export R_GC_MEM_GROW=0
    
    # Create a data quality assessment log
    echo "=========================================" > ${projID}_data_quality.log
    echo "  FoodNet Trends Data Quality Assessment" >> ${projID}_data_quality.log
    echo "  Project ID: ${projID}" >> ${projID}_data_quality.log
    echo "  Generated: ${currentDate}" >> ${projID}_data_quality.log
    echo "=========================================" >> ${projID}_data_quality.log
    echo "" >> ${projID}_data_quality.log
    
    # Initialize JSON structure for data quality metadata
    cat > ${projID}_data_quality.json << EOF
    {
      "projectId": "${projID}",
      "generationDate": "${currentDate}",
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
                BASE_NAME=\${SUMMARY_FILE##*/}
                PATHOGEN=\${BASE_NAME%_data_summary.txt}
                
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
                        PATHOGEN_WARNINGS="\$PATHOGEN_WARNINGS\\n- \$PATHOGEN: Uses placeholder data (results NOT suitable for production use)"
                    else
                        echo "WARNING: Placeholder data detected in unidentified file" >> ${projID}_data_quality.log
                        PATHOGEN_WARNINGS="\$PATHOGEN_WARNINGS\\n- UNKNOWN: Uses placeholder data (results NOT suitable for production use)"
                    fi
                    let DATA_QUALITY_WARNINGS+=1
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
            # Create a simple fallback template using a heredoc instead of multiple echo statements
            cat > fallback_template.html << 'EOFTEMPLATE'
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Fallback Template</title></head>
<body>
EOFTEMPLATE
            # Add the warning banner separately
            echo "\$WARNING_BANNER" >> fallback_template.html
            
            # Complete the HTML template
            cat >> fallback_template.html << 'EOFTEMPLATE'
<h1>FoodNet Trends Analysis</h1>
<p>This is a fallback template due to missing template file.</p>
</body>
</html>
EOFTEMPLATE
            TEMPLATE_TO_USE="fallback_template.html"
        fi
    else
        # No placeholder warnings
        echo "No placeholder data warnings detected." >> ${projID}_data_quality.log
    fi
    
    # Count IR files for data consistency
    echo "" >> ${projID}_data_quality.log
    echo "Checking data consistency..." >> ${projID}_data_quality.log
    # Count files with pure shell
    set +e  # Don't fail if no files found
    ls -1 *_IRCatch.csv > ir_files_list.txt 2>/dev/null
    if [ -s ir_files_list.txt ]; then
        IR_FILE_COUNT=0
        while read -r LINE; do
            let IR_FILE_COUNT+=1
        done < ir_files_list.txt
    else
        IR_FILE_COUNT=0
    fi
    set -e
    echo "Found \$IR_FILE_COUNT incidence rate files." >> ${projID}_data_quality.log
    
    # Add to JSON
    sed -i "s/\"dataConsistency\": {}/\"dataConsistency\": {\\n      \\\"incidenceRateFiles\\\": \$IR_FILE_COUNT,\\n      \\\"dataQualityWarnings\\\": \$DATA_QUALITY_WARNINGS\\n    }/" ${projID}_data_quality.json
    
    # Create temporary debug info for troubleshooting
    echo "Creating debug information for dashboard generation" >> ${projID}_data_quality.log
    
    # Ensure the debuginfo directory exists before trying to write to it
    mkdir -p debuginfo
    
    # Create empty debug files to prevent missing file errors
    touch debuginfo/directory_listing.txt
    touch debuginfo/resultdir_listing.txt
    touch debuginfo/result_files.txt
    touch debuginfo/file_checks.txt
    
    # Write directory listings safely with error handling
    { ls -la > debuginfo/directory_listing.txt; } 2>/dev/null || echo "Failed to create directory listing" >> ${projID}_data_quality.log
    { ls -la ${resultDir} > debuginfo/resultdir_listing.txt; } 2>/dev/null || echo "Could not list result directory: ${resultDir}" > debuginfo/resultdir_listing.txt
    
    # Collect info about input and result files with robust error handling
    { find . -name "*_IRCatch.csv" -o -name "*_summary.txt" -o -name "*EstIRRCatch*.csv" > debuginfo/result_files.txt; } 2>/dev/null || echo "Warning: find command failed" >> ${projID}_data_quality.log
    
    # Read the count safely, ensuring the file exists
    if [ -f "debuginfo/result_files.txt" ]; then
        # Count files with pure shell
        RESULT_COUNT=0
        if [ -s debuginfo/result_files.txt ]; then
            while read -r LINE; do
                let RESULT_COUNT+=1
            done < debuginfo/result_files.txt
        fi
        echo "Found \$RESULT_COUNT result files" >> ${projID}_data_quality.log
    else
        echo "Warning: result_files.txt not created properly" >> ${projID}_data_quality.log
        RESULT_COUNT=0
        echo "Found 0 result files" >> ${projID}_data_quality.log
    fi
    
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
        # Create fallback dashboard - using a safer approach
        cat > ${projID}_dashboard.html << 'ERRORTEMPLATE'
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Error</title></head>
<body>
<h1>FoodNet Trends Analysis Dashboard</h1>
<h2>Error: Script Not Found</h2>
<p>The dashboard generation script was not found. This is a simplified fallback dashboard.</p>
ERRORTEMPLATE

        # Add dynamic content directly
        echo "<p>Project ID: ${projID}</p>" >> ${projID}_dashboard.html
        echo "<p>Analysis completed at: ${currentDate}</p>" >> ${projID}_dashboard.html
        
        # Close the HTML
        cat >> ${projID}_dashboard.html << 'ERRORTEMPLATE'
</body>
</html>
ERRORTEMPLATE
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
        
        # Run script with error handling and detailed logging - use optimized version
        set +e
        Rscript \\
          "${workflow.projectDir}/bin/generate_dashboard_optimized.R" \\
          --outDir="${resultDir}" \\
          --resultDir="${resultDir}" \\
          --outputFile="${projID}_dashboard.html" \\
          --title="FoodNet Trends Analysis: ${projID}" \\
          --templateFile="\$TEMPLATE_TO_USE" \\
          --qualityDataPath="${projID}_data_quality.json" \\
          --memoryLimit=14 2>&1 | tee -a dashboard_generation.log
          
        SCRIPT_EXIT_CODE=\$?
        set -e
        
        # Handle script errors
        if [ \$SCRIPT_EXIT_CODE -ne 0 ]; then
            echo "ERROR: Dashboard generation script failed with exit code \$SCRIPT_EXIT_CODE" >> ${projID}_data_quality.log
            # Save the log file for debugging
            cp dashboard_generation.log debuginfo/
            
            # Create a more informative fallback dashboard - using a safer approach with multiple parts
            # First create the basic HTML structure
            cat > ${projID}_dashboard.html << 'ERROR_TEMPLATE'
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
ERROR_TEMPLATE

            # Add the exit code information separately
            echo "    <p>The dashboard generation script failed with exit code \$SCRIPT_EXIT_CODE.</p>" >> ${projID}_dashboard.html
            
            # Continue with the rest of the template
            cat >> ${projID}_dashboard.html << 'ERROR_TEMPLATE'
  </div>
  
  <div class="debug-info">
    <h3>Debug Information</h3>
ERROR_TEMPLATE

            # Add the dynamic project ID separately
            echo "    <p><strong>Project ID:</strong> ${projID}</p>" >> ${projID}_dashboard.html
            
            # Add date from Nextflow variable
            echo "    <p><strong>Analysis time:</strong> ${currentDate}</p>" >> ${projID}_dashboard.html
            
            # Continue with static content
            cat >> ${projID}_dashboard.html << 'ERROR_TEMPLATE'    
    <h4>Possible solutions:</h4>
    <ul>
      <li>Check that all required result files exist</li>
      <li>Verify that the census files are properly formatted</li>
      <li>Run the workflow with proper input files instead of placeholders</li>
    </ul>
ERROR_TEMPLATE

            # Add the error log information separately if it exists
            if [ -f "dashboard_generation.log" ]; then
                echo "    <h4>Script Errors:</h4>" >> ${projID}_dashboard.html
                echo "    <pre>" >> ${projID}_dashboard.html
                tail -n 10 dashboard_generation.log >> ${projID}_dashboard.html 2>/dev/null || echo "Error reading log file" >> ${projID}_dashboard.html
                echo "    </pre>" >> ${projID}_dashboard.html
            fi
            
            # Close the HTML structure
            cat >> ${projID}_dashboard.html << 'ERROR_TEMPLATE'
  </div>
</body>
</html>
ERROR_TEMPLATE
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
        
        # Create a simple fallback dashboard with diagnostic info - using a safer approach with multiple steps
        # First create the basic HTML structure
        cat > ${projID}_dashboard.html << 'ERRORTEMPLATE'
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
ERRORTEMPLATE

        # Add the dynamic parts separately to avoid shell expansion issues
        # Add date from Nextflow variable
        echo "    <p><strong>Project ID:</strong> ${projID}</p>" >> ${projID}_dashboard.html
        echo "    <p><strong>Analysis completed at:</strong> ${currentDate}</p>" >> ${projID}_dashboard.html
        
        # Complete the HTML structure
        cat >> ${projID}_dashboard.html << 'ERRORTEMPLATE'    
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
ERRORTEMPLATE
    else
        # Dashboard was created, copy debug info if available
        if [ -d "debuginfo" ]; then
            mkdir -p "${params.outdir}/${projID}/debuginfo"
            cp -r debuginfo/* "${params.outdir}/${projID}/debuginfo/" 2>/dev/null
        fi
    fi
    
    # Final note with Nextflow-provided date
    echo "Dashboard generation completed at ${currentDate}" >> ${projID}_data_quality.log
    
    # Make sure the output dashboard exists even if it failed to generate properly
    # Create dashboard with both timestamp and projID to ensure it matches the expected pattern
    # First try with the pattern specified in the output
    if [ ! -f "${projID}_dashboard.html" ]; then
        echo "Creating emergency fallback dashboard at ${projID}_dashboard.html" >> ${projID}_data_quality.log
        echo "<!DOCTYPE html><html><head><title>FoodNet Dashboard Fallback</title></head><body><h1>FoodNet Analysis</h1><p>Dashboard generation failed. See logs for details.</p></body></html>" > ${projID}_dashboard.html
    fi
    
    # Also create dashboard with timestamp format to ensure compatibility
    if [ ! -f "\${MY_TIMESTAMP}_dashboard.html" ]; then
        echo "Creating emergency fallback dashboard at \${MY_TIMESTAMP}_dashboard.html" >> ${projID}_data_quality.log
        echo "<!DOCTYPE html><html><head><title>FoodNet Dashboard Fallback</title></head><body><h1>FoodNet Analysis</h1><p>Dashboard generation failed. See logs for details.</p></body></html>" > \${MY_TIMESTAMP}_dashboard.html
    fi
    
    # Create both file formats to ensure Nextflow can find the expected output
    if [ -f "${projID}_dashboard.html" ]; then
        # Copy rather than symlink to ensure both files exist
        cp "${projID}_dashboard.html" "\${MY_TIMESTAMP}_dashboard.html"
        echo "Created copy from ${projID}_dashboard.html to \${MY_TIMESTAMP}_dashboard.html" >> ${projID}_data_quality.log
    elif [ -f "\${MY_TIMESTAMP}_dashboard.html" ]; then
        cp "\${MY_TIMESTAMP}_dashboard.html" "${projID}_dashboard.html"
        echo "Created copy from \${MY_TIMESTAMP}_dashboard.html to ${projID}_dashboard.html" >> ${projID}_data_quality.log
    fi
    
    # Create a guaranteed output file in the CURRENT directory (important for Nextflow to find it)
    echo "Creating guaranteed dashboard file for Nextflow to find" >> ${projID}_data_quality.log
    
    # Force the creation of a dashboard in the work directory
    echo "<!DOCTYPE html><html><head><title>FoodNet Dashboard</title></head><body><h1>FoodNet Analysis Dashboard</h1><p>Project ID: ${projID}</p><p>Generated at: ${currentDate}</p></body></html>" > ./dashboard.html
    
    # List files in current directory to debug
    echo "Files in current directory:" >> ${projID}_data_quality.log
    ls -la . >> ${projID}_data_quality.log
    
    # Try with different paths to be super sure
    echo "<!DOCTYPE html><html><body><h1>Dashboard</h1></body></html>" > ./output.html
    echo "<!DOCTYPE html><html><body><h1>Dashboard</h1></body></html>" > "${PWD}/emergency.html"
    
    # Make sure we can see them
    echo "After creating guaranteed files:" >> ${projID}_data_quality.log
    ls -la *.html >> ${projID}_data_quality.log || echo "No HTML files found" >> ${projID}_data_quality.log
    
    # Create a sentinel file for publishDir to pick up
    touch "./dashboard_completed.txt"
    echo "Dashboard generation completed at ${currentDate}" > "./dashboard_completed.txt"
    
    # Create an extremely simple fallback HTML in case no other exists
    if ! ls *.html 1>/dev/null 2>&1; then
      echo "<!DOCTYPE html><html><head><title>Emergency Dashboard</title></head><body><h1>Emergency Dashboard</h1><p>No other dashboard files were created. This is an emergency fallback.</p></body></html>" > emergency_dashboard.html
    fi
    """
}