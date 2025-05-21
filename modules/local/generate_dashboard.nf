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
    echo "  Generated: \\\$(date)" >> ${projID}_data_quality.log
    echo "=========================================" >> ${projID}_data_quality.log
    echo "" >> ${projID}_data_quality.log
    
    # Initialize JSON structure for data quality metadata
    cat > ${projID}_data_quality.json << EOF
    {
      "projectId": "${projID}",
      "generationDate": "\\\$(date -Iseconds)",
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
    
    # Look for placeholder warnings across all files (txt, log, summary files) - using safer find method
    placeholder_files=""
    for ext in log txt summary; do
        find . -maxdepth 1 -name "*.$ext" -type f -exec grep -l "PLACEHOLDER DATA\\|placeholder data\\|SYNTHETIC DATA" {} \\; >> placeholder_files.txt 2>/dev/null || true
    done
    if [ -s placeholder_files.txt ]; then
        placeholder_files=\\\$(cat placeholder_files.txt)
    fi
    data_quality_warnings=0
    
    # Check data summary files for detailed pathogen-specific warnings
    pathogen_warnings=""
    affected_pathogens=""
    # Use safer method to find summary files 
    find . -maxdepth 1 -name "*_data_summary.txt" > summary_files.txt
    while IFS= read -r summary_file; do
        if [ -f "\$summary_file" ]; then
            pathogen=\\\$(echo "\$summary_file" | sed 's/_data_summary.txt//')
            if grep -q "PLACEHOLDER DATA\\|placeholder data\\|SYNTHETIC DATA" "\$summary_file"; then
                if [ -z "\$affected_pathogens" ]; then
                    affected_pathogens="\\\"\$pathogen\\\""
                else
                    affected_pathogens="\$affected_pathogens, \\\"\$pathogen\\\""
                fi
                echo "WARNING: Placeholder data detected for pathogen: \$pathogen" >> ${projID}_data_quality.log
                pathogen_warnings="\$pathogen_warnings\\n- \$pathogen: Uses placeholder data (results NOT suitable for production use)"
                data_quality_warnings=\$((data_quality_warnings + 1))
            fi
        fi
    done < summary_files.txt
    
    # Update the JSON with pathogen-specific warnings
    if [ -n "\$affected_pathogens" ]; then
        # Update the JSON to indicate placeholder data is being used
        sed -i 's/"usesPlaceholderData": false/"usesPlaceholderData": true/' ${projID}_data_quality.json
        # Add the affected pathogens to the JSON
        sed -i "s/\"affectedPathogens\": \\[\\]/\"affectedPathogens\": [\$affected_pathogens]/" ${projID}_data_quality.json
    fi
    
    # Process overall placeholder data warnings
    if [ -n "\$placeholder_files" ] || [ -n "\$pathogen_warnings" ]; then
        echo "WARNING: Placeholder data detected in the following files:" >> ${projID}_data_quality.log
        if [ -n "\$placeholder_files" ]; then
            echo "\$placeholder_files" | tr ' ' '\\n' >> ${projID}_data_quality.log
        fi
        echo "" >> ${projID}_data_quality.log
        echo "*** RESULTS MAY NOT BE SUITABLE FOR PRODUCTION USE ***" >> ${projID}_data_quality.log
        echo "" >> ${projID}_data_quality.log
        
        # Add pathogen-specific warnings if available
        if [ -n "\$pathogen_warnings" ]; then
            echo "Pathogen-specific data quality issues:" >> ${projID}_data_quality.log
            echo -e "\$pathogen_warnings" >> ${projID}_data_quality.log
            echo "" >> ${projID}_data_quality.log
        fi
        
        # Create a warning banner with simplified approach
        if [ -n "\$pathogen_warnings" ]; then
            # Create a temporary HTML file for warnings
            echo "<ul style='margin-top:10px;text-align:left;'>" > warnings.html
            # For each warning line, convert to HTML list item
            echo -e "\$pathogen_warnings" | while read line; do
                if [ -n "\$line" ]; then
                    # Remove leading dash if present
                    clean_line=\${line#- }
                    # Add as list item
                    echo "<li>\$clean_line</li>" >> warnings.html
                fi
            done
            echo "</ul>" >> warnings.html
            # Read the formatted HTML
            warning_details=\\\$(cat warnings.html)
        else
            warning_details=""
        fi
        
        warning_banner="<div style='background-color: #f8d7da; color: #721c24; padding: 15px; margin: 20px 0; border: 1px solid #f5c6cb; border-radius: 5px;'><h3 style='margin-top:0'>⚠️ WARNING: Placeholder Data Detected</h3><p>This analysis contains placeholder data which may not represent real-world conditions.</p><p><strong>Results are NOT suitable for production or public health decision-making!</strong></p>\$warning_details</div>"
        
        # Create a modified template with the warning banner - with safer file handling
        if [ -f "${dashboardTemplate}" ]; then
            echo "Using dashboard template: ${dashboardTemplate}" >> ${projID}_data_quality.log
            # Create a safe copy with proper error handling
            if ! cp "${dashboardTemplate}" modified_template.html; then
                echo "ERROR: Failed to copy dashboard template" >> ${projID}_data_quality.log
                # Create a minimal fallback template if copy fails
                cat > modified_template.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Fallback Template</title></head>
<body>
<h1>FoodNet Trends Analysis</h1>
<p>This is a fallback template due to template file access issues.</p>
EOF
            else
                # Insert warning after opening body tag
                if ! sed -i "s/<body>/<body>\\n\$warning_banner/" modified_template.html; then
                    echo "WARNING: Failed to insert warning banner, continuing with unmodified template" >> ${projID}_data_quality.log
                fi
            fi
            dashboard_template="modified_template.html"
        else
            echo "ERROR: Could not find dashboard template: ${dashboardTemplate}" >> ${projID}_data_quality.log
            # Create a minimal fallback template
            cat > fallback_template.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Fallback Template</title></head>
<body>
<h1>FoodNet Trends Analysis</h1>
<p>This is a fallback template due to missing template file.</p>
EOF
            dashboard_template="fallback_template.html"
        fi
    else
        echo "No placeholder data warnings detected." >> ${projID}_data_quality.log
        dashboard_template="${dashboardTemplate}"
    fi
    
    # Check for data consistency
    echo "" >> ${projID}_data_quality.log
    echo "Checking data consistency..." >> ${projID}_data_quality.log
    
    # Count IR files
    ir_files=\\\$(ls -1 *_IRCatch.csv 2>/dev/null | wc -l)
    echo "Found \$ir_files incidence rate files." >> ${projID}_data_quality.log
    
    # Add to JSON
    sed -i "s/\"dataConsistency\": {}/\"dataConsistency\": {\\n      \\\"incidenceRateFiles\\\": \$ir_files,\\n      \\\"dataQualityWarnings\\\": \$data_quality_warnings\\n    }/" ${projID}_data_quality.json
    
    # Generate the dashboard with improved error handling
    echo "" >> ${projID}_data_quality.log
    echo "Generating dashboard with template: \$dashboard_template" >> ${projID}_data_quality.log
    
    # Verify dashboard script exists
    if [ ! -f "${dashboardScript}" ]; then
        echo "ERROR: Dashboard script not found: ${dashboardScript}" >> ${projID}_data_quality.log
        # Create a simple HTML dashboard as fallback
        cat > ${projID}_dashboard.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Error</title></head>
<body>
<h1>FoodNet Trends Analysis Dashboard</h1>
<h2>Error: Script Not Found</h2>
<p>The dashboard generation script was not found. This is a simplified fallback dashboard.</p>
<p>Analysis completed at: \\\$(date)</p>
<p>Project ID: ${projID}</p>
</body>
</html>
EOF
        echo "Created fallback dashboard due to missing script" >> ${projID}_data_quality.log
    else
        # Run script with error handling
        echo "Executing dashboard script: ${dashboardScript}" >> ${projID}_data_quality.log
        
        # Set -e off temporarily to handle errors manually
        set +e
        Rscript \\
          "${dashboardScript}" \\
          --outDir="${resultDir}" \\
          --resultDir="${resultDir}" \\
          --outputFile="${projID}_dashboard.html" \\
          --title="FoodNet Trends Analysis: ${projID}" \\
          --templateFile="\$dashboard_template" \\
          --qualityDataPath="${projID}_data_quality.json"
          
        r_exit_code=\$?
        set -e
        
        # Handle script errors
        if [ \$r_exit_code -ne 0 ]; then
            echo "ERROR: Dashboard generation script failed with exit code \$r_exit_code" >> ${projID}_data_quality.log
            # Create a simple HTML dashboard as fallback
            cat > ${projID}_dashboard.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Error</title></head>
<body>
<h1>FoodNet Trends Analysis Dashboard</h1>
<h2>Error: Script Failed</h2>
<p>The dashboard generation script failed with exit code \$r_exit_code.</p>
<p>Analysis completed at: \\\$(date)</p>
<p>Project ID: ${projID}</p>
</body>
</html>
EOF
            echo "Created fallback dashboard due to script failure" >> ${projID}_data_quality.log
        else
            echo "Dashboard generation script executed successfully" >> ${projID}_data_quality.log
        fi
    fi
    
    # Verify dashboard was created
    if [ ! -f "${projID}_dashboard.html" ]; then
        echo "ERROR: Dashboard file was not created" >> ${projID}_data_quality.log
        # Create a simple HTML dashboard as last resort
        cat > ${projID}_dashboard.html << EOF
<!DOCTYPE html>
<html>
<head><title>FoodNet Trends Dashboard - Error</title></head>
<body>
<h1>FoodNet Trends Analysis Dashboard</h1>
<h2>Error: Missing Output</h2>
<p>The dashboard output file was not created properly.</p>
<p>Analysis completed at: \\\$(date)</p>
<p>Project ID: ${projID}</p>
</body>
</html>
EOF
        echo "Created last-resort dashboard due to missing output file" >> ${projID}_data_quality.log
    fi
    
    # Add note to quality log
    echo "" >> ${projID}_data_quality.log
    echo "Dashboard generation completed at \\\$(date)" >> ${projID}_data_quality.log
    """
}