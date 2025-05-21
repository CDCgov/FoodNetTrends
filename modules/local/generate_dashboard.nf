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
    echo "  Generated: $(date)" >> ${projID}_data_quality.log
    echo "=========================================" >> ${projID}_data_quality.log
    echo "" >> ${projID}_data_quality.log
    
    # Initialize JSON structure for data quality metadata
    cat > ${projID}_data_quality.json << EOF
    {
      "projectId": "${projID}",
      "generationDate": "$(date -Iseconds)",
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
    
    # Look for placeholder warnings across all files (txt, log, summary files)
    placeholder_files=\$(grep -l "PLACEHOLDER DATA\\|placeholder data\\|SYNTHETIC DATA" *.{log,txt,summary} 2>/dev/null || echo "")
    data_quality_warnings=0
    
    # Check data summary files for detailed pathogen-specific warnings
    pathogen_warnings=""
    affected_pathogens=""
    for summary_file in *_data_summary.txt 2>/dev/null; do
        if [ -f "\$summary_file" ]; then
            pathogen=\$(echo "\$summary_file" | sed 's/_data_summary.txt//')
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
    done
    
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
        
        # Create a detailed warning banner for the dashboard with pathogen-specific info
        if [ -n "\$pathogen_warnings" ]; then
            warning_details="<ul style='margin-top:10px;text-align:left;'>\$(echo -e "\$pathogen_warnings" | sed 's/- /<li>/g' | sed 's/$/<\\/li>/g')</ul>"
        else
            warning_details=""
        fi
        
        warning_banner="<div style='background-color: #f8d7da; color: #721c24; padding: 15px; margin: 20px 0; border: 1px solid #f5c6cb; border-radius: 5px;'><h3 style='margin-top:0'>⚠️ WARNING: Placeholder Data Detected</h3><p>This analysis contains placeholder data which may not represent real-world conditions.</p><p><strong>Results are NOT suitable for production or public health decision-making!</strong></p>\$warning_details</div>"
        
        # Create a modified template with the warning banner
        if [ -f "${dashboardTemplate}" ]; then
            cp "${dashboardTemplate}" modified_template.html
            # Insert warning after opening body tag
            sed -i "s/<body>/<body>\\n\$warning_banner/" modified_template.html
            dashboard_template="modified_template.html"
        else
            echo "ERROR: Could not find dashboard template: ${dashboardTemplate}" >> ${projID}_data_quality.log
            dashboard_template="${dashboardTemplate}"
        fi
    else
        echo "No placeholder data warnings detected." >> ${projID}_data_quality.log
        dashboard_template="${dashboardTemplate}"
    fi
    
    # Check for data consistency
    echo "" >> ${projID}_data_quality.log
    echo "Checking data consistency..." >> ${projID}_data_quality.log
    
    # Count IR files
    ir_files=\$(ls -1 *_IRCatch.csv 2>/dev/null | wc -l)
    echo "Found \$ir_files incidence rate files." >> ${projID}_data_quality.log
    
    # Add to JSON
    sed -i "s/\"dataConsistency\": {}/\"dataConsistency\": {\\n      \\\"incidenceRateFiles\\\": \$ir_files,\\n      \\\"dataQualityWarnings\\\": \$data_quality_warnings\\n    }/" ${projID}_data_quality.json
    
    # Generate the dashboard
    echo "" >> ${projID}_data_quality.log
    echo "Generating dashboard with template: \$dashboard_template" >> ${projID}_data_quality.log
    
    Rscript \\
      ${dashboardScript} \\
      --outDir=${resultDir} \\
      --resultDir=${resultDir} \\
      --outputFile=${projID}_dashboard.html \\
      --title="FoodNet Trends Analysis: ${projID}" \\
      --templateFile=\$dashboard_template \\
      --qualityDataPath=${projID}_data_quality.json
    
    # Add note to quality log
    echo "" >> ${projID}_data_quality.log
    echo "Dashboard generation completed at $(date)" >> ${projID}_data_quality.log
    """
}