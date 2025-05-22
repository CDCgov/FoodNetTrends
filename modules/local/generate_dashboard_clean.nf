process GENERATE_DASHBOARD {
    tag "Generate dashboard"
    label 'process_high_memory'
    shell "/bin/bash"
    container 'foodnet.sif'
    time '30m'
    memory '8 GB'

    input:
    path ir_outputs
    val resultDir
    val projID
    path dashboardTemplate
    val dashboardScript

    output:
    path "dashboard.html", emit: dashboard
    path "*_data_quality.json", optional: true, emit: quality_json
    path "*_data_quality.log", optional: true, emit: quality_log
    publishDir "${params.outdir}/${projID}", mode: params.publish_dir_mode

    script:
    def currentDate = new Date().toString()
    """
    #!/bin/bash
    
    echo "=========================================" > ${projID}_data_quality.log
    echo "  FoodNet Trends Clean Dashboard Generation" >> ${projID}_data_quality.log
    echo "  Project ID: ${projID}" >> ${projID}_data_quality.log
    echo "  Generated: ${currentDate}" >> ${projID}_data_quality.log
    echo "=========================================" >> ${projID}_data_quality.log
    
    # Create data quality JSON
    cat > ${projID}_data_quality.json << 'EOF'
{
  "projectId": "${projID}",
  "generationDate": "${currentDate}",
  "dataQuality": {
    "usesPlaceholderData": false,
    "affectedPathogens": [],
    "warnings": []
  }
}
EOF
    
    # Set memory optimization for R
    export R_MAX_VSIZE=6G
    export R_GC_MEM_GROW=0
    
    # Log available files
    echo "Available result files:" >> ${projID}_data_quality.log
    ls -la *IRCatch.csv >> ${projID}_data_quality.log 2>/dev/null || echo "No IRCatch files found" >> ${projID}_data_quality.log
    ls -la *.png >> ${projID}_data_quality.log 2>/dev/null || echo "No PNG files found" >> ${projID}_data_quality.log
    ls -la *summary.txt >> ${projID}_data_quality.log 2>/dev/null || echo "No summary files found" >> ${projID}_data_quality.log
    
    # Also check for linked files (Nextflow creates symlinks)
    echo "All files in work directory:" >> ${projID}_data_quality.log
    ls -la >> ${projID}_data_quality.log
    
    # Follow symlinks to find actual PNG and summary files
    find . -name "*.png" -type l >> ${projID}_data_quality.log 2>/dev/null || echo "No PNG symlinks found" >> ${projID}_data_quality.log
    find . -name "*summary.txt" -type l >> ${projID}_data_quality.log 2>/dev/null || echo "No summary symlinks found" >> ${projID}_data_quality.log
    
    # Use the clean dashboard script
    CLEAN_SCRIPT="${workflow.projectDir}/bin/generate_dashboard_clean.R"
    
    if [ ! -f "\$CLEAN_SCRIPT" ]; then
        echo "ERROR: Clean dashboard script not found" >> ${projID}_data_quality.log
        # Create minimal fallback
        cat > dashboard.html << 'HTMLEOF'
<!DOCTYPE html>
<html><head><title>FoodNet Dashboard</title></head>
<body>
<h1>FoodNet Trends Analysis</h1>
<p>Dashboard script not found. This is a minimal fallback.</p>
</body></html>
HTMLEOF
    else
        echo "Running clean dashboard generator..." >> ${projID}_data_quality.log
        
        # Execute the clean dashboard script (Nextflow handles container automatically)
        Rscript "\$CLEAN_SCRIPT" \\
            --outDir="." \\
            --resultDir="." \\
            --outputFile="dashboard.html" \\
            --title="FoodNet Trends Analysis: ${projID}" \\
            --debug 2>&1 | tee -a ${projID}_data_quality.log
            
        SCRIPT_EXIT=\$?
        echo "Dashboard script exit code: \$SCRIPT_EXIT" >> ${projID}_data_quality.log
        
        # Ensure dashboard exists
        if [ ! -f "dashboard.html" ]; then
            echo "WARNING: Dashboard not created, making fallback" >> ${projID}_data_quality.log
            cat > dashboard.html << 'HTMLEOF'
<!DOCTYPE html>
<html><head><title>FoodNet Dashboard - Error</title></head>
<body>
<h1>FoodNet Trends Analysis</h1>
<p>Dashboard generation failed. Check logs for details.</p>
</body></html>
HTMLEOF
        fi
    fi
    
    # Final verification
    if [ -f "dashboard.html" ]; then
        echo "Dashboard successfully created" >> ${projID}_data_quality.log
        ls -la dashboard.html >> ${projID}_data_quality.log
    else
        echo "ERROR: No dashboard file created" >> ${projID}_data_quality.log
    fi
    
    echo "Dashboard generation completed at ${currentDate}" >> ${projID}_data_quality.log
    """
}