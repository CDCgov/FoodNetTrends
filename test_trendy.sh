#!/bin/bash
# Test script to verify trendy.nf fixes

# Create minimal test data
echo "Creating test data..."
mkdir -p test_data

# Create minimal MMWR data
echo "pathogen,state,year,count,population" > test_data/test_mmwr.csv
echo "CYCLOSPORA,CA,2020,5,10000000" >> test_data/test_mmwr.csv
echo "CYCLOSPORA,NY,2020,3,15000000" >> test_data/test_mmwr.csv
echo "CYCLOSPORA,GA,2021,4,8000000" >> test_data/test_mmwr.csv
echo "CYCLOSPORA,MD,2021,2,6000000" >> test_data/test_mmwr.csv

# Create minimal parasitic census data
echo "state,population,year,pathogentype" > test_data/test_census_para.csv
echo "CA,10000000,2020,Parasitic" >> test_data/test_census_para.csv
echo "NY,15000000,2020,Parasitic" >> test_data/test_census_para.csv
echo "GA,8000000,2021,Parasitic" >> test_data/test_census_para.csv
echo "MD,6000000,2021,Parasitic" >> test_data/test_census_para.csv

# Create minimal bacterial census data (not needed but included for completeness)
echo "state,population,year,pathogentype" > test_data/test_census_bact.csv
echo "CA,10000000,2020,Bacterial" >> test_data/test_census_bact.csv
echo "NY,15000000,2020,Bacterial" >> test_data/test_census_bact.csv
echo "GA,8000000,2021,Bacterial" >> test_data/test_census_bact.csv
echo "MD,6000000,2021,Bacterial" >> test_data/test_census_bact.csv

# Create test output directory
mkdir -p test_output

# Run a simple Nextflow process to test the module directly
echo "Running test process..."

# Create temporary workflow file
cat > test_workflow.nf << 'EOF'
#!/usr/bin/env nextflow

// Define process inputs
params.mmwrFile = "test_data/test_mmwr.csv"
params.censusFileB = "test_data/test_census_bact.csv"
params.censusFileP = "test_data/test_census_para.csv"
params.pathogen = "CYCLOSPORA"
params.projID = "test_run"
params.outdir = "test_output"
params.filter_travel = "NO,UNKNOWN,YES"
params.filter_cidt = "CIDT+,CX+,PARASITIC"
params.cores = 1
params.chains = 1
params.iterations = 10
params.seed = 123
params.publish_dir_mode = "copy"
params.output_suffix = "brm"

// Import the trendy process
include { TRENDY } from './modules/local/trendy'

// Create channels from inputs
mmwr_ch = Channel.fromPath(params.mmwrFile)
census_b_ch = Channel.fromPath(params.censusFileB)
census_p_ch = Channel.fromPath(params.censusFileP)

// Create combined channel
process_ch = Channel.of(params.pathogen)
    .combine(mmwr_ch)

// Create workflow
workflow {
    TRENDY(
        process_ch,
        census_b_ch,
        census_p_ch,
        params.projID,
        "${projectDir}/bin",
        params.filter_travel,
        params.filter_cidt,
        "${projectDir}/bin/templates"
    )
}
EOF

# Run the test workflow
echo "Starting Nextflow process..."
nextflow run test_workflow.nf -profile singularity

# Check for expected output file
if [ -f "test_output/test_run/CYCLOSPORA_brm.Rds" ]; then
    echo "SUCCESS: Test process completed and created expected output file"
    exit 0
else
    echo "ERROR: Test failed to create expected output file"
    # Check for error logs
    if [ -f "test_output/CYCLOSPORA_trendy.log" ]; then
        echo "Log file contents:"
        cat test_output/CYCLOSPORA_trendy.log
    fi
    exit 1
fi 