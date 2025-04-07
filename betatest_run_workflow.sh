#!/bin/bash

# Please set up paths and files

# Output directory (replace with your desired output path)
outDir="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/mmwroutput"
# e.g., outDir="/home/user/projects/my_project/output"

# Data directory containing input files (replace with your data directory path)
dataDir="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/data"
# e.g., dataDir="/home/user/data"

# Load necessary modules
module purge
module load nextflow/24.04.2
module load singularity
module load conda/24.3.0

# Housekeeping
if [[ ! -d $outDir ]]; then mkdir -p $outDir; fi
flag=$1

if [[ $flag == "run" ]]; then
    nextflow run main.nf \
        -entry CDC_SPLINE \
        -profile singularity,conda \
        -with-conda \
        -work-dir $outDir/work \
        --outdir $outDir \
        --mmwrFile "$dataDir/mmwr9623.sas7bdat" \
        # e.g., --mmwrFile "$dataDir/mmwr_sample_data.sas7bdat" \
        --censusFile_B "$dataDir/cen9623.sas7bdat" \
        # e.g., --censusFile_B "$dataDir/census_bacteria.sas7bdat" \
        --censusFile_P "$dataDir/cen9623_para.sas7bdat" \
        # e.g., --censusFile_P "$dataDir/census_parasite.sas7bdat" \
        --travel "NO, UNKNOWN,YES" \
        # e.g., --travel "NO,UNKNOWN,YES" \
        --cidt "CIDT+,CX+,PARASITIC" \
        # e.g., --cidt "CIDT+,CX+,PARASITIC" \
        --projID "BETA_2025APR7"
        # e.g., --projID "20240706"
fi

if [[ $flag == "full" ]]; then
    nextflow run main.nf \
        -entry CDC_SPLINE \
        -profile singularity,conda \
        -with-conda \
        -work-dir $outDir/work \
        --outdir $outDir \
        --mmwrFile "$dataDir/mmwr_sample_data.sas7bdat" \
        # e.g., --mmwrFile "$dataDir/mmwr_full_data.sas7bdat" \
        --censusFile_B "$dataDir/cen9623.sas7bdat" \
        # e.g., --censusFile_B "$dataDir/census_bacteria_full.sas7bdat" \
        --censusFile_P "$dataDir/cen9623_para.sas7bdat" \
        # e.g., --censusFile_P "$dataDir/census_parasite_full.sas7bdat" \
        --travel "NO,UNKNOWN,YES" \
        # e.g., --travel "NO,UNKNOWN,YES" \
        --cidt "CIDT+,CX+,PARASITIC" \
        # e.g., --cidt "CIDT+,CX+,PARASITIC" \
        --projID "BETA_2025APR7"
        # e.g., --projID "20240830"
fi
