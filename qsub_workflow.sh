#!/bin/bash
#$ -N FoodNetTrends
#$ -cwd
#$ -V
#$ -pe smp 8
#$ -l h_vmem=32G
#$ -l h_rt=01:00:00
#$ -S /bin/bash
#$ -o foodnettrends_hpc_run.out
#$ -e foodnettrends_hpc_run.err

# Set up paths, files
outDir="/scicomp/home-pure/smn9/projects/ticket/RSpline_test"
dataDir="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/"

# Load modules
module load nextflow/24.04.2
module load singularity
module load conda/24.3.0

# Housekeeping
if [[ ! -d $outDir ]]; then mkdir -p $outDir; fi
flag=$1

# Run the pipeline with different arguments depending on the flag
if [[ $flag == "run" ]]; then
    nextflow run main.nf \
    -entry CDC_SPLINE \
    -profile singularity,conda \
    -with-conda \
    -work-dir $outDir/work \
    --outdir $outDir \
    --mmwrFile "$dataDir/mmwr9623_Jan2024.sas7bdat" \
    --censusFile_B "$dataDir/cen9623.sas7bdat" \
    --censusFile_P "$dataDir/cen9623_para.sas7bdat" \
    --travel "NO,UNKNOWN,YES" \
    --cidt "CIDT+,CX+,PARASITIC" \
    --projID "20241114"
fi

if [[ $flag == "full" ]]; then
    nextflow run main.nf \
    -entry CDC_SPLINE \
    -profile singularity,conda \
    -with-conda \
    -work-dir $outDir/work \
    --outdir $outDir \
    --mmwrFile "$dataDir/mmwr9623.sas7bdat" \
    --censusFile_B "$dataDir/cen9623.sas7bdat" \
    --censusFile_P "$dataDir/cen9623_para.sas7bdat" \
    --travel "NO,UNKNOWN,YES" \
    --cidt "CIDT+,CX+,PARASITIC" \
    --projID "20241114"
fi

