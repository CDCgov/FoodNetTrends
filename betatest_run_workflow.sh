#!/bin/bash

outDir="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/mmwroutput"

dataDir="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data"

module purge
module load nextflow/24.04.2
module load singularity
module load miniconda

if [[ ! -d $outDir ]]; then mkdir -p $outDir; fi
flag=$1

if [[ $flag == "run" ]]; then
    nextflow run main.nf \
        -entry SPLINE \
        -profile singularity,conda \
        -with-conda \
        -work-dir $outDir/work \
        --outdir $outDir \
        --mmwrFile "$dataDir/mmwr9623.sas7bdat" \
        --censusFile_B "$dataDir/cen9623.sas7bdat" \
        --censusFile_P "$dataDir/cen9623_para.sas7bdat" \
        --travel "NO, UNKNOWN,YES" \
        --cidt "CIDT+,CX+,PARASITIC" \
        --projID "BETA_2025APR7"
fi

if [[ $flag == "full" ]]; then
    nextflow run main.nf \
        -entry SPLINE \
        -profile singularity,conda \
        -with-conda \
        -work-dir $outDir/work \
        --outdir $outDir \
        --mmwrFile "$dataDir/mmwr_sample_data.sas7bdat" \
        --censusFile_B "$dataDir/cen9623.sas7bdat" \
        --censusFile_P "$dataDir/cen9623_para.sas7bdat" \
        --travel "NO,UNKNOWN,YES" \
        --cidt "CIDT+,CX+,PARASITIC" \
        --projID "BETA_2025APR7"
fi
