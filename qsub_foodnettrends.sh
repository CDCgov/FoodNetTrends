#!/bin/bash
source /etc/profile

#$ -N FoodNetTrends
#$ -q all.q
#$ -pe smp 32
#$ -l h_vmem=128G
#$ -cwd

module load nextflow/24.10.4
module load singularity/4.1.4
module load java/17.0.6

export NCPUS=$NSLOTS

nextflow run main.nf -profile singularity -entry FoodNetTrends \
  --mmwrFile /scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/mmwr9623_Jan2024.sas7bdat \
  --censusFile_B /scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623.sas7bdat \
  --censusFile_P /scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623_para.sas7bdat

