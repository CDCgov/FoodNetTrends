#!/bin/bash
source /etc/profile
module load miniconda3

#$ -cwd
#$ -pe smp 12-32
#$ -N TEST
#$ -M qok9@cdc.gov
#$ -m e
#$ -q all.q
#$ -l h_vmem=256G

conda activate brms

R CMD BATCH /scicomp/home-pure/qok9/FN1_24_CxCIDT/Shigella/SplinesModel/Unspeciated_Shigella.R  /scicomp/home-pure/qok9/FN1_24_CxCIDT/Rout/Unspeciated_Shigella.Rout

conda deactivate
echo "Script ran to completion"
