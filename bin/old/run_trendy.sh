#!/usr/bin/env bash

singularity run ~/FoodNetTrends/foodnet.sif \
  Rscript $(pwd)/new_trendy.R \
  --mmwrFile "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/mmwr9623_Jan2024.sas7bdat" \
  --censusFile_B "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623.sas7bdat" \
  --censusFile_P "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623_para.sas7bdat" \
  --travel "NO,UNKNOWN,YES" \
  --cidt "CIDT+,CX+,PARASITIC" \
  --projID "20250204" \
  --debug FALSE

