#!/bin/bash

# set up paths, files
dataDir="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/"
outDir="output"  # Default to "output" directory in current location

# set up modules
module purge
module load nextflow/24.10.4
module load singularity/4.1.4
module load java/17.0.6

# Parse command line arguments
flag=$1
background=false

# Check for additional flags
for arg in "$@"; do
  if [[ "$arg" == "--bg" ]]; then
    background=true
  fi
  if [[ "$arg" == "--outdir="* ]]; then
    outDir="${arg#*=}"
  fi
done

# Create timestamp for automatic project ID
timestamp=$(date +%Y%m%d_%H%M%S)

# Build the base command
cmd="nextflow run main.nf -profile singularity -entry SPLINE \
  --mmwrFile \"$dataDir/mmwr9623_Jan2024.sas7bdat\" \
  --censusFileB \"$dataDir/cen9623.sas7bdat\" \
  --censusFileP \"$dataDir/cen9623_para.sas7bdat\" \
  --iterations 500 \
  --chains 2 \
  --adapt_delta 0.95 \
  --max_treedepth 10 \
  --seed 123 \
  --outdir \"$outDir\""

# Add background option if requested
if [[ $background == true ]]; then
  cmd="nextflow -bg $cmd > foodnet_run_${timestamp}.log"
fi

# Run with different options based on flag
if [[ $flag == "test" || $flag == "" ]]; then
  # Test run with minimal pathogens (default)
  # Use reduced parameters for faster testing
  test_cmd="${cmd/--iterations 500/--iterations 100}"
  test_cmd="${test_cmd/--chains 2/--chains 1}"
  test_cmd="${test_cmd/--adapt_delta 0.95/--adapt_delta 0.8}"
  test_cmd="${test_cmd/--max_treedepth 10/--max_treedepth 8}"
  eval $test_cmd
elif [[ $flag == "full" ]]; then
  # Full run with production parameters
  eval $cmd
elif [[ $flag == "resume" ]]; then
  # Resume a previous run
  eval $cmd" -resume"
else
  # Show usage
  echo "Usage: ./run_workflow.sh [test|full|resume] [--bg] [--outdir=path]"
  echo ""
  echo "Options:"
  echo "  test    Run with minimal parameters for testing (default)"
  echo "  full    Run with production parameters"
  echo "  resume  Resume previous run"
  echo "  --bg    Run in background"
  echo "  --outdir=path  Specify output directory (default: ./output)"
  exit 1
fi
