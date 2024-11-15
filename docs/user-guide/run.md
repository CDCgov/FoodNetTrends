## Running the Pipeline

The pipeline is implemented using Nextflow and can be run in two ways:

### Required Parameters

When running the pipeline, the following parameters are required:

- `--outdir`: Absolute path to the output directory.
  - **Example**: `/path/to/output/results`
- `--mmwrFile`: Absolute path to the MMWR data file.
  - **Example**: `/path/to/data/mmwr_data.sas7bdat`
- `--censusFile_B`: Absolute path to the bacterial census file.
  - **Example**: `/path/to/data/census_bacteria.sas7bdat`
- `--censusFile_P`: Absolute path to the parasitic census file.
  - **Example**: `/path/to/data/census_parasite.sas7bdat`
- `--travel`: List of travel statuses to include (comma-separated).
  - **Example**: `"NO,UNKNOWN,YES"`
- `--cidt`: List of CIDT variables (comma-separated).
  - **Example**: `"CIDT+,CX+,PARASITIC"`
- `--projID`: Unique project identifier.
  - **Example**: `"20240706"`

### Method 1: Running Directly with Nextflow

You can run the pipeline directly by invoking Nextflow with the required parameters:

```bash
module load nextflow singularity conda

nextflow run main.nf \
    -entry CDC_SPLINE \
    -profile singularity,conda \
    -with-conda \
    -work-dir /path/to/output/work \
    --outdir /path/to/output/results \
    --mmwrFile "/path/to/data/mmwr_data.sas7bdat" \
    --censusFile_B "/path/to/data/census_bacteria.sas7bdat" \
    --censusFile_P "/path/to/data/census_parasite.sas7bdat" \
    --travel "NO,UNKNOWN,YES" \
    --cidt "CIDT+,CX+,PARASITIC" \
    --projID "20240706"
```

**Explanation of the command:**

- **Module Loading**: Ensure that `nextflow`, `singularity`, and `conda` are loaded in your environment.
- **`-entry CDC_SPLINE`**: Specifies the entry workflow to run (defined in `main.nf`).
- **`-profile singularity,conda`**: Uses the Singularity container and Conda environment profiles.
- **`-with-conda`**: Enables the use of Conda environments specified in the pipeline.
- **`-work-dir`**: Specifies the working directory for Nextflow.
- **Parameter Flags (`--`)**: Provide the necessary parameters as described above.

### Method 2: Using the `run_workflow.sh` Script

Alternatively, you can use the provided `run_workflow.sh` script to execute the pipeline.

1. **Update the Script:**

   - Open the `run_workflow.sh` script.
   - Update the `outDir` and `dataDir` variables with the appropriate paths.
   - Ensure that the script includes the required parameters (`mmwrFile`, `censusFile_B`, `censusFile_P`, `travel`, `cidt`, `projID`).

2. **Run the Script:**

   ```bash
   bash run_workflow.sh run
   ```

   The script will set up the environment and execute the pipeline with the specified parameters.

## Interpreting Output

After the pipeline completes, you'll find several files and directories in your output folder (`--outdir`). These include:

- **SplineResults/**: Contains the results of the spline modeling, including `.Rds` files and plots (`.png` files).
- **EstIRRCatch_summary.csv**: A summary CSV file combining the estimation results from the spline models.
- **Logs**: Detailed logs of the pipeline execution for troubleshooting.

### Output Files and Directories

- `SplineResults/`: Directory containing:

  - `*.Rds`: R data files resulting from the spline modeling.
  - `*.png`: Plots generated from the modeling.

- `EstIRRCatch_summary.csv`: A combined CSV file summarizing the estimation of Incidence Rate Ratios (IRR) by catchment area.
