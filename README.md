# Background

The Foodborne Diseases Active Surveillance Network (FoodNet) monitors illnesses caused by enteric and foodborne pathogens across 10 U.S. sites. FoodNet data is used to track trends in these illnesses and to monitor progress toward federal disease reduction goals.

The original model for analyzing FoodNet data faced limitations, such as sensitivity to single-year aberrations and biases toward more populous sites. To address these issues, this enhanced model (FoodNetTrends) was developed using a Bayesian framework, incorporating thin-plate splines and site-specific interactions.

Key improvements include:
- Treating the year as a continuous variable.
- Including site-specific trends.
- Improved ability to handle uncertainty and noisy data.

# User Guide

## Table of Contents

- [Getting Started](#getting-started)
- [Preparing Files](#preparing-files)
- [Running the Pipeline](#running-the-pipeline)
  - [Required Parameters](#required-parameters)
  - [Method 1: Running Directly with Nextflow](#method-1-running-directly-with-nextflow)
  - [Method 2: Using the `run_workflow.sh` Script](#method-2-using-the-run_workflowsh-script)
- [Interpreting Output](#interpreting-output)
- [Advanced Configuration and Optimization](#advanced-configuration-and-optimization)

## Getting Started

To use the FoodNet enhanced model pipeline, you will need:

- **R (version 4.3.2)** with the `brms` and `tidybayes` packages.
- **Nextflow (version 24.04.2)** for running the pipeline.
- Access to the **FoodNet surveillance dataset**.
- **Singularity** or **Docker** for containerization (Singularity is used in the examples).

For a detailed list of software and dependencies, please refer to the [GitHub repository](https://github.com/CDCgov/FoodNetTrends).

We highly recommend using Docker or Singularity containers for full pipeline reproducibility. If these are not possible, Conda is also supported. See the [`-profile`](#-profile) section for more information.

## Preparing Files

Data preparation is crucial for accurate modeling. The FoodNet enhanced model requires the following inputs:

- **MMWR Data File** (`mmwrFile`): The MMWR (Morbidity and Mortality Weekly Report) data file in SAS format (`.sas7bdat`).
- **Census Data Files**: Two census data files in SAS format:
  - **Bacterial Census File** (`censusFile_B`)
  - **Parasitic Census File** (`censusFile_P`)
- **Parameters**:
  - **Travel Status** (`travel`): A list of travel statuses to include (e.g., `"NO,UNKNOWN,YES"`).
  - **CIDT Variables** (`cidt`): A list of CIDT (Culture-Independent Diagnostic Tests) variables (e.g., `"CIDT+,CX+,PARASITIC"`).
  - **Project ID** (`projID`): A unique project identifier (e.g., `"20240706"`).

Ensure that the data files are accessible and properly formatted.

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

### Understanding the Results

- **Spline Modeling Results**: The `.Rds` files can be loaded into R for further analysis or visualization.
- **Plots**: The `.png` files provide visual representations of the spline models, trends, and other relevant analyses.
- **Summary CSV**: `EstIRRCatch_summary.csv` contains aggregated results, which can be opened with any spreadsheet software or analyzed programmatically.

## Advanced Configuration and Optimization

### Core Nextflow Arguments

#### `-profile`

Use this parameter to choose a configuration profile. Profiles provide presets for different compute environments.

**Available Profiles:**

- `test`: Configuration for automated testing (if test data is available).
- `docker`: Uses [Docker](https://docker.com/).
- `singularity`: Uses [Singularity](https://sylabs.io/docs/).
- `conda`: Uses [Conda](https://conda.io/docs/).

**Note:**

- We recommend using Docker or Singularity for reproducibility.
- Multiple profiles can be loaded, e.g., `-profile singularity,conda`.
- If `-profile` is not specified, the pipeline runs locally, which is **not** recommended.

#### `-resume`

Resume a previous pipeline run:

```bash
nextflow run main.nf -resume
```

#### `-c`

Specify a custom Nextflow config file (for resource specifications or infrastructural tweaks):

```bash
nextflow run main.nf -c /path/to/custom.config
```

**Warning:** Do not use `-c <file>` to specify pipeline parameters. Use it only for resource configurations.

### Custom Configuration

#### Resource Requests

Customize compute resources by adjusting the Nextflow configuration. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for details.

#### Custom Containers

To use different containers or Conda environments for specific tools, adjust the profiles or configuration files accordingly.

#### Custom Tool Arguments

If you need to provide additional arguments to the R scripts or other tools within the pipeline, you may need to modify the pipeline scripts (`main.nf`, `spline.nf`, `trendy.nf`) accordingly.

### Running in the Background

Run Nextflow in the background:

```bash
nextflow run main.nf ... -bg
```

Alternatively, use `screen`, `tmux`, or submit Nextflow as a job to your scheduler.

### Nextflow Memory Requirements

Limit Nextflow's Java virtual machine memory usage by adding to your environment:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

---
## Credits

The Spline pipeline was largely developed by Samantha Sevilla and OAMD's SciComp Team with support Daniel Weller (CDC/DFWED/EDEB), based on R scripts developed by Daniel Weller (CDC/DFWED/EDEB) with support from Beau Bruce (CDC/DFWED/EDEB) and Erica Billig Rose (CDC/DFWED/EDEB). Detailed contributions can be found in our user-guides. 

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
