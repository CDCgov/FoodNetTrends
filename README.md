# FoodNet Trends Analysis Pipeline

A Nextflow-based pipeline for Bayesian modeling of foodborne disease surveillance data.

## Overview

The FoodNet Trends pipeline processes food-borne illness surveillance data from the CDC's Foodborne Diseases Active Surveillance Network (FoodNet), applies Bayesian hierarchical models with splines, and generates standardized incidence rates, trend analysis, and visualizations. The pipeline is designed for epidemiologists and statisticians analyzing surveillance data for public health decision-making.

## Features

- **Data Preprocessing**: Cleans and standardizes raw MMWR SAS data files
- **Bayesian Modeling**: Implements hierarchical models with splines for flexible trend analysis
- **Multi-pathogen Support**: Analyzes all pathogen types present in the supplied surveillance data
- **Trend Visualization**: Generates site-specific and overall trend plots
- **Comparative Analysis**: Calculates relative risks and percent changes against historical baselines
- **Containerized Execution**: Ensures reproducible analysis using Singularity containers

## Pipeline Structure

```
FoodNetTrends/
├── bin/                  # R scripts for analysis
│   ├── calcIR.R          # Data preprocessing
│   ├── trendy.R          # Bayesian modeling
│   └── functions.R       # Reusable functions
├── modules/local/        # Nextflow process modules
│   ├── preprocess.nf     # Preprocessing module
│   └── trendy.nf         # Analysis module
├── workflows/            # Nextflow workflows
│   ├── preprocess.nf     # Preprocessing workflow
│   └── spline.nf         # Main analysis workflow
├── main.nf               # Nextflow entry point
├── foodnet.def           # Singularity container definition
├── foodnet.yml           # Conda environment specification
└── run_workflow.sh       # User-friendly execution script
```

## Requirements

- Nextflow ≥ 24.10.4
- Singularity ≥ 4.1.4
- Java ≥ 17
- Access to CDC HPC environment (recommended)

## Initial Setup

### Building the Singularity Container

Before running the pipeline, you need to build the Singularity container from the definition file:

1. **Ensure Singularity is installed**:
   ```bash
   singularity --version
   ```
   If not installed, follow the [Singularity installation guide](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html).

2. **Build the container from the definition file**:
   ```bash
   singularity build foodnet.sif foodnet.def
   ```
   This process may take several minutes as it downloads the base container and installs all required R packages.

3. **Verify the container was built correctly**:
   ```bash
   singularity exec foodnet.sif Rscript -e "library(brms); library(dplyr); print('Container test successful!')"
   ```
   You should see "Container test successful!" if everything is working correctly.

4. **Update container permissions if needed**:
   ```bash
   chmod 755 foodnet.sif
   ```

### Environment Configuration

1. **Load required modules** (HPC environment):
   ```bash
   module load nextflow/24.10.4
   module load singularity/4.1.4
   module load java/17.0.6
   ```

2. **Make the workflow script executable**:
   ```bash
   chmod +x run_workflow.sh
   ```

## Quick Start

### Interactive Mode

The easiest way to run the pipeline is using the interactive script:

```bash
./run_workflow.sh
```

This guides you through configuring the pipeline with prompts for:
- Operation mode (preprocess, full analysis, or use existing preprocessed data)
- Input file selection (MMWR data, census data)
- Filtering options (states, pathogens, travel status)
- Model parameters (chains, iterations, sampling details)

### Command Line Execution

For automated workflows or advanced usage:

```bash
module load nextflow/24.10.4 singularity/4.1.4 java/17.0.6

nextflow run main.nf \
  -profile singularity \
  -entry SPLINE \
  --mmwrFile "path/to/mmwr9623_Jan2024.sas7bdat" \
  --censusFileB "path/to/cen9623.sas7bdat" \
  --censusFileP "path/to/cen9623_para.sas7bdat" \
  --travel "NO,UNKNOWN,YES" \
  --cidt "CIDT+,CX+,PARASITIC" \
  --pathogen "CAMPYLOBACTER,SALMONELLA" \
  --states "CA,CO,CT,GA,MD,MN,NM,NY,OR,TN" \
  --chains 2 \
  --iterations 500 \
  --adapt_delta 0.95 \
  --max_treedepth 10 \
  --outdir "results"
```

### Preprocessing Only

To preprocess data without performing analysis:

```bash
nextflow run main.nf \
  -profile singularity \
  -entry PREPROCESS_WORKFLOW \
  --mmwrFile "path/to/mmwr9623_Jan2024.sas7bdat" \
  --outdir "preprocessed" \
  --outputBase "foodnet_data" \
  --generateMetadata true
```

## Running Modes

The pipeline supports three operational modes:

1. **Mode 1: Preprocess data**
   - Clean raw data files and generate metadata
   - Creates standardized CSV files for downstream analysis

2. **Mode 2: Run complete pipeline**
   - Processes raw data files directly
   - Performs full modeling and analysis
   - Generates all outputs in a single workflow

3. **Mode 3: Use existing preprocessed data**
   - Skips preprocessing step
   - Uses previously generated clean CSV files
   - Useful for iterating on analysis parameters without repeating preprocessing

## Input Files

### Required Files

- **MMWR Data File**: SAS format surveillance data (`*.sas7bdat`)
  - Example: `/path/to/mmwr9623_Jan2024.sas7bdat`
- **Census Files**: Population data for incidence rate calculations
  - Bacterial pathogens: `/path/to/cen9623.sas7bdat`
  - Parasitic pathogens: `/path/to/cen9623_para.sas7bdat`

### Expected Data Structure

The MMWR file should contain the following key columns:
- `pathogen`: Type of foodborne pathogen
- `state`: State code
- `year`: Surveillance year
- `sero1`: Serotype information (for Salmonella analysis)
- `travel`: Travel status
- `method`: Diagnostic method

## Pipeline Parameters

### Core Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mmwrFile` | Path to MMWR SAS data file | (Required) |
| `--censusFileB` | Path to census file for bacterial pathogens | (Required) |
| `--censusFileP` | Path to census file for parasitic pathogens | (Required) |
| `--outdir` | Output directory | `output` |
| `--pathogen` | Comma-separated list of pathogens to analyze | `CAMPYLOBACTER,CYCLOSPORA` |

### Filtering Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--travel` | Travel types to include | `NO,UNKNOWN,YES` |
| `--cidt` | Diagnostic methods to include | `CIDT+,CX+,PARASITIC` |
| `--states` | States to include | All available states |
| `--salmonella_serotypes` | Salmonella serotypes to include | All serotypes |

### Model Parameters

| Parameter | Description | Default | 
|-----------|-------------|---------|
| `--chains` | Number of MCMC chains | `2` |
| `--iterations` | Number of iterations per chain | `500` |
| `--adapt_delta` | Adaptation parameter for MCMC | `0.95` |
| `--max_treedepth` | Maximum tree depth for MCMC | `10` |
| `--seed` | Random seed for reproducibility | `123` |

### Preprocessing Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--preprocessed` | Use preprocessed data | `false` |
| `--generateMetadata` | Generate metadata JSON | `true` |
| `--outputBase` | Base name for output files | Derived from input filename |

## Output Structure

```
<outdir>/
├── <projID>/
│   ├── preprocessed/                    # Preprocessed data
│   │   ├── <base>_clean.csv            # Cleaned data
│   │   └── <base>_metadata.json        # Data metadata
│   └── spline_results/                  # Analysis results by pathogen
│       ├── <pathogen>_brm.Rds          # Saved model (R object)
│       ├── <pathogen>_IRCatch.csv      # Incidence rate estimates
│       ├── <pathogen>_summary.txt      # Model summary statistics
│       ├── <pathogen>_site_trends.png  # Site-specific trend plots
│       ├── <pathogen>_overall_trend.png # Overall trend visualization
│       ├── <pathogen>_combined.png     # Combined visualization
│       └── <pathogen>_EstIRRCatch_*.csv # Relative risk comparisons
```

### Output Files

#### Incidence Rate Estimates (`<pathogen>_IRCatch.csv`)
- Contains annual incidence rate estimates by state
- Includes mean and median estimates with uncertainty intervals
- Used for trend analysis and reporting

#### Relative Risk Files (`<pathogen>_EstIRRCatch_*.csv`)
- Compare current incidence rates to historical baseline periods
- Include relative risk ratios and percent changes
- Multiple files represent different comparison periods:
  - 2016-2018: Healthy People 2030 baseline
  - 2020-2022: COVID-19 pandemic period
  - 2004-2006: Early FoodNet baseline
  - 2006-2008: Healthy People 2020 baseline

#### Visualization Files
- Site-specific trends: Individual state trends with uncertainty bands
- Overall trend: Combined national trend with uncertainty
- Combined visualization: Integrated view of site-specific and overall patterns

## Performance Considerations

- **Test Mode**: For validation, use 1 chain and ~100 iterations
- **Full Analysis**: For production, use 2+ chains and 500+ iterations
- **Memory Usage**: Increases with iterations, chains, and pathogen complexity
- **Runtime**: From ~30 minutes (test) to several hours (full analysis)
- **Parallelization**: Multiple pathogens are processed in parallel



## Troubleshooting

| Issue | Solution |
|-------|----------|
| **Model convergence warnings** | Increase `adapt_delta` (>0.95) and `iterations` |
| **Out of memory errors** | Reduce number of chains or run fewer pathogens simultaneously |
| **Missing serotype data** | Check preprocessing with `--generateMetadata true` to verify data |
| **Failed jobs** | Use `-resume` flag to continue from point of failure |
| **Container errors** | Rebuild container with `singularity build --force foodnet.sif foodnet.def` |
| **R package errors** | Check `foodnet.yml` for package compatibility |

### Common Error Messages

#### "Divergent transitions after warmup"
- **Cause**: MCMC sampler is having difficulty exploring the posterior
- **Solution**: Increase adapt_delta to 0.97-0.99 and rerun

#### "Maximum tree depth exceeded"
- **Cause**: Model complexity is causing inefficient sampling
- **Solution**: Increase max_treedepth parameter to 12-15

#### "Error in read_sas()"
- **Cause**: SAS file format issues or file not found
- **Solution**: Verify file paths and ensure SAS files are not corrupted

## Citations

This pipeline uses statistical methods from:

- R. McElreath. *Statistical Rethinking: A Bayesian Course with Examples in R and Stan*. Chapman and Hall/CRC, 2020.
- P. Bürkner. *brms: An R Package for Bayesian Multilevel Models Using Stan*. Journal of Statistical Software, 80(1), 1-28, 2017.

## License and Contributors

This software is released under the CDC Public Domain License.

**Contributors**:
- Joshua Forstedt (CDC/NCEZID/DIDRI/OAMD)
- Samantha Sevilla (CDC/NCEZID/DIDRI/OAMD)
- Daniel Weller (CDC/DFWED/EDEB)
- Beau Bruce (CDC/DFWED/EDEB)
- Erica Billig Rose (CDC/DFWED/EDEB)
