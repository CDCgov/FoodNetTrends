# FoodNet Trends Analysis Pipeline v1.0

A Nextflow-based pipeline for Bayesian modeling of foodborne disease surveillance data.

## Overview

The FoodNet Trends pipeline processes food-borne illness surveillance data from the CDC's Foodborne Diseases Active Surveillance Network (FoodNet), applies Bayesian hierarchical models with splines, and generates standardized incidence rates, trend analysis, and visualizations. The pipeline is designed for epidemiologists and statisticians analyzing surveillance data for public health decision-making.

## Features

- **Data Preprocessing**: Cleans and standardizes raw MMWR SAS data files
- **Bayesian Modeling**: Implements hierarchical models with splines for flexible trend analysis
- **Multi-pathogen Support**: Analyzes all pathogen types present in the supplied surveillance data
- **Travel-adjustment**: Includes option for travel-adjusted incidence rates
- **CIDT-adjustment**: Controls for culture-independent diagnostic test effects
- **Flexible Input**: Works with both raw SAS files and preprocessed CSV data
- **Interactive Dashboard**: Self-contained HTML visualization of results requiring no server setup
- **Reproducible Results**: Consistent model settings and complete logging
- **Low Code Required**: Interactive shell script guides analysis configuration

## Version Information

- **Current Version**: 1.0 (May 2025)
- **Release Status**: Production/Stable
- **Dependencies**:
  - Nextflow ≥ 24.10.4
  - Singularity ≥ 4.1.4
  - R ≥ 4.4.0 (via container)
  - brms R package ≥ 3.0.0 (via container)
  - R visualization packages: plotly, DT, htmlwidgets (for dashboard)

## Pipeline Structure

The FoodNet Trends pipeline consists of the following primary components:

```
FoodNetTrends/
├── assets/                   # Assets for visualization
│   └── dashboard_template.html # Dashboard HTML template
├── bin/                      # Core R scripts for analysis
│   ├── trendy.R              # Main Bayesian modeling script
│   ├── functions.R           # Shared statistical functions
│   ├── calcIR.R              # Data preprocessing script
│   └── generate_dashboard.R  # Interactive dashboard generator
├── conf/                     # Configuration profiles
├── modules/                  # Nextflow processes
│   └── local/
│       ├── trendy.nf         # Bayesian modeling process
│       └── preprocess.nf     # Data preprocessing process
├── workflows/                # Workflow definitions
│   ├── spline.nf             # Main spline analysis workflow
│   └── preprocess.nf         # Data preprocessing workflow
└── run_workflow.sh           # Interactive execution script
```

## Output Files

The pipeline generates output files in a standardized folder structure:

```
results/
└── <project_id>/
    ├── preprocessed/         # Preprocessed data files
    │   ├── <base>.csv            # Cleaned data
    │   ├── metadata/             # Dataset metadata directory
    │   |   └── <base>_metadata.json  # Dataset metadata 
    │   └── logs/                 # Preprocessing logs
    ├── spline_results/       # Modeling results for each pathogen
    │   ├── <pathogen>_brm.Rds          # Saved model (R object)
    │   ├── <pathogen>_IRCatch.csv      # Incidence rate estimates
    │   ├── <pathogen>_summary.txt      # Model summary statistics
    │   ├── <pathogen>_site_trends.png  # Site-specific trend plots
    │   ├── <pathogen>_overall_trend.png # Overall trend visualization
    │   ├── <pathogen>_combined.png     # Combined visualization
    │   ├── <pathogen>_EstIRRCatch_<comparisonYears>.csv # Relative risk comparisons
    │   └── logs/                       # Per-pathogen log directory
    ├── <project_id>_YYYYMMDD_dashboard.html # Interactive dashboard
    └── pipeline_info/        # Pipeline execution information
```

### Standardized File Naming Convention

The pipeline uses a consistent file naming convention for all outputs to improve organization and discovery:

```
<identifier>_<filetype>[_<subtype>].<extension>
```

Where:
- `<identifier>` is either the pathogen name (e.g., "CAMPYLOBACTER") or the dataset base name
- `<filetype>` indicates the content type (e.g., "brm", "IRCatch", "site_trends")
- `<subtype>` (optional) provides additional context for specialized files (e.g., comparison years)
- `<extension>` is the standard file extension (e.g., "csv", "Rds", "png", "txt")

#### Incidence Rate Estimates (`<pathogen>_IRCatch.csv`)

The main result file contains incidence rate estimates by year and state, with the following columns:

- `state`: State or site abbreviation
- `year`: Year of estimate
- `median_incidence`: Median incidence rate per 100,000 population
- `lower_hdi`, `upper_hdi`: Lower and upper 95% credible interval bounds
- `pathogen`: Pathogen name
- `travel`: Travel adjustment status
- `culture`: Culture confirmation status

#### Relative Risk Files (`<pathogen>_EstIRRCatch_<years>.csv`)

Comparison files show relative risks and percent changes between the most recent year and historical periods:

- `state`: State or site abbreviation 
- `year`: Current year
- `comparison_period`: Historical years used for comparison (e.g., "2016-2018")
- `current_incidence`: Current incidence rate
- `period_incidence`: Average incidence for comparison period
- `relative_risk`: Ratio of current to historical incidence 
- `percent_change`: Percent change from historical period

#### Interactive Dashboard (`<project_id>_YYYYMMDD_dashboard.html`)

The dashboard is a self-contained HTML file that provides an interactive visualization of all analysis results:

- **No server required**: Open directly in any modern web browser
- **Embedded data**: All data is contained within the HTML file
- **Interactive visualization**: Filter, zoom, and explore results
- **Multiple views**: Trend analysis, geographic distribution, and relative risk comparisons
- **Data tables**: Sortable and searchable data tables for detailed exploration

## Interactive Dashboard

The FoodNet Trends pipeline automatically generates an interactive HTML dashboard that allows users to explore and visualize results without requiring any server setup or additional software.

### Dashboard Features

- **Pathogen filtering**: View results for specific pathogens or compare multiple pathogens
- **Time period selection**: Focus on specific years or examine long-term trends
- **Geographic visualization**: View state-by-state distribution of incidence rates
- **Trend analysis**: Interactive time series plots with confidence intervals
- **Relative risk comparison**: Compare current rates with historical baselines
- **Responsive design**: Works on desktop and tablet devices
- **Data tables**: Sortable and filterable data tables for detailed analysis
- **Export capability**: Download visualizations and data for presentations or reports

### Accessing the Dashboard

After pipeline completion, the dashboard can be found in the project output directory with the naming pattern `<project_id>_YYYYMMDD_dashboard.html`. To view:

1. Navigate to the output directory
2. Open the dashboard HTML file with any modern web browser:
   ```bash
   firefox results/<project_id>/<project_id>_YYYYMMDD_dashboard.html
   ```

### Dashboard Customization

The dashboard appearance and behavior can be customized through pipeline parameters:

- `--enable_dashboard`: Enable/disable dashboard generation (default: true)
- `--dashboard_title`: Custom title for the dashboard
- `--dashboard_logo`: Path to custom logo image for branding

## Requirements

- Nextflow ≥ 24.10.4
- Singularity ≥ 4.1.4
- Java ≥ 17
- Access to HPC environment

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
  --mmwrFile "path/to/mmwr9624_May2025.sas7bdat" \
  --censusFileB "path/to/cen9624.sas7bdat" \
  --censusFileP "path/to/cen9624_para.sas7bdat" \
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
  --mmwrFile "path/to/mmwr9624_May2025.sas7bdat" \
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
  - Example: `/path/to/mmwr9624_May2025.sas7bdat`
- **Census Files**: Population data for incidence rate calculations
  - Bacterial pathogens: `/path/to/cen9624.sas7bdat`
  - Parasitic pathogens: `/path/to/cen9624_para.sas7bdat`

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

## Performance Considerations

- **Test Mode**: For validation, use 1 chain and ~100 iterations
- **Full Analysis**: For production, use 2+ chains and 500+ iterations
- **Memory Usage**: Increases with iterations, chains, and pathogen complexity
- **Runtime**: From ~30 minutes (test) to several hours (full analysis)
- **Parallelization**: Multiple pathogens are processed in parallel

## Monitoring Progress

The pipeline includes a robust, environment-agnostic progress tracking system that works reliably in any environment (local, HPC, containers, etc.):

### Using the Progress Monitor Script

A dedicated monitoring script provides real-time progress visualization:

```bash
# Basic usage - monitor current directory
./bin/monitor_progress.sh

# Monitor specific directory
./bin/monitor_progress.sh -d /path/to/results

# Monitor continuously, updating every 5 seconds
./bin/monitor_progress.sh --follow

# Monitor only a specific pathogen
./bin/monitor_progress.sh --pathogen SALMONELLA

# All options
./bin/monitor_progress.sh --dir /path/to/results --pathogen CAMPYLOBACTER --follow --interval 10
```

**Script Options:**

| Option | Description |
|--------|-------------|
| `-d, --dir DIR` | Directory to monitor for progress files (default: current dir) |
| `-p, --pathogen NAME` | Only show progress for specific pathogen |
| `-f, --follow` | Continuously update (like 'tail -f') |
| `-i, --interval SEC` | Refresh interval in seconds (default: 5) |
| `-h, --help` | Show help message |

The monitor script provides a detailed, color-coded view of:
- Overall progress percentage for each pathogen with visual progress bar
- Current analysis stage and milestone
- Elapsed time and estimated time remaining
- Recent activity logs for context
- Files generated so far

### Progress File System

The progress tracking system generates three types of files during pipeline execution:

1. **`[PATHOGEN]_progress.txt`**: Structured progress information including:
   - Current percentage completion
   - Processing stage
   - Status message
   - Time elapsed and remaining
   - Current milestone

2. **`[PATHOGEN]_progress_log.txt`**: Comprehensive log with timestamps for all stages and milestones, useful for debugging or understanding the analysis process.

3. **`[PATHOGEN]_percent.txt`**: Simple file containing only the percentage completion number, designed for easy polling by other monitoring systems.

### Integration with Custom Monitoring Systems

The progress files use a simple format that can be easily integrated with custom monitoring systems, dashboards, or notification scripts:

```bash
# Example: Simple script to send notification when analysis reaches 75%
watch -n 60 'for f in *_percent.txt; do
  pct=$(cat $f);
  if [ "$pct" -ge 75 ] && [ "$pct" -lt 76 ]; then
    pathogen=${f%_percent.txt};
    echo "$pathogen analysis is 75% complete" | mail -s "Analysis Progress" user@example.com;
  fi;
done'
```

### Monitoring in HPC Environments

In HPC environments where direct terminal access might be limited, the progress tracking system provides several options:

1. **Scheduled Status Checks**: Set up a periodic job to run `monitor_progress.sh` and output to a file
2. **File-Based Monitoring**: Check the progress files directly from a shared filesystem
3. **Email Notifications**: Create a simple script to send email alerts at specific milestones

The file-based progress tracking ensures visibility into analysis progress regardless of the execution environment.

### Troubleshooting Progress Monitoring

If you're having trouble seeing progress information:

1. **Check Work Directories**: Progress files are created in the Nextflow work directories. Use:
   ```bash
   find $(pwd) -name "*_progress.txt" -o -name "*_percent.txt"
   ```
   
2. **Monitor Multiple Directories**: On HPC systems, specify multiple directories to monitor:
   ```bash
   ./bin/monitor_progress.sh -d /scicomp/scratch/username/nextflow/work
   ```
   
3. **Increase Search Depth**: Use the recursive option in find command:
   ```bash
   find /scicomp/scratch/username/nextflow/work -maxdepth 5 -name "*_progress.txt"
   ```

4. **Check Permissions**: Ensure that progress files are being created with appropriate permissions:
   ```bash
   ls -la $(find . -name "*_progress.txt" | head -1)
   ```

5. **Terminal Compatibility**: If progress bars appear garbled, use the `--no-color` option to disable colored output.
   ```bash
   ./bin/monitor_progress.sh --dir /path/to/work --no-color
   ```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| **Model convergence warnings** | Increase `adapt_delta` (>0.95) and `iterations` |
| **Out of memory errors** | Reduce number of chains or run fewer pathogens simultaneously |
| **Missing serotype data** | Check preprocessing with `--generateMetadata true` to verify data |
| **Failed jobs** | Use `-resume` flag to continue from point of failure |
| **Container errors** | Rebuild container with `singularity build --force foodnet.sif foodnet.def` |
| **R package errors** | Check `foodnet.yml` for package compatibility |
| **No progress visible** | Use `./bin/monitor_progress.sh` to check progress files |

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
