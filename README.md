# FoodNet Trends Pipeline

## Introduction
**FoodNet Trends** is a bioinformatics pipeline that performs spline-based modeling of foodborne illness surveillance data. The pipeline processes FoodNet MMWR data and applies Bayesian hierarchical models to estimate incidence rates and trends across different pathogens and sites.

## Features
1. Preprocesses raw MMWR surveillance data
2. Applies Bayesian hierarchical models with splines for flexible trend analysis
3. Generates incidence rate estimates with uncertainty intervals
4. Creates visualizations of pathogen-specific trends
5. Calculates comparative statistics across different time periods
6. Organizes results in a structured output directory
7. Supports analysis of multiple pathogens in a single run

## Requirements
- Nextflow (version 23.04.0 or later)
- Singularity (version 4.0.0 or later)
- Java (version 17 or later)
- SGE cluster environment (for distributed computing)

## Input Data
The pipeline requires the following input data:

- **mmwrFile**: Path to FoodNet MMWR SAS data file
  - Example: `/path/to/mmwr9623_Jan2024.sas7bdat`
- **censusFileB**: Path to census data for bacterial pathogens
  - Example: `/path/to/cen9623.sas7bdat`
- **censusFileP**: Path to census data for parasitic pathogens
  - Example: `/path/to/cen9623_para.sas7bdat`

## Parameters
### Input/Output Parameters
- **outdir**: Output directory (default: `output`)
- **projID**: Project identifier (default: auto-generated timestamp)
- **pathogen**: Comma-separated list of pathogens to analyze (default: `CAMPYLOBACTER,CYCLOSPORA`)
  - Supported pathogens: `CAMPYLOBACTER`, `CYCLOSPORA`, `SALMONELLA`, `SHIGELLA`, `STEC`, `VIBRIO`, `YERSINIA`

### Data Filtering Parameters
- **travel**: Travel types to include (default: `NO,UNKNOWN,YES`)
- **cidt**: CIDT types to include (default: `CIDT+,CX+,PARASITIC`)
- **preprocessed**: Whether to use preprocessed data (default: `false`)
- **cleanFile**: Path to cleaned CSV file when using preprocessed data

### Model Parameters
- **chains**: Number of MCMC chains (default: `2`)
- **iterations**: Number of MCMC iterations per chain (default: `500`)
- **adapt_delta**: Adaptation parameter for HMC (default: `0.95`)
- **max_treedepth**: Maximum tree depth for HMC (default: `10`)
- **seed**: Random seed for reproducibility (default: `123`)

## Running the Pipeline

### 1. Interactive Mode (Recommended)
The easiest way to run the pipeline is using the interactive script:

```bash
./run_workflow.sh
```

The interactive script will guide you through:
- Selecting run mode (Test, Full analysis, or Resume)
- Choosing pathogens to analyze (individual or all available)
- Setting MCMC parameters (chains, iterations, etc.)
- Running in background or foreground

### 2. Direct Nextflow Command
For advanced users or automated workflows, you can run the pipeline directly:

```bash
module load nextflow/24.10.4 singularity/4.1.4 java/17.0.6

nextflow run main.nf \
  -profile singularity \
  -entry SPLINE \
  --mmwrFile "/path/to/mmwr9623_Jan2024.sas7bdat" \
  --censusFileB "/path/to/cen9623.sas7bdat" \
  --censusFileP "/path/to/cen9623_para.sas7bdat" \
  --pathogen "CAMPYLOBACTER,SALMONELLA" \
  --chains 2 \
  --iterations 500 \
  --adapt_delta 0.95 \
  --max_treedepth 10 \
  --outdir "output"
```

### 3. Run Profiles
The pipeline includes preconfigured profiles:

- **Test Profile**: For quick validation with minimal resources
  ```bash
  ./run_workflow.sh test
  ```

- **Production Profile**: For full analysis with robust settings
  ```bash
  ./run_workflow.sh full
  ```

- **Resume Execution**: To continue an interrupted run
  ```bash
  ./run_workflow.sh resume
  ```

## Output
The pipeline generates a structured output directory:

```
output/
└── [projID]/
    ├── pipeline_info/            # Execution reports and logs
    │   ├── execution_report.html # Pipeline execution report
    │   ├── execution_trace.txt   # Detailed execution trace
    │   └── pipeline_dag.html     # Execution graph
    ├── preprocessed/             # Preprocessed data files
    │   └── clean_mmwr.csv        # Cleaned MMWR data
    └── spline_results/           # Model results for each pathogen
        ├── [pathogen]_brm.Rds           # Saved model object
        ├── [pathogen]_IRCatch.csv       # Incidence rate estimates
        ├── [pathogen]_summary.txt       # Model summary statistics
        ├── [pathogen]_site_trends.png   # Site-specific trend plots
        ├── [pathogen]_overall_trend.png # Overall trend plot
        ├── [pathogen]_combined.png      # Combined visualization
        └── [pathogen]_EstIRRCatch_*.csv # Relative risk comparisons
```

## Visualizations
The pipeline generates several visualizations:

1. **Site-specific trends**: Incidence trends for each surveillance site
2. **Overall trend**: Combined trend across all sites
3. **Combined visualization**: Integrated view of site-specific and overall trends

## Relative Risk Analysis
The pipeline calculates relative risks and percent changes compared to baseline periods:
- 2016-2018 (federal goals baseline)
- 2020-2022 (most recent 3 years)
- 2004-2006 (earliest years)
- 2006-2008 (historical baseline)
- 2010-2012 (historical baseline)

## Performance Considerations
- For test runs, use `test` mode with reduced parameters (chains=1, iterations=100)
- For production runs, consider using at least 2 chains with 500+ iterations
- Running with multiple pathogens will launch parallel jobs on the cluster
- Background execution is recommended for long-running analyses

## Troubleshooting
Common issues and solutions:

1. **Convergence Warnings**: Increase adapt_delta (0.95+) and iterations
2. **Memory Errors**: Reduce the number of chains or run with fewer pathogens
3. **Failed Jobs**: Use the resume feature to continue from the point of failure

## Credits
The FoodNet Trends pipeline was developed by Samantha Sevilla, Josh Forstedt, and OAMD's SciComp Team with support from Daniel Weller (CDC/DFWED/EDEB), based on R scripts developed by Daniel Weller (CDC/DFWED/EDEB) with support from Beau Bruce (CDC/DFWED/EDEB) and Erica Billig Rose (CDC/DFWED/EDEB).

## Contributions and Support
If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

## License
This software is released under the CDC Public Domain License.