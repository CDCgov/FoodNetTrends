# FoodNet Trends Pipeline

## Introduction
**FoodNet Trends** is a bioinformatics pipeline that performs spline-based modeling of foodborne illness surveillance data. The pipeline processes FoodNet MMWR data and applies Bayesian hierarchical models to estimate incidence rates and trends across different pathogens and sites.

## Features
1. Preprocesses raw MMWR surveillance data with pathogen name standardization
2. Reuses existing preprocessed data to save computation time
3. Flexible pathogen grouping options for STEC and Salmonella analyses
4. Applies Bayesian hierarchical models with splines for flexible trend analysis
5. Generates incidence rate estimates with uncertainty intervals
6. Creates visualizations of pathogen-specific trends
7. Calculates comparative statistics across different time periods
8. Configurable serotype recoding and catchment area definitions
9. Supports analysis of multiple pathogens in a single run

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
  - Additional pathogens automatically detected from data in pathogen-agnostic mode

### Data Filtering Parameters
- **travel**: Travel types to include (default: `NO,UNKNOWN,YES`)
- **cidt**: CIDT types to include (default: `CIDT+,CX+,PARASITIC`)
- **states**: Comma-separated list of states to include (default: `null` = all states)
- **preprocessed**: Whether to use preprocessed data (default: `false`)
- **cleanFile**: Path to cleaned CSV file when using preprocessed data
- **skip_preprocessing**: Skip the preprocessing step entirely (default: `false`)

### Data Processing Parameters
- **matching_sensitivity**: Pathogen name standardization sensitivity (default: `MEDIUM`)
  - `STRICT`: Exact matches only
  - `MEDIUM`: Exact + prefix matching + limited fuzzy matching
  - `RELAXED`: All matching methods + broader fuzzy matching
- **pathogen_grouping**: Configuration for pathogen subgroup analysis
  - Automatically configured through interactive prompts

### Model Parameters
- **chains**: Number of MCMC chains (default: `2`)
- **iterations**: Number of MCMC iterations per chain (default: `500`)
- **adapt_delta**: Adaptation parameter for HMC (default: `0.95`)
- **max_treedepth**: Maximum tree depth for HMC (default: `10`)
- **seed**: Random seed for reproducibility (default: `123`)

### Configuration Parameters (Optional)
- **serotype_config**: Path to CSV file with custom serotype recoding rules
- **catchment_config**: Path to CSV file with custom catchment area definitions
  - See example configurations in `analysis_configs/examples/`

## Running the Pipeline

### 1. Interactive Mode (Recommended)
The easiest way to run the pipeline is using the interactive script:

```bash
./run_workflow.sh
```

The interactive script will guide you through:
- Selecting run mode (Test, Full analysis, or Resume)
- Checking for and optionally reusing existing preprocessed data
- Choosing pathogens to analyze (individual or all available)
- Configuring pathogen subgroup analysis for STEC and Salmonella
- Setting serotype and catchment configurations
- Selecting pathogen name matching sensitivity
- Setting MCMC parameters (chains, iterations, etc.)
- Running in background or foreground

#### Preprocessed Data Reuse
The pipeline automatically searches for existing preprocessed data from previous runs. When found, you can:
- View file metadata (creation date, size)
- Check if preprocessing reports are available
- Choose to reuse the data or run preprocessing again
- Automatically extract available pathogens from the preprocessed file

#### Pathogen Grouping Options
For STEC analysis:
- Combined: Analyze all STEC together
- O157 only: Focus on STEC O157
- Non-O157 only: Focus on non-O157 STEC
- Both separately: Run separate analyses for O157 and non-O157

For Salmonella analysis:
- Combined: Analyze all serotypes together
- Custom selection: View ranked list of serotypes by frequency and select specific ones

Serotypes are automatically grouped with non-informative values (NOT SPECIATED, UNKNOWN, etc.) classified as "Missing"

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

### 3. Preprocessing Only
To run only the preprocessing step without the full analysis:

```bash
nextflow run main.nf \
  -profile singularity \
  -entry PREPROCESS_ONLY \
  --mmwrFile "/path/to/mmwr9623_Jan2024.sas7bdat" \
  --censusFileB "/path/to/cen9623.sas7bdat" \
  --censusFileP "/path/to/cen9623_para.sas7bdat" \
  --outdir "output"
```

This generates the preprocessed data files without running the Bayesian modeling, useful for data validation or when you want to inspect the cleaned data before full analysis.

### 4. Run Profiles
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

### Using Configuration Files
To use custom serotype or catchment configurations:

```bash
nextflow run main.nf \
  --mmwrFile /path/to/mmwr.sas7bdat \
  --censusFileB /path/to/census_b.sas7bdat \
  --censusFileP /path/to/census_p.sas7bdat \
  --serotype_config analysis_configs/custom_serotypes.csv \
  --catchment_config analysis_configs/custom_catchment.csv
```

Example configuration files are provided in `analysis_configs/examples/`.

### Configuration File Formats

#### Serotype Configuration (CSV)
Controls how serotypes are recoded. Format:
```csv
serotype_value,replacement,pathogen,match_type
NOT SPECIATED,Missing,all,exact
Typhimurium var. 5-,Typhimurium,SALMONELLA,contains
```

#### Catchment Configuration (CSV)
Defines which states participate in surveillance by year. Format:
```csv
state,start_year,end_year,pathogen_type
CA,1996,9999,both
CO,2001,9999,both
```

## Output
The pipeline generates a structured output directory:

```
output/
└── [projID]/
    ├── pipeline_info/                           # Execution reports and logs
    │   ├── execution_report_*.html              # Pipeline execution report
    │   ├── execution_trace_*.txt                # Detailed execution trace
    │   ├── execution_timeline_*.html            # Execution timeline visualization
    │   └── pipeline_dag_*.html                  # Execution graph
    ├── preprocessed/                            # Preprocessed data files
    │   ├── clean_mmwr.csv                       # Cleaned MMWR data
    │   ├── clean_mmwr_preprocessing_report.csv  # Pathogen name standardization report
    │   ├── resource_profile.csv                 # Data complexity metrics
    │   ├── metadata_states.csv                  # State-level summary
    │   ├── metadata_cidt.csv                    # Diagnostic method summary
    │   └── metadata_travel.csv                  # Travel status summary
    └── spline_results/                          # Model results for each pathogen
        ├── [pathogen]_brm.Rds                   # Saved model object
        ├── [pathogen]_IRCatch.csv               # Incidence rate estimates
        ├── [pathogen]_summary.txt               # Model summary statistics
        ├── [pathogen]_site_trends.png           # Site-specific trend plots
        ├── [pathogen]_overall_trend.png         # Overall trend plot
        ├── [pathogen]_combined.png              # Combined visualization
        └── [pathogen]_EstIRRCatch_*.csv         # Relative risk comparisons
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
- For production runs, use `-profile production` which sets optimized parameters:
  - chains=4, iterations=2000, adapt_delta=0.99, max_treedepth=15
- Running with multiple pathogens will launch parallel jobs on the cluster
- Background execution is recommended for long-running analyses

## Troubleshooting
Common issues and solutions:

1. **Convergence Warnings**: Increase adapt_delta (0.95+) and iterations
2. **Memory Errors**: Reduce the number of chains or run with fewer pathogens
3. **Failed Jobs**: Use the resume feature to continue from the point of failure
4. **No Data Found for Pathogen**: Check pathogen spelling matches preprocessed data exactly
5. **Pathogen Name Variations**: Adjust matching_sensitivity from STRICT to MEDIUM or RELAXED
6. **Missing Serotype Data**: Verify serotypesummary column exists in data
7. **STEC Grouping Issues**: Ensure stec_class column is present for O157/non-O157 splitting

## Credits
The FoodNet Trends pipeline was developed by Samantha Sevilla, Josh Forstedt, and OAMD's SciComp Team with support from Daniel Weller (CDC/DFWED/EDEB), based on R scripts developed by Daniel Weller (CDC/DFWED/EDEB) with support from Beau Bruce (CDC/DFWED/EDEB) and Erica Billig Rose (CDC/DFWED/EDEB).

## Contributions and Support
If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

## License
This software is released under the CDC Public Domain License.

## SHARE IT Act Compliance
 
```
Organization: NCEZID/AMD
 
Contact email: ncezid_shareit@cdc.gov
 
Exemption: NA
 
Exemption Justification: NA
 
Description fields:
```
