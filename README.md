# FoodNet Trends Pipeline

## Introduction
**FoodNet Trends** is a bioinformatics pipeline that performs spline-based modeling of foodborne illness surveillance data. The pipeline processes FoodNet MMWR data and applies Bayesian hierarchical models to estimate incidence rates and trends across different pathogens and sites.

## Features
1. Preprocesses raw MMWR surveillance data
2. Applies Bayesian hierarchical models with splines for flexible trend analysis
3. Generates incidence rate estimates with uncertainty intervals
4. Creates visualizations of pathogen-specific trends
5. Organizes results in a structured output directory

## Usage
The pipeline requires the following input data:

- **mmwrFile**: Path to FoodNet MMWR SAS data file
  - Example: `/path/to/mmwr9623_Jan2024.sas7bdat`
- **censusFile_B**: Path to census data for bacterial pathogens
  - Example: `/path/to/cen9623.sas7bdat`
- **censusFile_P**: Path to census data for parasitic pathogens
  - Example: `/path/to/cen9623_para.sas7bdat`

Optional parameters with defaults:
- **outdir**: Output directory (default: `output`)
- **travel**: Travel types to include (default: `"NO,UNKNOWN,YES"`)
- **cidt**: CIDT types to include (default: `"CIDT+,CX+,PARASITIC"`)
- **projID**: Project identifier (default: auto-generated timestamp)

### Running the Pipeline

There are two ways to run the pipeline:

#### 1. Direct Nextflow Command
module load nextflow/24.10.4 singularity/4.1.4 java/17.0.6

nextflow run main.nf
-profile singularity
-entry FoodNetTrends
--mmwrFile "/path/to/mmwr9623_Jan2024.sas7bdat"
--censusFile_B "/path/to/cen9623.sas7bdat"
--censusFile_P "/path/to/cen9623_para.sas7bdat"
--outdir "output"


#### 2. Using the `run_workflow.sh` Script
Run test mode (default)
./run_workflow.sh

Run in background
./run_workflow.sh test --bg

Resume a previous run
./run_workflow.sh resume

Specify custom output directory
./run_workflow.sh test --outdir=/path/to/output


## Output
The pipeline generates a structured output directory:
output/
└── [projID]/
├── pipeline_info/ # Execution reports and logs
├── preprocess_results/ # Preprocessed data files
└── spline_results/ # Model results and visualizations
├── [pathogen]_brm.Rds # Saved model object
├── [pathogen]_IRCatch.csv # Incidence rate estimates
├── [pathogen]_site_trends.png # Site-specific trend plots
├── [pathogen]_overall_trend.png # Overall trend plot
└── [pathogen]_combined.png # Combined visualization


## Credits
The FoodNet Trends pipeline was developed by Samantha Sevilla, Josh Forstedt,  and OAMD's SciComp Team with support from Daniel Weller (CDC/DFWED/EDEB), based on R scripts developed by Daniel Weller (CDC/DFWED/EDEB) with support from Beau Bruce (CDC/DFWED/EDEB) and Erica Billig Rose (CDC/DFWED/EDEB).

## Contributions and Support
If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

