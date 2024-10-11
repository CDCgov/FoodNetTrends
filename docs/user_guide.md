# Getting Started
To use the FoodNet enhanced model pipeline, you will need:
- R (version 4.3.2) with `brms` and `tidybayes` packages.
- Nextflow (version 24.04.2) for running the pipeline.
- Access to the FoodNet surveillance dataset.

For a detailed list of software and dependencies, please refer to the [GitHub repository](https://github.com/CDCgov/FoodNetTrends).

# Preparing Files
Data preparation is crucial for accurate modeling. The FoodNet enhanced model requires the following inputs:
- Site-level data aggregated annually.
- Data formatted for Bayesian analysis using the `brms` package.
- Metadata regarding county-year levels for accurate site identification.

Ensure that data files follow the specified format to prevent errors during pipeline execution.

# Running the Pipeline
The enhanced FoodNet model is implemented using Nextflow. To run the pipeline:
1. Clone the repository:
   ```bash
   git clone https://github.com/CDCgov/FoodNetTrends
   cd FoodNetTrends
   
3. Execute the pipeline:
   ```bash
   nextflow run main.nf
   
# Interpreting Output
The output includes:
- Posterior predictive distributions for site-specific trends.
- Median incidence estimates with 50%, 75%, and 95% credibility intervals.
- Summary statistics comparing the original and enhanced models using R² and ELPD.

Refer to the visualizations in the output folder for graphical comparisons between models.

# Test Data
A set of test data is included to validate your setup:
- Download the test data from [here](https://github.com/CDCgov/FoodNetTrends/test_data).
