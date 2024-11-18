## Getting Started

To use the FoodNet enhanced model pipeline, you will need:

- **R (version 4.3.2)** with the `brms` and `tidybayes` packages.
- **Nextflow (version 24.04.2)** for running the pipeline.
- Access to the **FoodNet surveillance dataset**.
- **Singularity** or **Docker** for containerization (Singularity is used in the examples).

For a detailed list of software and dependencies, please refer to the [GitHub repository](https://github.com/CDCgov/FoodNetTrends).

We highly recommend using Docker or Singularity containers for full pipeline reproducibility. If these are not possible, Conda is also supported. See the [`-profile`](#-profile) section for more information.
