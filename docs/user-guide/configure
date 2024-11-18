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
For additional documentation and support, please refer to the [GitHub repository](https://github.com/CDCgov/FoodNetTrends).
