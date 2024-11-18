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
