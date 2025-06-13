process RESOURCE_PROFILER {
    tag "Profiling resources for pathogen analysis"
    label 'process_low'
    shell "/bin/bash"
    container 'foodnet.sif'

    publishDir "${params.outdir}/${params.projID}/preprocessed", mode: 'copy'

    input:
    path cleanFile

    output:
    path "resource_profile.csv", emit: profile

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
    })

    # Read the cleaned data
    data <- read_csv("${cleanFile}", show_col_types = FALSE)
    
    # Ensure column names are lowercase (matching preprocess.R output)
    names(data) <- tolower(names(data))
    
    # Debug: Print column names and first few rows
    cat("\\nColumn names:", paste(names(data), collapse=", "), "\\n")
    cat("Number of rows:", nrow(data), "\\n")
    if (nrow(data) > 0) {
        cat("First few pathogens:", paste(head(unique(data\$pathogen), 10), collapse=", "), "\\n")
    }
    
    # Calculate metrics for each pathogen
    pathogen_metrics <- data %>%
        filter(!is.na(pathogen)) %>%  # Filter out NA values only
        group_by(pathogen) %>%
        summarise(
            rows = n(),
            sites = n_distinct(state),
            years = n_distinct(year),
            .groups = 'drop'
        ) %>%
        mutate(
            # Complexity score: rows * sites * years
            complexity = rows * sites * years,
            # Data size category for logging
            size_category = case_when(
                rows > 50000 ~ "extra_large",
                rows > 20000 ~ "large",
                rows > 10000 ~ "medium",
                rows > 5000 ~ "small",
                TRUE ~ "tiny"
            )
        )
    
    # Check if we have any valid pathogens
    if (nrow(pathogen_metrics) == 0) {
        stop("No valid pathogens found in the data. Check if the pathogen column contains proper pathogen names.")
    }
    
    # Write to CSV
    write_csv(pathogen_metrics, "resource_profile.csv")
    
    # Print summary for logging
    cat("\\nResource Profile Summary:\\n")
    cat("========================\\n")
    for (i in 1:nrow(pathogen_metrics)) {
        cat(sprintf("%-15s: %6d rows, %2d sites, %2d years (complexity: %d, category: %s)\\n", 
                    pathogen_metrics\$pathogen[i], 
                    pathogen_metrics\$rows[i], 
                    pathogen_metrics\$sites[i], 
                    pathogen_metrics\$years[i], 
                    pathogen_metrics\$complexity[i], 
                    pathogen_metrics\$size_category[i]))
    }
    """
}