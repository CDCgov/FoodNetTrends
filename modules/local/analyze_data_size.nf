process ANALYZE_DATA_SIZE {
    tag "Analyzing data size for resource allocation"
    label 'process_low'
    shell "/bin/bash"
    container 'foodnet.sif'

    input:
    path cleanFile

    output:
    path "pathogen_metrics.json", emit: metrics

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(dplyr)
        library(jsonlite)
        library(readr)
    })

    # Read the cleaned data
    data <- read_csv("${cleanFile}", show_col_types = FALSE)
    
    # Ensure column names are lowercase (matching preprocess.R output)
    names(data) <- tolower(names(data))
    
    # Calculate metrics for each pathogen
    pathogen_metrics <- data %>%
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
    
    # Convert to named list for JSON output
    metrics_list <- list()
    for (i in 1:nrow(pathogen_metrics)) {
        p <- pathogen_metrics\$pathogen[i]
        metrics_list[[p]] <- list(
            rows = pathogen_metrics\$rows[i],
            sites = pathogen_metrics\$sites[i],
            years = pathogen_metrics\$years[i],
            complexity = pathogen_metrics\$complexity[i],
            size_category = pathogen_metrics\$size_category[i]
        )
    }
    
    # Write to JSON
    write_json(metrics_list, "pathogen_metrics.json", pretty = TRUE)
    
    # Also print summary for logging
    cat("\\nPathogen Data Size Analysis:\\n")
    cat("============================\\n")
    for (p in names(metrics_list)) {
        m <- metrics_list[[p]]
        cat(sprintf("%-15s: %6d rows, %2d sites, %2d years (complexity: %d, category: %s)\\n", 
                    p, m\$rows, m\$sites, m\$years, m\$complexity, m\$size_category))
    }
    """
}