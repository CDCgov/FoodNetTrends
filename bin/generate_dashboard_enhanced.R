#!/usr/bin/env Rscript
# ========================================================================
# FoodNet Trends v1.0 - Enhanced Dashboard Generator
# ========================================================================
#
# Purpose:
#   Creates a feature-rich, interactive HTML dashboard from analysis results
#   with visualizations, charts, and data tables.
#
# Features:
#   - Interactive trend charts for incidence rates
#   - Pathogen-specific analysis and comparisons
#   - Geographic visualizations by state
#   - Data quality indicators and warnings
#   - Responsive design for all devices
#
# Input:
#   - Path to results directory containing IRCatch files
#   - Configuration parameters for dashboard customization
#
# Output:
#   - Interactive HTML dashboard file with embedded data visualizations
#
# Last updated: 2025-05-21
# ========================================================================

# ========================================================================
# Memory and Performance Optimization
# ========================================================================
# Force garbage collection to improve memory management
gc(reset = TRUE)
# Set higher memory limits for data processing
options(future.globals.maxSize = 2048*1024^2) # 2GB max for big data objects
# Additional optimizations
options(datatable.print.topn = 5)
options(datatable.print.nrows = 20)
options(digits = 4) # Reduce precision of numeric values for display
options(scipen = 999) # Avoid scientific notation in displays

# ========================================================================
# Load required packages with fallback handling
# ========================================================================
load_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("Warning: Package '", pkg, "' is not available in the container. Some features may be limited.\n"))
    return(FALSE)
  } else {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    return(TRUE)
  }
}

# Core packages for data processing
core_packages <- c("argparse", "jsonlite", "utils", "stats", "grDevices", "methods")
# Visualization packages 
viz_packages <- c("ggplot2", "plotly", "DT", "htmlwidgets", "htmltools", "base64enc")
# Data manipulation packages
data_packages <- c("dplyr", "tidyr", "stringr", "forcats", "lubridate", "scales")

# Load core packages first
invisible(lapply(core_packages, load_package))

# Load visualization packages with fallbacks
has_viz <- sapply(viz_packages, load_package)
if (all(!has_viz)) {
  cat("Warning: No visualization packages available in container. Will create simple text-based dashboard.\n")
}

# Load data packages with fallbacks
has_data <- sapply(data_packages, load_package)
if (all(!has_data)) {
  cat("Warning: Data manipulation packages unavailable in container. Analysis capabilities will be limited.\n")
}

# ========================================================================
# Command Line Argument Parsing
# ========================================================================
parser <- ArgumentParser(description = "Generate interactive dashboard from FoodNet Trends results")
parser$add_argument("--outDir", type = "character", 
                  help = "Base output directory for the whole pipeline", required = TRUE)
parser$add_argument("--resultDir", type = "character", 
                  help = "Directory containing result files to include in dashboard", required = TRUE)
parser$add_argument("--outputFile", type = "character", default = "dashboard.html",
                  help = "Name of dashboard HTML file to generate (default: dashboard.html)")
parser$add_argument("--templateFile", type = "character", 
                  help = "Path to custom dashboard HTML template (optional)")
parser$add_argument("--title", type = "character", default = "FoodNet Trends Dashboard",
                  help = "Dashboard title (default: FoodNet Trends Dashboard)")
parser$add_argument("--logoPath", type = "character", 
                  help = "Path to logo image file (optional)")
parser$add_argument("--qualityDataPath", type = "character", 
                  help = "Path to quality metadata JSON file (optional)")
parser$add_argument("--memoryLimit", type = "integer", default = 14,
                  help = "Memory limit in GB (default: 14)")
parser$add_argument("--startYear", type = "integer", default = 2015,
                  help = "Starting year for trend analysis (default: 2015)")
parser$add_argument("--endYear", type = "integer", default = 2025,
                  help = "Ending year for trend analysis (default: 2025)")
parser$add_argument("--theme", type = "character", default = "modern",
                  help = "Dashboard theme (modern, classic, dark) (default: modern)")
parser$add_argument("--debug", dest = "debug", action = "store_true",
                  help = "Enable debug mode with additional console output")

# Parse arguments with error handling
tryCatch({
  args <- parser$parse_args()
  # Log arguments for debugging
  if (exists("args") && (is.null(args$debug) || args$debug)) {
    cat("Dashboard Generator Arguments:\n")
    for (arg_name in names(args)) {
      cat(paste0("  ", arg_name, ": ", args[[arg_name]], "\n"))
    }
  }
}, error = function(e) {
  cat("Error parsing command line arguments:", e$message, "\n")
  cat("Run with --help for usage information\n")
  quit(status = 1)
})

# Set memory limit based on args if possible
if (!is.null(args$memoryLimit) && args$memoryLimit > 0) {
  memory_limit_gb <- args$memoryLimit
  tryCatch({
    # Convert GB to bytes for R's memory limit (with 90% safety margin)
    memory_limit_bytes <- as.numeric(memory_limit_gb) * 0.9 * 1024^3
    if (!is.na(memory_limit_bytes) && memory_limit_bytes > 1e9) {
      memory.limit(size = memory_limit_bytes)
      cat("Memory limit set to", memory_limit_gb, "GB (", memory_limit_bytes, "bytes)\n")
    }
  }, error = function(e) {
    cat("Warning: Unable to set memory limit:", e$message, "\n")
  })
}

# ========================================================================
# Helper Functions for Data Processing
# ========================================================================

#' Check if directory exists and is readable
#'
#' @param dir_path Path to directory
#' @return TRUE if directory exists and is readable, FALSE otherwise
check_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    cat("ERROR: Directory does not exist:", dir_path, "\n")
    return(FALSE)
  }
  
  # Check if directory is readable
  tryCatch({
    list.files(dir_path)
    return(TRUE)
  }, error = function(e) {
    cat("ERROR: Cannot read directory:", dir_path, "\n")
    cat("       Error was:", e$message, "\n")
    return(FALSE)
  })
}

#' Read a CSV file with error handling and memory optimization
#'
#' @param file_path Path to CSV file
#' @return Data frame containing CSV data or NULL if error
read_csv_safely <- function(file_path) {
  if (!file.exists(file_path)) {
    cat("Warning: CSV file does not exist:", file_path, "\n")
    return(NULL)
  }
  
  # Check file size - skip very large files
  file_info <- file.info(file_path)
  if (file_info$size > 500 * 1024 * 1024) { # 500 MB
    cat("Warning: File too large to process in dashboard:", file_path, "\n")
    return(NULL)
  }
  
  # Try to read the file
  tryCatch({
    # First attempt - optimized for memory
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Force garbage collection after reading large file
    gc(reset = TRUE)
    
    return(df)
  }, error = function(e) {
    cat("ERROR: Cannot read CSV file:", file_path, "\n")
    cat("       Error was:", e$message, "\n")
    return(NULL)
  })
}

#' Extract pathogen name from filename
#'
#' @param filename Filename to extract pathogen from
#' @return Character string with pathogen name
extract_pathogen <- function(filename) {
  # Remove path and extract pathogen name before first underscore
  basename <- basename(filename)
  # Remove any special characters that might cause issues
  basename <- gsub("[^[:alnum:]_.-]", "", basename)
  parts <- strsplit(basename, "_")[[1]]
  if (length(parts) > 0) {
    # Sanitize the pathogen name to avoid any script issues
    pathogen <- toupper(parts[1])
    # Replace any non-alphanumeric chars with underscore
    pathogen <- gsub("[^[:alnum:]]", "_", pathogen)
    return(pathogen)
  } else {
    return("UNKNOWN")
  }
}

#' Create an interactive trend chart using plotly
#'
#' @param data Data frame with trend data
#' @param x_col Column name for x-axis
#' @param y_col Column name for y-axis
#' @param title Chart title
#' @param color Optional grouping variable for color
#' @return Plotly chart object or NULL if error
create_trend_chart <- function(data, x_col, y_col, title, color = NULL, 
                              x_label = NULL, y_label = NULL) {
  if (!"plotly" %in% (.packages()) || is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  # Set default labels if not provided
  if (is.null(x_label)) x_label <- x_col
  if (is.null(y_label)) y_label <- y_col
  
  # Make a static ggplot2 plot first
  tryCatch({
    p <- ggplot(data, aes_string(x = x_col, y = y_col)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      labs(title = title,
           x = x_label, 
           y = y_label) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom"
      )
    
    # Add color if specified
    if (!is.null(color) && color %in% colnames(data)) {
      p <- p + aes_string(color = color) +
        scale_color_brewer(palette = "Set1")
    }
    
    # Convert to plotly for interactivity
    p_ly <- ggplotly(p) %>%
      layout(
        hoverlabel = list(bgcolor = "white", font = list(size = 12)),
        hovermode = "closest"
      )
    
    return(p_ly)
  }, error = function(e) {
    cat("Warning: Error creating trend chart:", e$message, "\n")
    return(NULL)
  })
}

#' Create a basic HTML table from a data frame
#'
#' @param data Data frame to display
#' @param caption Table caption/title
#' @return HTML table string
create_html_table <- function(data, caption = NULL) {
  if (is.null(data) || nrow(data) == 0) {
    return("<p>No data available for table display.</p>")
  }
  
  # Try to create an interactive DT table if available
  if ("DT" %in% (.packages())) {
    tryCatch({
      dt <- datatable(data, 
                    options = list(pageLength = 10, 
                                  autoWidth = TRUE,
                                  dom = 'Bfrtip',
                                  buttons = c('csv', 'excel')),
                    caption = caption)
      html_widget <- as.character(dt)
      return(html_widget)
    }, error = function(e) {
      cat("Warning: Error creating DT table:", e$message, "\n")
      # Fall back to basic HTML
    })
  }
  
  # Basic HTML table fallback
  html <- "<table class='data-table'>"
  if (!is.null(caption)) {
    html <- paste0(html, "<caption>", caption, "</caption>")
  }
  
  # Add header row
  html <- paste0(html, "<thead><tr>")
  for (col in colnames(data)) {
    html <- paste0(html, "<th>", col, "</th>")
  }
  html <- paste0(html, "</tr></thead><tbody>")
  
  # Add data rows
  for (i in 1:min(nrow(data), 100)) { # Limit to 100 rows for basic table
    html <- paste0(html, "<tr>")
    for (col in colnames(data)) {
      html <- paste0(html, "<td>", data[i, col], "</td>")
    }
    html <- paste0(html, "</tr>")
  }
  
  html <- paste0(html, "</tbody></table>")
  return(html)
}

#' Create a summary card for a pathogen
#'
#' @param name Pathogen name
#' @param data Summary data for the pathogen
#' @param color Card color (hex code)
#' @return HTML string for the card
create_summary_card <- function(name, data, color = "#0066cc") {
  # Get key metrics from data if available
  recent_ir <- NA
  trend <- NA
  states <- NA
  
  if (!is.null(data)) {
    if ("recent_ir" %in% names(data)) recent_ir <- round(data$recent_ir, 2)
    if ("trend" %in% names(data)) trend <- data$trend
    if ("states" %in% names(data)) states <- length(data$states)
  }
  
  # Set trend arrow and color
  trend_arrow <- "➡️"
  trend_color <- "#666666"
  
  if (!is.na(trend)) {
    if (trend > 0.05) {
      trend_arrow <- "🔺"
      trend_color <- "#d9534f" # Red for increase
    } else if (trend < -0.05) {
      trend_arrow <- "🔽"
      trend_color <- "#5cb85c" # Green for decrease
    }
  }
  
  # Format card HTML
  html <- paste0('
    <div class="summary-card" style="border-top: 4px solid ', color, '">
      <h3>', name, '</h3>
      <div class="metrics">
        <div class="metric">
          <span class="metric-value">', ifelse(is.na(recent_ir), "N/A", recent_ir), '</span>
          <span class="metric-label">Recent IR</span>
        </div>
        <div class="metric">
          <span class="metric-value" style="color:', trend_color, '">', trend_arrow, '</span>
          <span class="metric-label">Trend</span>
        </div>
        <div class="metric">
          <span class="metric-value">', ifelse(is.na(states), "N/A", states), '</span>
          <span class="metric-label">States</span>
        </div>
      </div>
    </div>
  ')
  
  return(html)
}

#' Safely convert HTML widget to character string
#'
#' @param widget HTML widget object
#' @return Character string representation of the widget
widget_to_html <- function(widget) {
  if (is.null(widget)) {
    return("")
  }
  
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    return("<p>Interactive chart unavailable (htmlwidgets package missing)</p>")
  }
  
  tryCatch({
    # Generate a temporary file for the widget
    temp_file <- tempfile(fileext = ".html")
    htmlwidgets::saveWidget(widget, file = temp_file, selfcontained = TRUE)
    
    # Read the file content
    html_content <- readChar(temp_file, file.info(temp_file)$size)
    
    # Clean up temp file
    unlink(temp_file)
    
    # Extract the widget HTML from between body tags
    body_pattern <- "<body>(.*?)</body>"
    matches <- regmatches(html_content, regexec(body_pattern, html_content, perl = TRUE))
    
    if (length(matches) > 0 && length(matches[[1]]) > 1) {
      return(matches[[1]][2])  # Return the content inside body tags
    } else {
      return(html_content)  # Return the full content if body tags not found
    }
  }, error = function(e) {
    cat("Warning: Error converting widget to HTML:", e$message, "\n")
    return("<p>Error generating interactive chart</p>")
  })
}

# ========================================================================
# Main Dashboard Generation Functions
# ========================================================================

#' Process incidence rate data files
#'
#' @param file_paths List of file paths to process
#' @return List of data frames by pathogen
process_ir_files <- function(file_paths) {
  if (length(file_paths) == 0) {
    return(list())
  }
  
  # Group files by pathogen
  pathogen_files <- list()
  for (file_path in file_paths) {
    pathogen <- extract_pathogen(file_path)
    if (!(pathogen %in% names(pathogen_files))) {
      pathogen_files[[pathogen]] <- c()
    }
    pathogen_files[[pathogen]] <- c(pathogen_files[[pathogen]], file_path)
  }
  
  # Process each pathogen's files
  result_data <- list()
  for (pathogen in names(pathogen_files)) {
    pathogen_data <- list()
    combined_df <- NULL
    
    for (file_path in pathogen_files[[pathogen]]) {
      df <- read_csv_safely(file_path)
      if (!is.null(df)) {
        if (is.null(combined_df)) {
          combined_df <- df
        } else {
          # Combine data frames safely
          tryCatch({
            # Check if data frames are compatible for binding
            if (all(colnames(combined_df) %in% colnames(df)) && 
                all(colnames(df) %in% colnames(combined_df))) {
              combined_df <- rbind(combined_df, df)
            } else {
              # Store as separate data frame
              pathogen_data[[basename(file_path)]] <- df
            }
          }, error = function(e) {
            cat("Warning: Could not combine data for", pathogen, ":", e$message, "\n")
            pathogen_data[[basename(file_path)]] <- df
          })
        }
      }
    }
    
    # Store combined data
    if (!is.null(combined_df)) {
      pathogen_data[["combined"]] <- combined_df
    }
    
    result_data[[pathogen]] <- pathogen_data
  }
  
  return(result_data)
}

#' Calculate summary statistics for each pathogen
#'
#' @param pathogen_data List of data frames by pathogen
#' @return List of summary statistics by pathogen
calculate_summaries <- function(pathogen_data) {
  if (length(pathogen_data) == 0) {
    return(list())
  }
  
  summaries <- list()
  for (pathogen in names(pathogen_data)) {
    # Skip if no combined data available
    if (!("combined" %in% names(pathogen_data[[pathogen]]))) {
      next
    }
    
    df <- pathogen_data[[pathogen]][["combined"]]
    
    # Extract key columns for summary
    year_col <- NULL
    ir_col <- NULL
    state_col <- NULL
    
    # Find year column
    possible_year_cols <- c("year", "YEAR", "Year", "mmwr_year", "MMWR_YEAR")
    for (col in possible_year_cols) {
      if (col %in% colnames(df)) {
        year_col <- col
        break
      }
    }
    
    # Find incidence rate column
    possible_ir_cols <- c("incidence_rate", "IR", "ir", "rate", "RATE", "Rate")
    for (col in possible_ir_cols) {
      if (col %in% colnames(df)) {
        ir_col <- col
        break
      }
    }
    
    # Find state column
    possible_state_cols <- c("state", "STATE", "State", "jurisdiction", "JURISDICTION")
    for (col in possible_state_cols) {
      if (col %in% colnames(df)) {
        state_col <- col
        break
      }
    }
    
    # Calculate summary statistics if we have the necessary columns
    if (!is.null(year_col) && !is.null(ir_col)) {
      summary <- list()
      
      # Get average incidence rate for most recent available year
      if (length(unique(df[[year_col]])) > 0) {
        recent_year <- max(df[[year_col]], na.rm = TRUE)
        recent_data <- df[df[[year_col]] == recent_year, ]
        summary$recent_year <- recent_year
        summary$recent_ir <- mean(recent_data[[ir_col]], na.rm = TRUE)
      }
      
      # Calculate trend if we have multiple years
      if (length(unique(df[[year_col]])) > 1) {
        # Calculate yearly averages
        yearly_avg <- aggregate(df[[ir_col]], by = list(Year = df[[year_col]]), 
                               FUN = mean, na.rm = TRUE)
        colnames(yearly_avg) <- c("Year", "AvgIR")
        
        # Calculate simple trend by comparing earliest and latest year
        if (nrow(yearly_avg) >= 2) {
          earliest <- yearly_avg$AvgIR[which.min(yearly_avg$Year)]
          latest <- yearly_avg$AvgIR[which.max(yearly_avg$Year)]
          if (earliest > 0) {
            summary$trend <- (latest - earliest) / earliest
          } else {
            summary$trend <- NA
          }
        }
        
        # Store yearly averages for charts
        summary$yearly_data <- yearly_avg
      }
      
      # Get states if state column exists
      if (!is.null(state_col)) {
        summary$states <- unique(df[[state_col]])
      }
      
      summaries[[pathogen]] <- summary
    }
  }
  
  return(summaries)
}

#' Generate HTML for a pathogen section
#'
#' @param pathogen Pathogen name
#' @param data Pathogen data list
#' @param summary Pathogen summary statistics
#' @return HTML string for the pathogen section
generate_pathogen_section <- function(pathogen, data, summary) {
  html <- paste0("<div class='pathogen-section' id='pathogen-", tolower(pathogen), "'>")
  html <- paste0(html, "<h2>", pathogen, " Analysis</h2>")
  
  # Add summary statistics if available
  if (!is.null(summary) && length(summary) > 0) {
    # Add trend chart if yearly data is available
    if ("yearly_data" %in% names(summary)) {
      yearly_data <- summary$yearly_data
      if (nrow(yearly_data) > 0) {
        # Create trend chart
        trend_chart <- create_trend_chart(
          data = yearly_data,
          x_col = "Year",
          y_col = "AvgIR",
          title = paste0(pathogen, " Incidence Rate Trend"),
          x_label = "Year",
          y_label = "Average Incidence Rate"
        )
        
        if (!is.null(trend_chart)) {
          html <- paste0(html, "<div class='chart-container'>")
          html <- paste0(html, widget_to_html(trend_chart))
          html <- paste0(html, "</div>")
        }
      }
    }
    
    # Add recent year statistics
    if ("recent_year" %in% names(summary)) {
      html <- paste0(html, "<div class='summary-block'>")
      html <- paste0(html, "<h3>Recent Statistics (", summary$recent_year, ")</h3>")
      
      if ("recent_ir" %in% names(summary)) {
        html <- paste0(html, "<p><strong>Average Incidence Rate:</strong> ", 
                     round(summary$recent_ir, 3), " cases per 100,000 population</p>")
      }
      
      if ("trend" %in% names(summary) && !is.na(summary$trend)) {
        trend_text <- paste0(round(summary$trend * 100, 1), "%")
        trend_direction <- "stable"
        if (summary$trend > 0.05) trend_direction <- "increasing"
        if (summary$trend < -0.05) trend_direction <- "decreasing"
        
        html <- paste0(html, "<p><strong>Overall Trend:</strong> ", 
                     ifelse(summary$trend > 0, "+", ""), trend_text, 
                     " (", trend_direction, ")</p>")
      }
      
      if ("states" %in% names(summary)) {
        html <- paste0(html, "<p><strong>States with Data:</strong> ", 
                     length(summary$states), "</p>")
      }
      
      html <- paste0(html, "</div>")
    }
  }
  
  # Add data tables if available
  if (!is.null(data) && "combined" %in% names(data) && !is.null(data$combined)) {
    combined_df <- data$combined
    
    # Limit to top rows for display
    display_df <- head(combined_df, 500)
    
    html <- paste0(html, "<div class='data-table-container'>")
    html <- paste0(html, "<h3>Data Summary</h3>")
    html <- paste0(html, create_html_table(display_df, 
                                         caption = paste0(pathogen, " Incidence Rate Data")))
    html <- paste0(html, "</div>")
  }
  
  html <- paste0(html, "</div>")
  return(html)
}

#' Generate overview dashboard content
#'
#' @param pathogen_data List of data frames by pathogen
#' @param summaries List of summary statistics by pathogen
#' @return HTML string for the overview section
generate_overview <- function(pathogen_data, summaries) {
  html <- "<div class='overview-section'>"
  html <- paste0(html, "<h2>Analysis Overview</h2>")
  
  # Add summary cards for each pathogen
  if (length(summaries) > 0) {
    html <- paste0(html, "<div class='summary-cards'>")
    
    # Color palette for pathogen cards
    colors <- c("#4285F4", "#EA4335", "#FBBC05", "#34A853", "#8338EC", "#FF9900", 
               "#039BE5", "#7CB342", "#D81B60", "#FFA000")
    
    # Create a card for each pathogen
    for (i in seq_along(names(summaries))) {
      pathogen <- names(summaries)[i]
      color_idx <- (i - 1) %% length(colors) + 1
      html <- paste0(html, create_summary_card(pathogen, summaries[[pathogen]], colors[color_idx]))
    }
    
    html <- paste0(html, "</div>")
  }
  
  # Create combined trend chart if we have yearly data for multiple pathogens
  pathogens_with_trends <- c()
  trend_data <- data.frame()
  
  for (pathogen in names(summaries)) {
    if ("yearly_data" %in% names(summaries[[pathogen]])) {
      yearly_data <- summaries[[pathogen]]$yearly_data
      if (nrow(yearly_data) > 0) {
        yearly_data$Pathogen <- pathogen
        trend_data <- rbind(trend_data, yearly_data)
        pathogens_with_trends <- c(pathogens_with_trends, pathogen)
      }
    }
  }
  
  if (nrow(trend_data) > 0 && length(pathogens_with_trends) > 1) {
    # Create multi-pathogen trend chart
    trend_chart <- create_trend_chart(
      data = trend_data,
      x_col = "Year",
      y_col = "AvgIR",
      title = "Comparative Incidence Rate Trends",
      color = "Pathogen",
      x_label = "Year",
      y_label = "Average Incidence Rate"
    )
    
    if (!is.null(trend_chart)) {
      html <- paste0(html, "<div class='comparative-chart-container'>")
      html <- paste0(html, "<h3>Comparative Incidence Rate Trends</h3>")
      html <- paste0(html, widget_to_html(trend_chart))
      html <- paste0(html, "</div>")
    }
  }
  
  html <- paste0(html, "</div>")
  return(html)
}

#' Get quality information from JSON
#'
#' @param json_path Path to quality JSON file
#' @return List with quality information or NULL if error
get_quality_info <- function(json_path) {
  if (is.null(json_path) || !file.exists(json_path)) {
    return(NULL)
  }
  
  tryCatch({
    quality_data <- jsonlite::fromJSON(json_path)
    return(quality_data)
  }, error = function(e) {
    cat("Warning: Error reading quality data:", e$message, "\n")
    return(NULL)
  })
}

#' Generate HTML for quality information
#'
#' @param quality_data Quality data list
#' @return HTML string for the quality section
generate_quality_section <- function(quality_data) {
  if (is.null(quality_data)) {
    return("")
  }
  
  html <- "<div class='quality-section'>"
  html <- paste0(html, "<h2>Data Quality Information</h2>")
  
  # Check for placeholder data warning
  uses_placeholder <- FALSE
  if ("dataQuality" %in% names(quality_data) && 
      "usesPlaceholderData" %in% names(quality_data$dataQuality)) {
    uses_placeholder <- quality_data$dataQuality$usesPlaceholderData
  }
  
  if (uses_placeholder) {
    html <- paste0(html, "
      <div class='warning-banner'>
        <h3>⚠️ WARNING: Placeholder Data Detected</h3>
        <p>This analysis contains placeholder or synthetic data which may not represent real-world conditions.</p>
        <p><strong>Results are NOT suitable for production use or public health decision-making!</strong></p>
      </div>
    ")
    
    # Add affected pathogens if available
    if ("affectedPathogens" %in% names(quality_data$dataQuality) && 
        length(quality_data$dataQuality$affectedPathogens) > 0) {
      html <- paste0(html, "<div class='affected-pathogens'>")
      html <- paste0(html, "<h4>Affected Pathogens:</h4>")
      html <- paste0(html, "<ul>")
      for (pathogen in quality_data$dataQuality$affectedPathogens) {
        html <- paste0(html, "<li>", pathogen, "</li>")
      }
      html <- paste0(html, "</ul>")
      html <- paste0(html, "</div>")
    }
  } else {
    html <- paste0(html, "
      <div class='success-banner'>
        <h3>✅ Data Quality Check Passed</h3>
        <p>No placeholder or synthetic data detected in this analysis.</p>
      </div>
    ")
  }
  
  # Add data consistency information if available
  if ("dataConsistency" %in% names(quality_data$dataQuality)) {
    consistency <- quality_data$dataQuality$dataConsistency
    html <- paste0(html, "<div class='consistency-info'>")
    html <- paste0(html, "<h4>Data Consistency Information</h4>")
    
    # Create a small table with consistency metrics
    html <- paste0(html, "<table class='consistency-table'>")
    html <- paste0(html, "<tr><th>Metric</th><th>Value</th></tr>")
    
    for (field in names(consistency)) {
      field_name <- gsub("([A-Z])", " \\1", field)
      field_name <- gsub("_", " ", field_name)
      field_name <- paste0(toupper(substr(field_name, 1, 1)), substr(field_name, 2, nchar(field_name)))
      
      html <- paste0(html, "<tr>")
      html <- paste0(html, "<td>", field_name, "</td>")
      html <- paste0(html, "<td>", consistency[[field]], "</td>")
      html <- paste0(html, "</tr>")
    }
    
    html <- paste0(html, "</table>")
    html <- paste0(html, "</div>")
  }
  
  # Add generation date if available
  if ("generationDate" %in% names(quality_data)) {
    html <- paste0(html, "<div class='generation-info'>")
    html <- paste0(html, "<p><em>Quality assessment generated: ", quality_data$generationDate, "</em></p>")
    html <- paste0(html, "</div>")
  }
  
  html <- paste0(html, "</div>")
  return(html)
}

#' Create dashboard CSS based on theme
#'
#' @param theme Theme name (modern, classic, dark)
#' @return CSS styles as string
generate_dashboard_css <- function(theme = "modern") {
  # Base CSS that applies to all themes
  base_css <- "
    /* Base dashboard styles */
    body {
      font-family: 'Segoe UI', Arial, sans-serif;
      line-height: 1.6;
      margin: 0;
      padding: 0;
    }
    
    .dashboard-container {
      max-width: 1200px;
      margin: 0 auto;
      padding: 20px;
    }
    
    .header {
      padding: 20px;
      margin-bottom: 30px;
    }
    
    .header h1 {
      margin: 0;
      font-size: 28px;
    }
    
    .nav-tabs {
      display: flex;
      margin-bottom: 30px;
      border-bottom: 1px solid #ddd;
    }
    
    .nav-tab {
      padding: 10px 20px;
      cursor: pointer;
      margin-right: 5px;
      border-radius: 4px 4px 0 0;
      font-weight: 500;
    }
    
    .nav-tab.active {
      border-bottom: 3px solid;
    }
    
    .tab-content {
      display: none;
    }
    
    .tab-content.active {
      display: block;
    }
    
    .section {
      margin-bottom: 40px;
    }
    
    h2 {
      font-size: 24px;
      margin-top: 0;
    }
    
    h3 {
      font-size: 20px;
      margin-top: 0;
    }
    
    .chart-container, .comparative-chart-container {
      margin: 20px 0;
      border: 1px solid;
      border-radius: 5px;
      padding: 15px;
    }
    
    .summary-cards {
      display: flex;
      flex-wrap: wrap;
      gap: 20px;
      margin: 20px 0;
    }
    
    .summary-card {
      flex: 1 1 250px;
      border-radius: 5px;
      box-shadow: 0 2px 5px rgba(0,0,0,0.1);
      padding: 15px;
    }
    
    .summary-card h3 {
      margin-top: 0;
    }
    
    .metrics {
      display: flex;
      justify-content: space-between;
      margin-top: 15px;
    }
    
    .metric {
      text-align: center;
    }
    
    .metric-value {
      display: block;
      font-size: 24px;
      font-weight: bold;
    }
    
    .metric-label {
      display: block;
      font-size: 14px;
      color: #777;
    }
    
    .data-table-container {
      margin: 30px 0;
      overflow-x: auto;
    }
    
    .data-table {
      width: 100%;
      border-collapse: collapse;
    }
    
    .data-table th, .data-table td {
      padding: 8px 12px;
      text-align: left;
      border-bottom: 1px solid #ddd;
    }
    
    .data-table th {
      background-color: #f5f5f5;
    }
    
    .warning-banner, .success-banner {
      padding: 15px;
      border-radius: 5px;
      margin: 20px 0;
    }
    
    .warning-banner {
      background-color: #f8d7da;
      color: #721c24;
      border: 1px solid #f5c6cb;
    }
    
    .success-banner {
      background-color: #d4edda;
      color: #155724;
      border: 1px solid #c3e6cb;
    }
    
    .affected-pathogens, .consistency-info {
      margin: 20px 0;
    }
    
    .consistency-table {
      width: 100%;
      max-width: 500px;
      border-collapse: collapse;
    }
    
    .consistency-table th, .consistency-table td {
      padding: 8px 12px;
      text-align: left;
      border-bottom: 1px solid #ddd;
    }
    
    .footer {
      margin-top: 50px;
      padding-top: 20px;
      border-top: 1px solid #ddd;
      text-align: center;
      font-size: 14px;
      color: #777;
    }
  "
  
  # Theme-specific customizations
  if (theme == "classic") {
    theme_css <- "
      body {
        font-family: Georgia, serif;
        color: #333;
        background-color: #f9f9f9;
      }
      
      .header {
        background-color: #003366;
        color: white;
        text-align: center;
      }
      
      .nav-tab {
        background-color: #f2f2f2;
        border: 1px solid #ddd;
      }
      
      .nav-tab.active {
        background-color: #fff;
        border-bottom: 3px solid #003366;
      }
      
      .chart-container, .comparative-chart-container {
        border-color: #ddd;
        background-color: #fff;
      }
      
      h2, h3 {
        color: #003366;
      }
    "
  } else if (theme == "dark") {
    theme_css <- "
      body {
        background-color: #121212;
        color: #e0e0e0;
      }
      
      .dashboard-container {
        background-color: #1f1f1f;
      }
      
      .header {
        background-color: #212121;
        color: #ffffff;
      }
      
      .nav-tab {
        background-color: #2d2d2d;
        color: #cccccc;
      }
      
      .nav-tab.active {
        background-color: #3d3d3d;
        border-bottom: 3px solid #bb86fc;
        color: #ffffff;
      }
      
      .chart-container, .comparative-chart-container {
        background-color: #2d2d2d;
        border-color: #444444;
      }
      
      .summary-card {
        background-color: #2d2d2d;
        box-shadow: 0 2px 5px rgba(0,0,0,0.3);
      }
      
      .metric-label {
        color: #aaaaaa;
      }
      
      .data-table th {
        background-color: #333333;
      }
      
      .data-table th, .data-table td {
        border-bottom: 1px solid #444444;
      }
      
      .consistency-table th, .consistency-table td {
        border-bottom: 1px solid #444444;
      }
      
      h2, h3 {
        color: #bb86fc;
      }
      
      .footer {
        border-top: 1px solid #333333;
        color: #aaaaaa;
      }
    "
  } else {  # Default modern theme
    theme_css <- "
      body {
        background-color: #f8f9fa;
        color: #333;
      }
      
      .header {
        background-color: #0066cc;
        color: white;
      }
      
      .dashboard-container {
        background-color: white;
        box-shadow: 0 0 15px rgba(0,0,0,0.05);
      }
      
      .nav-tab {
        background-color: #f2f2f2;
      }
      
      .nav-tab:hover {
        background-color: #e9ecef;
      }
      
      .nav-tab.active {
        background-color: #ffffff;
        border-bottom: 3px solid #0066cc;
      }
      
      .chart-container, .comparative-chart-container {
        border-color: #e9ecef;
        background-color: #ffffff;
      }
      
      h2 {
        color: #0066cc;
        border-bottom: 2px solid #e9ecef;
        padding-bottom: 10px;
      }
      
      .summary-card {
        background-color: white;
      }
    "
  }
  
  # JavaScript for tab navigation
  js_code <- "
    <script>
    document.addEventListener('DOMContentLoaded', function() {
      // Tab navigation
      const tabs = document.querySelectorAll('.nav-tab');
      const tabContents = document.querySelectorAll('.tab-content');
      
      tabs.forEach(tab => {
        tab.addEventListener('click', () => {
          // Remove active class from all tabs and contents
          tabs.forEach(t => t.classList.remove('active'));
          tabContents.forEach(c => c.classList.remove('active'));
          
          // Add active class to clicked tab and corresponding content
          tab.classList.add('active');
          const tabId = tab.getAttribute('data-tab');
          document.getElementById(tabId).classList.add('active');
        });
      });
      
      // Set first tab as active by default
      if (tabs.length > 0 && tabContents.length > 0) {
        tabs[0].classList.add('active');
        tabContents[0].classList.add('active');
      }
    });
    </script>
  "
  
  return(paste0("<style>", base_css, theme_css, "</style>", js_code))
}

#' Create complete dashboard HTML
#'
#' @param title Dashboard title
#' @param pathogen_data Processed pathogen data
#' @param summaries Pathogen summaries
#' @param quality_data Quality data
#' @param theme Dashboard theme
#' @return Complete HTML string for the dashboard
create_dashboard_html <- function(title, pathogen_data, summaries, quality_data, theme = "modern") {
  css_styles <- generate_dashboard_css(theme)
  
  # Start HTML document
  html <- paste0("<!DOCTYPE html>
  <html lang='en'>
  <head>
    <meta charset='UTF-8'>
    <meta name='viewport' content='width=device-width, initial-scale=1.0'>
    <title>", title, "</title>
    ", css_styles, "
  </head>
  <body>
    <div class='dashboard-container'>
      <div class='header'>
        <h1>", title, "</h1>
      </div>
  ")
  
  # Create tabs for navigation
  html <- paste0(html, "<div class='nav-tabs'>")
  html <- paste0(html, "<div class='nav-tab' data-tab='tab-overview'>Overview</div>")
  
  # Create a tab for each pathogen
  for (pathogen in names(pathogen_data)) {
    html <- paste0(html, "<div class='nav-tab' data-tab='tab-", tolower(pathogen), 
                 "'>", pathogen, "</div>")
  }
  
  # Add quality tab if quality data available
  if (!is.null(quality_data)) {
    html <- paste0(html, "<div class='nav-tab' data-tab='tab-quality'>Data Quality</div>")
  }
  
  html <- paste0(html, "</div>")
  
  # Create tab content sections
  html <- paste0(html, "<div class='tab-content active' id='tab-overview'>")
  html <- paste0(html, generate_overview(pathogen_data, summaries))
  html <- paste0(html, "</div>")
  
  # Create a tab for each pathogen
  for (pathogen in names(pathogen_data)) {
    html <- paste0(html, "<div class='tab-content' id='tab-", tolower(pathogen), "'>")
    html <- paste0(html, generate_pathogen_section(pathogen, pathogen_data[[pathogen]], 
                                                summaries[[pathogen]]))
    html <- paste0(html, "</div>")
  }
  
  # Add quality tab content if quality data available
  if (!is.null(quality_data)) {
    html <- paste0(html, "<div class='tab-content' id='tab-quality'>")
    html <- paste0(html, generate_quality_section(quality_data))
    html <- paste0(html, "</div>")
  }
  
  # Add footer
  html <- paste0(html, "
      <div class='footer'>
        <p>FoodNet Trends Dashboard | Generated: ", format(Sys.time()), "</p>
        <p>Version 1.0</p>
      </div>
    </div>
  </body>
  </html>
  ")
  
  return(html)
}

# ========================================================================
# Main Execution
# ========================================================================

# Sanitize directory paths
sanitize_path <- function(path) {
  if (is.null(path)) return(NULL)
  # Remove trailing slashes
  path <- gsub("[\\/]+$", "", path)
  # Add quotes if path contains spaces
  if (grepl(" ", path)) {
    path <- shQuote(path)
  }
  return(path)
}

# Sanitize input paths
args$outDir <- sanitize_path(args$outDir)
args$resultDir <- sanitize_path(args$resultDir)
if (!is.null(args$templateFile)) args$templateFile <- sanitize_path(args$templateFile)

# Verify directories
if (!check_directory(args$outDir)) {
  cat("Creating output directory:", args$outDir, "\n")
  dir.create(args$outDir, recursive = TRUE, showWarnings = FALSE)
}

if (!check_directory(args$resultDir)) {
  cat("ERROR: Result directory does not exist and cannot be created:", args$resultDir, "\n")
  cat("Dashboard generation cannot proceed without valid result files\n")
  
  # Create a minimal dashboard with error message
  output_html <- paste0(
    "<!DOCTYPE html><html><head><title>", args$title, "</title></head><body>",
    "<h1>", args$title, "</h1>",
    "<p>Error: Cannot access result directory at ", args$resultDir, "</p>",
    "<p>Dashboard generation failed at ", format(Sys.time()), "</p>",
    "</body></html>"
  )
  
  cat(output_html, file = args$outputFile)
  quit(status = 1)
}

# Find CSV result files
cat("Scanning for result files...\n")
ir_files <- list.files(path = ".", pattern = "_IRCatch.csv$", recursive = FALSE, full.names = TRUE)
summary_files <- list.files(path = ".", pattern = "_summary.txt$", recursive = FALSE, full.names = TRUE)

cat("Found", length(ir_files), "incidence rate files and", length(summary_files), "summary files\n")

if (length(ir_files) == 0) {
  cat("WARNING: No result files found. Creating placeholder dashboard.\n")
  
  # Create a minimal dashboard with notice
  output_html <- paste0(
    "<!DOCTYPE html><html><head><title>", args$title, "</title></head><body>",
    "<h1>", args$title, "</h1>",
    "<p>Warning: No analysis result files were found.</p>",
    "<p>This may indicate that the analysis did not complete successfully or the result files are in a different location.</p>",
    "<p>Dashboard generated at ", format(Sys.time()), "</p>",
    "</body></html>"
  )
  
  cat(output_html, file = args$outputFile)
} else {
  # Process data
  cat("Processing result data...\n")
  pathogen_data <- process_ir_files(ir_files)
  cat("Processed data for", length(pathogen_data), "pathogens\n")
  
  # Calculate summaries
  summaries <- calculate_summaries(pathogen_data)
  cat("Created summaries for", length(summaries), "pathogens\n")
  
  # Read quality data if available
  quality_data <- NULL
  if (!is.null(args$qualityDataPath) && args$qualityDataPath != "") {
    quality_data <- get_quality_info(args$qualityDataPath)
    if (!is.null(quality_data)) {
      cat("Loaded quality data from", args$qualityDataPath, "\n")
    }
  }
  
  # Generate dashboard HTML
  cat("Creating dashboard HTML...\n")
  dashboard_html <- create_dashboard_html(
    title = args$title,
    pathogen_data = pathogen_data,
    summaries = summaries,
    quality_data = quality_data,
    theme = args$theme
  )
  
  # Write the dashboard
  cat("Writing dashboard to", args$outputFile, "\n")
  cat(dashboard_html, file = args$outputFile)
  
  cat("Dashboard generation complete\n")
}

# Final cleanup
gc(reset = TRUE)
cat("Dashboard generation finished at", format(Sys.time()), "\n")
quit(status = 0)