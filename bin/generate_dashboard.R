#!/usr/bin/env Rscript
# =========================================================================
# FoodNet Trends v1.0 - Dashboard Generator
# =========================================================================
#
# Purpose:
#   Creates a self-contained HTML dashboard from analysis results
#
# Input:
#   - Path to results directory
#   - Optional configuration parameters
#
# Output:
#   - Self-contained HTML dashboard file with embedded data
#
# Last updated: 2025-05-18
# =========================================================================

# =========================================================================
# Load required packages with error handling
# =========================================================================
required_packages <- c("argparse", "jsonlite", "htmlwidgets", "plotly", 
                      "dplyr", "ggplot2", "DT", "htmltools", "base64enc")

for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# =========================================================================
# Parse command-line arguments
# =========================================================================
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

# Parse arguments with error handling
tryCatch({
  args <- parser$parse_args()
}, error = function(e) {
  cat("Error parsing command line arguments:", e$message, "\n")
  cat("Run with --help for usage information\n")
  quit(status = 1)
})

# =========================================================================
# Helper functions
# =========================================================================

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

#' Find all results for a given file type
#'
#' @param result_dir Directory containing result files
#' @param pattern Pattern to match (e.g., "_IRCatch.csv")
#' @return Data frame with file paths and metadata
find_result_files <- function(result_dir, pattern) {
  # Find all files matching pattern
  files <- list.files(result_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  
  if (length(files) == 0) {
    return(NULL)
  }
  
  # Extract pathogen names from file paths
  pathogen_names <- sub(paste0("^(.+)_", sub("\\.", "\\\\.", pattern), "$"), "\\1", basename(files))
  
  # Create a data frame with file information
  result <- data.frame(
    file_path = files,
    pathogen = pathogen_names,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

#' Read all incidence rate files and combine them
#'
#' @param ir_files Data frame with file paths from find_result_files()
#' @return Combined data frame with all incidence rates
read_ir_files <- function(ir_files) {
  if (is.null(ir_files) || nrow(ir_files) == 0) {
    cat("WARNING: No incidence rate files found\n")
    return(NULL)
  }
  
  # Initialize empty list for data frames
  df_list <- list()
  
  # Read each file and add to list
  for (i in 1:nrow(ir_files)) {
    tryCatch({
      df <- read.csv(ir_files$file_path[i], stringsAsFactors = FALSE)
      df_list[[i]] <- df
    }, error = function(e) {
      cat("WARNING: Failed to read file:", ir_files$file_path[i], "\n")
      cat("         Error was:", e$message, "\n")
    })
  }
  
  # Combine all data frames
  if (length(df_list) == 0) {
    return(NULL)
  }
  
  combined_df <- do.call(rbind, df_list)
  return(combined_df)
}

#' Read result summary files
#'
#' @param summary_files Data frame with file paths from find_result_files()
#' @return List of summary text content
read_summary_files <- function(summary_files) {
  if (is.null(summary_files) || nrow(summary_files) == 0) {
    cat("WARNING: No summary files found\n")
    return(NULL)
  }
  
  # Initialize empty list for summaries
  summary_list <- list()
  
  # Read each file and add to list
  for (i in 1:nrow(summary_files)) {
    tryCatch({
      text <- readLines(summary_files$file_path[i])
      summary_list[[summary_files$pathogen[i]]] <- text
    }, error = function(e) {
      cat("WARNING: Failed to read file:", summary_files$file_path[i], "\n")
      cat("         Error was:", e$message, "\n")
    })
  }
  
  return(summary_list)
}

#' Create a time series plot for incidence rates
#'
#' @param ir_data Data frame with incidence rates
#' @param pathogens Vector of pathogens to include (NULL for all)
#' @return A plotly object
create_ir_plot <- function(ir_data, pathogens = NULL) {
  if (is.null(ir_data) || nrow(ir_data) == 0) {
    # Return a simple empty plot with a message
    empty_df <- data.frame(x = 1, y = 1, label = "No data available")
    p <- ggplot(empty_df, aes(x = x, y = y, label = label)) +
      geom_text() +
      theme_void() +
      labs(title = "No Data Available")
    return(ggplotly(p))
  }
  
  # Check if we have the required columns
  required_cols <- c("year", "state", "pathogen")
  if (!all(required_cols %in% colnames(ir_data))) {
    # Create a message about missing columns
    missing_cols <- required_cols[!required_cols %in% colnames(ir_data)]
    empty_df <- data.frame(x = 1, y = 1, 
                         label = paste("Missing required columns:", 
                                     paste(missing_cols, collapse = ", ")))
    p <- ggplot(empty_df, aes(x = x, y = y, label = label)) +
      geom_text() +
      theme_void() +
      labs(title = "Data Format Error")
    return(ggplotly(p))
  }
  
  # Determine which column to use for y-axis (incidence rate)
  y_col <- NULL
  if ("median_incidence" %in% colnames(ir_data)) {
    y_col <- "median_incidence"
    lower_col <- if ("lower_hdi" %in% colnames(ir_data)) "lower_hdi" else NULL
    upper_col <- if ("upper_hdi" %in% colnames(ir_data)) "upper_hdi" else NULL
  } else if ("ir" %in% colnames(ir_data)) {
    y_col <- "ir"
    lower_col <- if ("ir_lower" %in% colnames(ir_data)) "ir_lower" else NULL
    upper_col <- if ("ir_upper" %in% colnames(ir_data)) "ir_upper" else NULL
  } else {
    # No recognizable incidence rate column
    empty_df <- data.frame(x = 1, y = 1, 
                         label = "No recognized incidence rate column found")
    p <- ggplot(empty_df, aes(x = x, y = y, label = label)) +
      geom_text() +
      theme_void() +
      labs(title = "Data Format Error")
    return(ggplotly(p))
  }
  
  # Filter by pathogens if provided
  if (!is.null(pathogens) && length(pathogens) > 0) {
    ir_data <- ir_data[ir_data$pathogen %in% pathogens, ]
  }
  
  # If after filtering we have no data, return empty plot
  if (nrow(ir_data) == 0) {
    empty_df <- data.frame(x = 1, y = 1, label = "No data available for selected pathogens")
    p <- ggplot(empty_df, aes(x = x, y = y, label = label)) +
      geom_text() +
      theme_void() +
      labs(title = "No Data Available")
    return(ggplotly(p))
  }
  
  # Get unique pathogens for color scale
  unique_pathogens <- unique(ir_data$pathogen)
  
  # Create tooltip text based on available columns
  if (!is.null(lower_col) && !is.null(upper_col)) {
    tooltip_text <- paste0("State: ", ir_data$state, 
                          "<br>Year: ", ir_data$year,
                          "<br>Pathogen: ", ir_data$pathogen,
                          "<br>Incidence: ", round(ir_data[[y_col]], 2),
                          "<br>95% CI: ", round(ir_data[[lower_col]], 2), " - ", 
                          round(ir_data[[upper_col]], 2))
  } else {
    tooltip_text <- paste0("State: ", ir_data$state, 
                          "<br>Year: ", ir_data$year,
                          "<br>Pathogen: ", ir_data$pathogen,
                          "<br>Incidence: ", round(ir_data[[y_col]], 2))
  }
  
  # Create base plot with ggplot2
  p <- ggplot(ir_data, aes_string(x = "year", y = y_col, color = "pathogen", text = "tooltip_text")) +
    geom_line(aes(group = interaction(pathogen, state)), alpha = 0.5) +
    geom_point(size = 1) +
    facet_wrap(~ state, scales = "free_y") +
    labs(title = "Incidence Rates by State and Pathogen",
         subtitle = "Incidence per 100,000 population",
         x = "Year",
         y = "Incidence Rate",
         color = "Pathogen") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  # Convert to plotly for interactivity
  plt <- ggplotly(p, tooltip = "text")
  
  # Add custom hover behavior
  plt <- plt %>% layout(
    hovermode = "closest",
    legend = list(orientation = "h", y = 1.1),
    margin = list(t = 100)
  )
  
  return(plt)
}

#' Create a map visualization of incidence rates
#'
#' @param ir_data Data frame with incidence rates
#' @param year_selected Selected year for map
#' @param pathogen_selected Selected pathogen for map
#' @return A plotly map object
create_map_plot <- function(ir_data, year_selected, pathogen_selected) {
  if (is.null(ir_data)) {
    return(NULL)
  }
  
  # Filter data for selected year and pathogen
  map_data <- ir_data[ir_data$year == year_selected & ir_data$pathogen == pathogen_selected, ]
  
  # Create a basic US map plot
  # This is a simplified version - in production, use proper US state boundaries and geojson
  # For now, create a placeholder that would be replaced with actual map
  p <- plot_ly(map_data, 
              type = "choropleth",
              z = ~median_incidence,
              text = ~paste("State:", state, 
                           "<br>Incidence:", round(median_incidence, 2),
                           "<br>95% CI:", round(lower_hdi, 2), "-", round(upper_hdi, 2)),
              colorscale = "YlOrRd",
              marker = list(line = list(color = "rgb(255,255,255)", width = 1))) %>%
    layout(
      title = paste("Incidence Rates for", pathogen_selected, "in", year_selected),
      geo = list(scope = "usa")
    )
  
  return(p)
}

#' Create a comparison plot for relative risk
#'
#' @param irr_data Data frame with relative risk data
#' @param comparison_period Selected comparison period
#' @return A plotly object
create_rr_plot <- function(irr_data, comparison_period) {
  if (is.null(irr_data)) {
    return(NULL)
  }
  
  # Filter for selected comparison period
  rr_data <- irr_data[irr_data$comparison_period == comparison_period, ]
  
  # Create plot
  p <- ggplot(rr_data, aes(x = reorder(state, relative_risk), 
                          y = relative_risk, 
                          fill = relative_risk > 1,
                          text = paste("State:", state,
                                      "<br>Relative Risk:", round(relative_risk, 2),
                                      "<br>Current:", round(current_incidence, 2),
                                      "<br>Baseline:", round(period_incidence, 2)))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    facet_wrap(~ pathogen, scales = "free_x") +
    scale_fill_manual(values = c("blue", "red"), guide = "none") +
    labs(title = paste("Relative Risk Compared to", comparison_period),
         x = "State",
         y = "Relative Risk (Current / Baseline)") +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Convert to plotly
  plt <- ggplotly(p, tooltip = "text")
  
  return(plt)
}

#' Embed an image file directly in HTML with base64 encoding
#'
#' @param img_path Path to image file
#' @return HTML-ready string with base64-encoded image
embed_image <- function(img_path) {
  if (is.null(img_path) || !file.exists(img_path)) {
    return(NULL)
  }
  
  # Get file extension
  ext <- tolower(tools::file_ext(img_path))
  
  # Define MIME type based on extension
  mime_type <- switch(ext,
                     "png" = "image/png",
                     "jpg" = "image/jpeg",
                     "jpeg" = "image/jpeg",
                     "gif" = "image/gif",
                     "svg" = "image/svg+xml",
                     "image/png")  # Default to PNG
  
  # Read and encode the image
  img_data <- readBin(img_path, "raw", file.info(img_path)$size)
  img_encoded <- base64encode(img_data)
  
  # Return data URL
  return(paste0("data:", mime_type, ";base64,", img_encoded))
}

# =========================================================================
# Main function to generate dashboard
# =========================================================================
generate_dashboard <- function() {
  cat("===============================================================\n")
  cat("FoodNet Trends Dashboard Generator\n")
  cat("===============================================================\n")
  
  # Check directories
  cat("Checking output directory:", args$outDir, "\n")
  if (!check_directory(args$outDir)) {
    stop("Invalid output directory")
  }
  
  cat("Checking results directory:", args$resultDir, "\n")
  if (!check_directory(args$resultDir)) {
    stop("Invalid results directory")
  }
  
  # Ensure output directory exists (option 2 fix)
  if (!dir.exists(args$outDir)) {
    dir.create(args$outDir, recursive = TRUE)
  }
  
  # Find pathogen result files
  cat("Finding incidence rate files...\n")
  ir_files <- find_result_files(args$resultDir, "_IRCatch.csv")
  if (is.null(ir_files)) {
    cat("WARNING: No incidence rate files found\n")
  } else {
    cat("Found", nrow(ir_files), "incidence rate files\n")
  }
  
  cat("Finding relative risk files...\n")
  rr_files <- find_result_files(args$resultDir, "_EstIRRCatch_.+\\.csv")
  if (is.null(rr_files)) {
    cat("WARNING: No relative risk files found\n")
  } else {
    cat("Found", nrow(rr_files), "relative risk files\n")
  }
  
  cat("Finding summary files...\n")
  summary_files <- find_result_files(args$resultDir, "_summary.txt")
  if (is.null(summary_files)) {
    cat("WARNING: No summary files found\n")
  } else {
    cat("Found", nrow(summary_files), "summary files\n")
  }
  
  # Read data files
  cat("Reading incidence rate data...\n")
  ir_data <- read_ir_files(ir_files)
  
  cat("Reading relative risk data...\n")
  rr_data <- read_ir_files(rr_files)
  
  cat("Reading summary data...\n")
  summary_data <- read_summary_files(summary_files)
  
  # Extract metadata for dashboard configuration
  pathogens <- c() # Initialize as empty vector
  
  # Get pathogens from data if available
  if (!is.null(ir_data) && nrow(ir_data) > 0 && "pathogen" %in% colnames(ir_data)) {
    pathogens <- unique(ir_data$pathogen)
  } 
  # Fallback: get pathogens from file names
  else if (!is.null(ir_files) && nrow(ir_files) > 0) {
    pathogens <- unique(ir_files$pathogen)
  }
  # Final fallback: check if summary files have pathogen info
  else if (!is.null(summary_files) && nrow(summary_files) > 0) {
    pathogens <- unique(summary_files$pathogen)
  }
  # Last resort: use a default
  else {
    # Check if data directory contains any pathogen-named files to guess from
    potential_pathogen_files <- list.files(args$resultDir, pattern="^[A-Z]+_.*", full.names=FALSE)
    if (length(potential_pathogen_files) > 0) {
      extracted_pathogens <- unique(sub("^([A-Z]+)_.*$", "\\1", potential_pathogen_files))
      if (length(extracted_pathogens) > 0) {
        pathogens <- extracted_pathogens
      } else {
        pathogens <- c("UNKNOWN")  # Default if nothing else works
      }
    } else {
      pathogens <- c("UNKNOWN")  # Default if nothing else works
    }
  }
  
  # Get states from data if available, or provide defaults
  states <- c()
  if (!is.null(ir_data) && nrow(ir_data) > 0 && "state" %in% colnames(ir_data)) {
    states <- unique(ir_data$state)
  } else {
    states <- c("CA", "CO", "CT", "GA", "NY")  # FoodNet default states
  }
  
  # Get years from data if available, or provide defaults
  years <- c()
  if (!is.null(ir_data) && nrow(ir_data) > 0 && "year" %in% colnames(ir_data)) {
    years <- sort(unique(ir_data$year))
  } else {
    years <- 2020:2022  # Default year range
  }
  
  # Get comparison periods from data if available
  comparison_periods <- c()
  if (!is.null(rr_data) && nrow(rr_data) > 0 && "comparison_period" %in% colnames(rr_data)) {
    comparison_periods <- unique(rr_data$comparison_period)
  }
  
  cat("Dashboard will include:\n")
  cat("  Pathogens:", paste(pathogens, collapse=", "), "\n")
  cat("  States:", paste(states, collapse=", "), "\n")
  cat("  Years:", paste(range(years), collapse=" to "), "\n")
  cat("  Comparison periods:", paste(comparison_periods, collapse=", "), "\n")
  
  # Process logo if provided
  logo_data_url <- NULL
  if (!is.null(args$logoPath) && file.exists(args$logoPath)) {
    cat("Processing logo image...\n")
    logo_data_url <- embed_image(args$logoPath)
  }
  
  # Create template variables
  template_vars <- list(
    title = args$title,
    pathogens = if (is.null(pathogens) || length(pathogens) == 0) "[]" else jsonlite::toJSON(pathogens, auto_unbox = TRUE),
    states = if (is.null(states) || length(states) == 0) "[]" else jsonlite::toJSON(states, auto_unbox = TRUE),
    years = if (is.null(years) || length(years) == 0) "[]" else jsonlite::toJSON(years, auto_unbox = TRUE),
    comparison_periods = if (is.null(comparison_periods) || length(comparison_periods) == 0) "[]" else jsonlite::toJSON(comparison_periods, auto_unbox = TRUE),
    ir_data = if (is.null(ir_data) || nrow(ir_data) == 0) "[]" else jsonlite::toJSON(ir_data, auto_unbox = TRUE),
    rr_data = if (is.null(rr_data) || nrow(rr_data) == 0) "[]" else jsonlite::toJSON(rr_data, auto_unbox = TRUE),
    summary_data = if (is.null(summary_data) || length(summary_data) == 0) "[]" else jsonlite::toJSON(summary_data, auto_unbox = TRUE),
    logo_data_url = if (is.null(logo_data_url)) "" else logo_data_url,
    generation_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  # Debug output for each variable
  cat("[DEBUG] Dashboard JSON variables to be embedded:\n")
  for (name in names(template_vars)) {
    cat(sprintf("  %s: %s\n", name, substr(template_vars[[name]], 1, 120)))
  }
  
  # Create plots
  cat("Creating visualizations...\n")
  if (!is.null(ir_data)) {
    # Create main time series plot
    time_series_plot <- create_ir_plot(ir_data)
    
    # Create map visualization for first pathogen and most recent year
    latest_year <- max(years)
    first_pathogen <- pathogens[1]
    map_plot <- create_map_plot(ir_data, latest_year, first_pathogen)
    
    # Create relative risk comparison
    if (!is.null(rr_data) && length(comparison_periods) > 0) {
      first_period <- comparison_periods[1]
      rr_plot <- create_rr_plot(rr_data, first_period)
    } else {
      rr_plot <- NULL
    }
  } else {
    time_series_plot <- NULL
    map_plot <- NULL
    rr_plot <- NULL
  }
  
  # Embed plots in HTML widgets
  dashboard_widgets <- tryCatch({
    tagList(
      # Add title and header
      tags$div(class = "dashboard-header",
              tags$h1(args$title),
              tags$p(paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))),
      
      # Add filters and controls
      tags$div(class = "dashboard-controls",
              tags$div(class = "control-group",
                      tags$label("Select Pathogens:"),
                      tags$select(id = "pathogen-select", multiple = TRUE,
                                if (length(pathogens) > 0) {
                                  lapply(pathogens, function(p) tags$option(value = p, p))
                                } else {
                                  tags$option(value = "NONE", "No pathogens detected")
                                })),
              tags$div(class = "control-group",
                      tags$label("Select States:"),
                      tags$select(id = "state-select", multiple = TRUE,
                                if (length(states) > 0) {
                                  lapply(states, function(s) tags$option(value = s, s))
                                } else {
                                  tags$option(value = "NONE", "No states detected")
                                }))),
      
      # Add time series visualization
      tags$div(class = "dashboard-widget",
              tags$h2("Incidence Rate Trends"),
              if (!is.null(time_series_plot)) {
                tryCatch({
                  as_widget(time_series_plot)
                }, error = function(e) {
                  tags$div(class = "error-message", 
                          paste("Error rendering time series plot:", e$message))
                })
              } else {
                tags$div(class = "error-message", "No incidence rate data available")
              }),
      
      # Add map visualization
      tags$div(class = "dashboard-widget",
              tags$h2("Geographic Distribution"),
              if (!is.null(map_plot)) {
                tryCatch({
                  as_widget(map_plot)
                }, error = function(e) {
                  tags$div(class = "error-message", 
                          paste("Error rendering map plot:", e$message))
                })
              } else {
                tags$div(class = "error-message", "No geographic data available")
              }),
      
      # Add relative risk comparison
      tags$div(class = "dashboard-widget",
              tags$h2("Relative Risk Comparison"),
              if (!is.null(rr_plot)) {
                tryCatch({
                  as_widget(rr_plot)
                }, error = function(e) {
                  tags$div(class = "error-message", 
                          paste("Error rendering relative risk plot:", e$message))
                })
              } else {
                tags$div(class = "error-message", "No relative risk data available")
              }),
      
      # Add data tables
      tags$div(class = "dashboard-widget",
              tags$h2("Data Tables"),
              tags$div(class = "tab-container",
                      tags$div(class = "tab-headers",
                              tags$div(class = "tab-header active", "data-tab" = "incidence", "Incidence Rates"),
                              tags$div(class = "tab-header", "data-tab" = "relative-risk", "Relative Risks")),
                      tags$div(class = "tab-content active", "data-tab" = "incidence",
                              if (!is.null(ir_data) && nrow(ir_data) > 0) {
                                tryCatch({
                                  DT::datatable(ir_data)
                                }, error = function(e) {
                                  tags$div(class = "error-message", 
                                          paste("Error rendering data table:", e$message))
                                })
                              } else {
                                tags$div(class = "error-message", "No incidence rate data available")
                              }),
                      tags$div(class = "tab-content", "data-tab" = "relative-risk",
                              if (!is.null(rr_data) && nrow(rr_data) > 0) {
                                tryCatch({
                                  DT::datatable(rr_data)
                                }, error = function(e) {
                                  tags$div(class = "error-message", 
                                          paste("Error rendering data table:", e$message))
                                })
                              } else {
                                tags$div(class = "error-message", "No relative risk data available")
                              }))),
      
      # Add footer
      tags$div(class = "dashboard-footer",
              tags$p("FoodNet Trends Analysis Pipeline v1.0"),
              tags$p("Centers for Disease Control and Prevention"))
    )
  }, error = function(e) {
    # Return a simple error message if there's any issue with widget creation
    tagList(
      tags$div(class = "dashboard-header",
              tags$h1(args$title),
              tags$p(paste("Generated (with errors):", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))),
      tags$div(class = "dashboard-widget error-message",
              tags$h2("Dashboard Generation Error"),
              tags$p(paste("Error creating dashboard widgets:", e$message)),
              tags$p("Try running the pipeline again with complete data.")),
      tags$div(class = "dashboard-footer",
              tags$p("FoodNet Trends Analysis Pipeline v1.0"),
              tags$p("Centers for Disease Control and Prevention"))
    )
  })
  
  # Create HTML widgets - with error trapping
  widget_html <- tryCatch({
    htmltools::renderTags(dashboard_widgets)$html
  }, error = function(e) {
    # Return a simple HTML error message
    paste0('<div class="dashboard-header">',
           '<h1>', args$title, '</h1>',
           '<p>Generated (with errors): ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>',
           '</div>',
           '<div class="dashboard-widget error-message">',
           '<h2>Dashboard Rendering Error</h2>',
           '<p>Error rendering HTML: ', e$message, '</p>',
           '<p>Try running the pipeline again with complete data.</p>',
           '</div>',
           '<div class="dashboard-footer">',
           '<p>FoodNet Trends Analysis Pipeline v1.0</p>',
           '<p>Centers for Disease Control and Prevention</p>',
           '</div>')
  })
  
  html_template <- gsub('<div id="dashboard-container"></div>', widget_html, html_template, fixed = TRUE)
  
  # Use custom template if provided
  if (!is.null(args$templateFile) && file.exists(args$templateFile)) {
    cat("Using custom template:", args$templateFile, "\n")
    html_template <- readLines(args$templateFile, warn = FALSE)
    html_template <- paste(html_template, collapse = "\n")
  }
  
  # Replace template variables (always, even if NULL)
  for (name in names(template_vars)) {
    placeholder2 <- paste0("{{", name, "}}")
    placeholder3 <- paste0("{{{", name, "}}}")
    value <- template_vars[[name]]
    html_template <- gsub(placeholder2, value, html_template, fixed = TRUE)
    html_template <- gsub(placeholder3, value, html_template, fixed = TRUE)
  }

  # Remove logo block from HTML template
  html_template <- gsub("{{#if logo_data_url}}.*?{{/if}}", "", html_template, perl=TRUE)

  cat("Pathogens:", paste(pathogens, collapse=", "), "\n")
  cat("States:", paste(states, collapse=", "), "\n")
  cat("Years:", paste(years, collapse=", "), "\n")
  cat("IR data rows:", if (!is.null(ir_data)) nrow(ir_data) else 0, "\n")
  cat("RR data rows:", if (!is.null(rr_data)) nrow(rr_data) else 0, "\n")
  
  # Create HTML widgets
  widget_html <- htmltools::renderTags(dashboard_widgets)$html
  html_template <- gsub('<div id="dashboard-container"></div>', widget_html, html_template, fixed = TRUE)
  
  # Write HTML to file
  output_path <- args$outputFile
  cat("Writing dashboard to:", output_path, "\n")
  writeLines(html_template, output_path)
  
  cat("Dashboard generation complete!\n")
  return(output_path)
}

# Run the dashboard generator
result <- tryCatch({
  generate_dashboard()
}, error = function(e) {
  cat("ERROR: Dashboard generation failed\n")
  cat("       ", e$message, "\n")
  quit(status = 1)
})

cat("===============================================================\n")
cat("Dashboard successfully created at:", result, "\n")
cat("Open this file in a web browser to view the dashboard\n")
cat("===============================================================\n") 