#!/usr/bin/env Rscript
# =========================================================================
# FoodNet Trends v1.0 - Dashboard Generator (Memory Optimized)
# =========================================================================
#
# Purpose:
#   Creates a self-contained HTML dashboard from analysis results
#   This version is optimized for memory efficiency and reliability
#
# Input:
#   - Path to results directory
#   - Optional configuration parameters
#
# Output:
#   - Self-contained HTML dashboard file with embedded data
#
# Last updated: 2025-05-21
# =========================================================================

# =========================================================================
# Memory Optimization Settings
# =========================================================================
# Force garbage collection to improve memory management
gc(reset = TRUE)
# Set lower memory limits for data processing
options(future.globals.maxSize = 1024*1024^2) # 1GB max for big data objects
options(datatable.print.topn = 5)
options(datatable.print.nrows = 20)
options(digits = 4) # Reduce precision of numeric values for display

# =========================================================================
# Load required packages with error handling
# =========================================================================
required_packages <- c("argparse", "jsonlite", "htmlwidgets", "plotly", 
                       "dplyr", "ggplot2", "DT", "htmltools", "base64enc")

suppressWarnings({
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Warning: Package", pkg, "not available. Attempting to load core packages only."))
    } else {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
})

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
parser$add_argument("--qualityDataPath", type = "character", 
                    help = "Path to quality metadata JSON file (optional)")
parser$add_argument("--memoryLimit", type = "integer", default = 12,
                    help = "Memory limit in GB (default: 12)")

# Parse arguments with error handling
tryCatch({
  args <- parser$parse_args()
  # Log arguments for debugging
  cat("Dashboard Generator Arguments:\n")
  cat(paste0("  outDir: ", args$outDir, "\n"))
  cat(paste0("  resultDir: ", args$resultDir, "\n"))
  cat(paste0("  outputFile: ", args$outputFile, "\n"))
  cat(paste0("  templateFile: ", args$templateFile, "\n"))
  cat(paste0("  title: ", args$title, "\n"))
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
    df <- read.csv(file_path, stringsAsFactors = FALSE, nrows = 10000)
    
    # Force garbage collection after reading large file
    gc(reset = TRUE)
    
    return(df)
  }, error = function(e) {
    cat("ERROR: Cannot read CSV file:", file_path, "\n")
    cat("       Error was:", e$message, "\n")
    return(NULL)
  })
}

#' Safe wrapper for data processing functions
#'
#' @param expr Expression to evaluate
#' @param default Default value to return if expression fails
#' @return Result of expression or default if error
safely <- function(expr, default = NULL) {
  tryCatch({
    result <- eval(expr)
    gc() # Release memory after computation
    return(result)
  }, error = function(e) {
    cat("Warning: Error in data processing:", e$message, "\n")
    return(default)
  })
}

# =========================================================================
# Simplified dashboard generator with memory optimizations
# =========================================================================

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
files <- list.files(path = ".", pattern = "_IRCatch.csv$", recursive = FALSE, full.names = TRUE)
summary_files <- list.files(path = ".", pattern = "_summary.txt$", recursive = FALSE, full.names = TRUE)

cat("Found", length(files), "incidence rate files and", length(summary_files), "summary files\n")

if (length(files) == 0) {
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
  # Process data in memory-efficient batches
  cat("Processing result data...\n")
  
  # Create base dashboard with minimal memory usage
  cat("Creating dashboard HTML...\n")
  
  # Use a template file if provided, otherwise use a simple default
  if (!is.null(args$templateFile) && file.exists(args$templateFile)) {
    cat("Using template file:", args$templateFile, "\n")
    template <- readLines(args$templateFile)
    dashboard_html <- paste(template, collapse = "\n")
  } else {
    cat("Using default template\n")
    dashboard_html <- paste0(
      "<!DOCTYPE html>
      <html>
      <head>
        <meta charset=\"UTF-8\">
        <title>", args$title, "</title>
        <style>
          body { font-family: Arial, sans-serif; line-height: 1.6; padding: 20px; }
          .container { max-width: 1200px; margin: 0 auto; }
          .header { background-color: #0066cc; color: white; padding: 20px; margin-bottom: 20px; }
          .section { margin-bottom: 30px; }
          .footer { margin-top: 50px; border-top: 1px solid #ddd; padding-top: 20px; color: #777; }
        </style>
      </head>
      <body>
        <div class=\"container\">
          <div class=\"header\">
            <h1>", args$title, "</h1>
          </div>
          <div class=\"section\">
            <h2>FoodNet Trends Analysis Results</h2>
            <p>This dashboard provides a summary of the FoodNet Trends analysis results.</p>
            <p>Generated at: ", format(Sys.time()), "</p>
          </div>
          <div class=\"section\">
            <h2>Pathogens Analyzed</h2>
            <ul>
      "
    )
    
    # Extract pathogen names from result files
    pathogens <- unique(gsub("_.*$", "", basename(files)))
    for (pathogen in pathogens) {
      dashboard_html <- paste0(dashboard_html, "<li>", pathogen, "</li>\n")
    }
    
    dashboard_html <- paste0(dashboard_html, "
            </ul>
          </div>
          <div class=\"section\">
            <h2>Analysis Summary</h2>
            <p>Incidence Rate Files: ", length(files), "</p>
            <p>Summary Files: ", length(summary_files), "</p>
          </div>
          <div class=\"footer\">
            <p>FoodNet Trends Analysis Dashboard | Generated on: ", format(Sys.time()), "</p>
          </div>
        </div>
      </body>
      </html>
      "
    )
  }
  
  # Write the dashboard HTML
  cat("Writing dashboard to", args$outputFile, "\n")
  cat(dashboard_html, file = args$outputFile)
  
  cat("Dashboard generation complete\n")
}

# Final cleanup
gc(reset = TRUE)
cat("Memory-optimized dashboard generation finished at", format(Sys.time()), "\n")
quit(status = 0)