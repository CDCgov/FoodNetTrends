#!/usr/bin/env Rscript
# =========================================================================
# FoodNetTrends v1.0 - Interactive Dashboard Generator
# =========================================================================
#
# Purpose:
#   Generates interactive HTML dashboards from completed analysis results.
#   Creates self-contained reports with embedded visualizations, data tables,
#   and responsive design for comprehensive result exploration.
#
# Key Features:
#   - Interactive image gallery with modal expansion
#   - Tabbed data interface for multi-pathogen results
#   - Dark/light mode accessibility options
#   - Responsive design for multiple device types
#   - Embedded visualizations for offline viewing
#
# Last updated: 2025-05-22
# =========================================================================

# Performance optimizations
gc(reset = TRUE)
options(
  scipen = 999,
  digits = 4,
  warn = 1
)

# ========================================================================
# Load Essential Packages Only
# ========================================================================
load_package <- function(pkg) {
  suppressWarnings(suppressMessages({
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE)
      return(TRUE)
    }
    return(FALSE)
  }))
}

# Core packages only - no heavy dependencies
essential_packages <- c("argparse", "jsonlite", "utils")
loaded <- sapply(essential_packages, load_package)

# Load base64enc package if available, use fallback encoding otherwise
has_base64enc <- load_package("base64enc")

if (!all(loaded)) {
  cat("Error: Required packages not available in container\n")
  quit(status = 1)
}

# ========================================================================
# Command Line Arguments
# ========================================================================
parser <- ArgumentParser(description = "Generate clean performance dashboard")
parser$add_argument("--outDir", type = "character", required = TRUE)
parser$add_argument("--resultDir", type = "character", required = TRUE)
parser$add_argument("--outputFile", type = "character", default = "dashboard.html")
parser$add_argument("--title", type = "character", default = "FoodNetTrends Dashboard")
parser$add_argument("--debug", dest = "debug", action = "store_true")

args <- parser$parse_args()

if (args$debug) {
  cat("Clean Dashboard Arguments:\n")
  for (arg_name in names(args)) {
    cat(paste0("  ", arg_name, ": ", args[[arg_name]], "\n"))
  }
}

# ========================================================================
# Data Collection Functions
# ========================================================================

#' Scan directory for result files and images
scan_results <- function(result_dir) {
  if (!dir.exists(result_dir)) {
    cat("Creating output directory:", result_dir, "\n")
    dir.create(result_dir, recursive = TRUE)
  }
  
  files <- list.files(result_dir, full.names = TRUE, recursive = FALSE)
  
  result_data <- list(
    ir_files = files[grepl("_IRCatch\\.csv$", files)],
    png_files = files[grepl("\\.png$", files)],
    summary_files = files[grepl("_summary\\.txt$", files)],
    log_files = files[grepl("\\.log$", files)]
  )
  
  # Also check for symlinked files (Nextflow creates symlinks)
  symlinked_files <- list.files(result_dir, full.names = TRUE, recursive = FALSE)
  symlinked_files <- symlinked_files[file.info(symlinked_files)$islink]
  
  if (length(symlinked_files) > 0) {
    # Add symlinked files to appropriate categories
    result_data$png_files <- c(result_data$png_files, symlinked_files[grepl("\\.png$", symlinked_files)])
    result_data$summary_files <- c(result_data$summary_files, symlinked_files[grepl("_summary\\.txt$", symlinked_files)])
    result_data$ir_files <- c(result_data$ir_files, symlinked_files[grepl("_IRCatch\\.csv$", symlinked_files)])
  }
  
  cat("Found files:\n")
  cat("  IR files:", length(result_data$ir_files), "\n")
  cat("  PNG files:", length(result_data$png_files), "\n")
  cat("  Summary files:", length(result_data$summary_files), "\n")
  
  return(result_data)
}

#' Extract pathogen from filename
extract_pathogen <- function(filename) {
  basename <- basename(filename)
  pathogen <- gsub("_.*$", "", basename)
  return(toupper(pathogen))
}

#' Read IR data with basic processing
read_ir_data <- function(file_path) {
  tryCatch({
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    data$pathogen <- extract_pathogen(file_path)
    return(data)
  }, error = function(e) {
    cat("Warning: Could not read", file_path, "\n")
    return(NULL)
  })
}

#' Read summary data
read_summary_data <- function(file_path) {
  tryCatch({
    lines <- readLines(file_path)
    pathogen <- extract_pathogen(file_path)
    return(list(pathogen = pathogen, content = lines))
  }, error = function(e) {
    return(NULL)
  })
}

#' Convert image to base64 for embedding
embed_image <- function(image_path) {
  tryCatch({
    # Skip embedding if file is too large (>5MB) to prevent memory issues
    file_size <- file.info(image_path)$size
    if (file_size > 5 * 1024 * 1024) {
      cat("Skipping large image:", image_path, "(", round(file_size/1024/1024, 1), "MB)\n")
      return("")
    }
    
    # Read image file
    image_data <- readBin(image_path, "raw", file_size)
    
    # Convert to base64 - use base64enc if available, otherwise skip
    if (has_base64enc) {
      base64_data <- base64encode(image_data)
    } else {
      # If base64enc not available, skip image embedding
      cat("Warning: base64enc package not available, skipping image:", image_path, "\n")
      return("")
    }
    
    # Determine MIME type
    ext <- tolower(tools::file_ext(image_path))
    mime_type <- switch(ext,
                       "png" = "image/png",
                       "jpg" = "image/jpeg", 
                       "jpeg" = "image/jpeg",
                       "gif" = "image/gif",
                       "image/png")
    
    return(paste0("data:", mime_type, ";base64,", base64_data))
  }, error = function(e) {
    cat("Warning: Could not embed image:", image_path, ":", e$message, "\n")
    return("")
  })
}

# ========================================================================
# Dashboard HTML Generation
# ========================================================================

#' Generate clean CSS with muted color palette
generate_css <- function() {
  return('
<style>
/* Clean, Professional Styling with Muted Colors */
* {
  margin: 0;
  padding: 0;
  box-sizing: border-box;
}

body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  background: #fafafa;
  color: #2c3e50;
  line-height: 1.6;
  transition: background-color 0.3s ease, color 0.3s ease;
}

/* Dark mode styles */
body.dark-mode {
  background: #1a1a1a;
  color: #e0e0e0;
}

body.dark-mode .header {
  background: linear-gradient(135deg, #2c2c2c 0%, #1a1a1a 100%);
}

body.dark-mode .section {
  background: #2c2c2c;
  border-left-color: #555;
}

body.dark-mode .stat-card {
  background: #333;
  border-left-color: #4a90e2;
}

body.dark-mode .stat-number {
  color: #5dade2;
}

body.dark-mode .pathogen-card {
  background: #2c2c2c;
  border-color: #555;
}

body.dark-mode .pathogen-header {
  background: linear-gradient(45deg, #555, #444);
}

body.dark-mode .image-item {
  background: #2c2c2c;
}

body.dark-mode .data-table th {
  background: #555;
}

body.dark-mode .data-table tr:nth-child(even) {
  background: #333;
}

body.dark-mode .summary-text {
  background: #333;
  border-left-color: #555;
  color: #e0e0e0;
}

body.dark-mode .tab {
  background: #333;
  border-color: #555;
  color: #e0e0e0;
}

body.dark-mode .tab:hover {
  background: #444;
}

body.dark-mode .tab.active {
  background: #2c2c2c;
  border-color: #666;
  color: #fff;
}

body.dark-mode .tab-content {
  background: #2c2c2c;
  border-color: #555;
}

.container {
  max-width: 1200px;
  margin: 0 auto;
  padding: 20px;
}

.header {
  background: linear-gradient(135deg, #34495e 0%, #2c3e50 100%);
  color: white;
  padding: 30px 0;
  margin-bottom: 30px;
  border-radius: 8px;
  box-shadow: 0 2px 10px rgba(0,0,0,0.1);
  position: relative;
}

.header h1 {
  text-align: center;
  font-size: 2.5em;
  font-weight: 300;
  letter-spacing: -1px;
}

.header .subtitle {
  text-align: center;
  font-size: 1.1em;
  opacity: 0.9;
  margin-top: 10px;
}

/* Dark mode toggle button */
.dark-mode-toggle {
  position: absolute;
  top: 20px;
  right: 20px;
  background: rgba(255,255,255,0.2);
  border: 1px solid rgba(255,255,255,0.3);
  border-radius: 20px;
  padding: 8px 16px;
  color: white;
  cursor: pointer;
  font-size: 0.9em;
  transition: all 0.3s ease;
}

.dark-mode-toggle:hover {
  background: rgba(255,255,255,0.3);
  transform: translateY(-1px);
}

body.dark-mode .dark-mode-toggle {
  background: rgba(0,0,0,0.3);
  border-color: rgba(255,255,255,0.2);
}

body.dark-mode .dark-mode-toggle:hover {
  background: rgba(0,0,0,0.5);
}

.section {
  background: white;
  border-radius: 8px;
  padding: 25px;
  margin-bottom: 25px;
  box-shadow: 0 2px 5px rgba(0,0,0,0.05);
  border-left: 4px solid #95a5a6;
}

.section h2 {
  color: #34495e;
  font-size: 1.8em;
  margin-bottom: 20px;
  font-weight: 400;
}

.section h3 {
  color: #7f8c8d;
  font-size: 1.3em;
  margin: 20px 0 15px 0;
  font-weight: 500;
}

.stats-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
  gap: 20px;
  margin: 20px 0;
}

.stat-card {
  background: #ecf0f1;
  padding: 20px;
  border-radius: 6px;
  text-align: center;
  border-left: 4px solid #3498db;
}

.stat-number {
  font-size: 2.2em;
  font-weight: 600;
  color: #2980b9;
  margin-bottom: 5px;
}

.stat-label {
  color: #7f8c8d;
  font-size: 0.9em;
  text-transform: uppercase;
  letter-spacing: 1px;
}

.pathogen-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
  gap: 25px;
}

.pathogen-card {
  background: white;
  border: 1px solid #bdc3c7;
  border-radius: 8px;
  overflow: hidden;
  transition: transform 0.2s ease;
}

.pathogen-card:hover {
  transform: translateY(-2px);
  box-shadow: 0 4px 15px rgba(0,0,0,0.1);
}

.pathogen-header {
  background: linear-gradient(45deg, #95a5a6, #7f8c8d);
  color: white;
  padding: 15px 20px;
  font-weight: 600;
  font-size: 1.1em;
}

.pathogen-content {
  padding: 20px;
}

.image-gallery {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
  gap: 20px;
  margin: 20px 0;
}

.image-item {
  text-align: center;
  background: white;
  border-radius: 8px;
  padding: 15px;
  box-shadow: 0 2px 5px rgba(0,0,0,0.05);
}

.image-item img {
  max-width: 100%;
  height: auto;
  border-radius: 4px;
  border: 1px solid #ecf0f1;
  cursor: pointer;
  transition: transform 0.2s ease;
}

.image-item img:hover {
  transform: scale(1.02);
  box-shadow: 0 4px 12px rgba(0,0,0,0.15);
}

/* Modal for full-size image viewing */
.modal {
  display: none;
  position: fixed;
  z-index: 1000;
  left: 0;
  top: 0;
  width: 100%;
  height: 100%;
  background-color: rgba(0,0,0,0.9);
  cursor: pointer;
}

.modal-content {
  display: block;
  margin: auto;
  max-width: 90%;
  max-height: 90%;
  margin-top: 2.5%;
  background: white;
  padding: 20px;
  border-radius: 8px;
  box-shadow: 0 4px 20px rgba(0,0,0,0.5);
}

body.dark-mode .modal-content {
  background: #666;
}

.close {
  position: absolute;
  top: 20px;
  right: 35px;
  color: #fff;
  font-size: 40px;
  font-weight: bold;
  cursor: pointer;
  z-index: 1001;
}

.close:hover {
  color: #ccc;
}

.image-caption {
  margin-top: 10px;
  color: #7f8c8d;
  font-size: 0.9em;
}

.data-table {
  width: 100%;
  border-collapse: collapse;
  margin: 20px 0;
  font-size: 0.9em;
}

.data-table th {
  background: #95a5a6;
  color: white;
  padding: 12px 8px;
  text-align: left;
  font-weight: 500;
}

.data-table td {
  padding: 10px 8px;
  border-bottom: 1px solid #ecf0f1;
}

.data-table tr:nth-child(even) {
  background: #f8f9fa;
}

.summary-text {
  background: #f8f9fa;
  padding: 15px;
  border-radius: 4px;
  border-left: 3px solid #95a5a6;
  font-family: Monaco, "Lucida Console", monospace;
  font-size: 0.85em;
  line-height: 1.4;
  color: #2c3e50;
  overflow-x: auto;
}

.footer {
  text-align: center;
  color: #95a5a6;
  font-size: 0.9em;
  margin-top: 40px;
  padding: 20px;
  border-top: 1px solid #ecf0f1;
}

/* Tabbed data table styling */
.tabs {
  display: flex;
  border-bottom: 2px solid #ecf0f1;
  margin-bottom: 20px;
  flex-wrap: wrap;
}

.tab {
  background: #f8f9fa;
  border: 1px solid #ecf0f1;
  border-bottom: none;
  padding: 12px 20px;
  cursor: pointer;
  transition: all 0.2s ease;
  margin-right: 2px;
  font-weight: 500;
  border-radius: 4px 4px 0 0;
}

.tab:hover {
  background: #e9ecef;
}

.tab.active {
  background: white;
  border-color: #95a5a6;
  border-bottom: 2px solid white;
  margin-bottom: -2px;
  color: #2c3e50;
}

.tab-content {
  display: none;
  background: white;
  border-radius: 0 4px 4px 4px;
  padding: 20px;
  border: 1px solid #ecf0f1;
}

.tab-content.active {
  display: block;
}

.tab-content h4 {
  color: #34495e;
  margin-bottom: 15px;
  font-size: 1.2em;
}

@media (max-width: 768px) {
  .container { padding: 10px; }
  .header h1 { font-size: 2em; }
  .stats-grid { grid-template-columns: 1fr; }
  .pathogen-grid { grid-template-columns: 1fr; }
  .tabs { flex-direction: column; }
  .tab { margin-right: 0; margin-bottom: 2px; border-radius: 4px; }
}
</style>
')
}

#' Generate overview statistics
generate_overview <- function(result_data, ir_data_list) {
  n_pathogens <- length(unique(sapply(result_data$ir_files, extract_pathogen)))
  n_images <- length(result_data$png_files)
  
  # Calculate total observations
  total_obs <- 0
  if (length(ir_data_list) > 0) {
    total_obs <- sum(sapply(ir_data_list, function(x) if(!is.null(x)) nrow(x) else 0))
  }
  
  return(paste0('
<div class="section">
  <h2>Analysis Overview</h2>
  <div class="stats-grid">
    <div class="stat-card">
      <div class="stat-number">', n_pathogens, '</div>
      <div class="stat-label">Pathogens Analyzed</div>
    </div>
    <div class="stat-card">
      <div class="stat-number">', n_images, '</div>
      <div class="stat-label">Visualizations Generated</div>
    </div>
    <div class="stat-card">
      <div class="stat-number">', total_obs, '</div>
      <div class="stat-label">Total Observations</div>
    </div>
    <div class="stat-card">
      <div class="stat-number">', format(Sys.Date(), "%Y-%m-%d"), '</div>
      <div class="stat-label">Analysis Date</div>
    </div>
  </div>
</div>
'))
}

#' Generate image gallery section
generate_image_gallery <- function(png_files) {
  if (length(png_files) == 0) {
    return('<div class="section"><h2>Visualizations</h2><p>No visualization images found.</p></div>')
  }
  
  # Check if we can embed images
  if (!has_base64enc) {
    # If no base64enc, create file list instead of embedded images
    html <- '<div class="section"><h2>Visualizations</h2><p><strong>Note:</strong> Images cannot be embedded (base64enc package unavailable), but the following visualization files were generated:</p><ul>'
    for (png_file in png_files) {
      pathogen <- extract_pathogen(png_file)
      basename <- basename(png_file)
      html <- paste0(html, '<li><strong>', pathogen, ':</strong> ', basename, '</li>')
    }
    html <- paste0(html, '</ul></div>')
    return(html)
  }
  
  # Group images by pathogen and type
  image_groups <- list()
  for (png_file in png_files) {
    pathogen <- extract_pathogen(png_file)
    basename <- basename(png_file)
    
    if (grepl("_overall\\.png$", basename)) {
      type <- "Overall Trends"
    } else if (grepl("_state_trends\\.png$", basename)) {
      type <- "State Trends"
    } else if (grepl("_trend\\.png$", basename)) {
      type <- "Trend Analysis"
    } else {
      type <- "Analysis"
    }
    
    if (is.null(image_groups[[pathogen]])) {
      image_groups[[pathogen]] <- list()
    }
    image_groups[[pathogen]][[type]] <- png_file
  }
  
  html <- '<div class="section"><h2>Visualizations</h2>'
  
  for (pathogen in names(image_groups)) {
    html <- paste0(html, '<h3>', pathogen, '</h3><div class="image-gallery">')
    
    for (type in names(image_groups[[pathogen]])) {
      image_path <- image_groups[[pathogen]][[type]]
      base64_img <- embed_image(image_path)
      
      if (base64_img != "") {
        image_id <- paste0("img_", pathogen, "_", gsub("[^A-Za-z0-9]", "_", type))
        html <- paste0(html, '
<div class="image-item">
  <img src="', base64_img, '" alt="', pathogen, ' - ', type, '" onclick="openModal(\'', image_id, '\')">
  <div class="image-caption">', pathogen, ' - ', type, ' (click to expand)</div>
</div>')
      }
    }
    html <- paste0(html, '</div>')
  }
  
  html <- paste0(html, '</div>')
  return(html)
}

#' Generate tabbed pathogen data interface
generate_pathogen_summaries <- function(ir_data_list, summary_data_list) {
  if (length(ir_data_list) == 0) {
    return('<div class="section"><h2>Pathogen Data Analysis</h2><p>No analysis data found.</p></div>')
  }
  
  # Get unique pathogens from both IR data and summary data
  ir_pathogens <- unique(sapply(ir_data_list, function(x) if(!is.null(x)) x$pathogen[1] else NULL))
  ir_pathogens <- ir_pathogens[!is.null(ir_pathogens)]
  
  summary_pathogens <- unique(sapply(summary_data_list, function(x) if(!is.null(x)) x$pathogen else NULL))
  summary_pathogens <- summary_pathogens[!is.null(summary_pathogens)]
  
  # Combine and get all unique pathogens
  pathogens <- unique(c(ir_pathogens, summary_pathogens))
  pathogens <- sort(pathogens[!is.null(pathogens) & pathogens != ""])
  
  if (args$debug) {
    cat("Debug - IR pathogens found:", paste(ir_pathogens, collapse=", "), "\n")
    cat("Debug - Summary pathogens found:", paste(summary_pathogens, collapse=", "), "\n")
    cat("Debug - All pathogens for tabs:", paste(pathogens, collapse=", "), "\n")
  }
  
  html <- '<div class="section"><h2>Pathogen Data Analysis</h2>'
  
  # Generate tabs
  html <- paste0(html, '<div class="tabs">')
  for (i in seq_along(pathogens)) {
    pathogen <- pathogens[i]
    active_class <- if (i == 1) " active" else ""
    html <- paste0(html, '<div class="tab', active_class, '" onclick="showTab(\'', pathogen, '\')">', pathogen, '</div>')
  }
  html <- paste0(html, '</div>')
  
  # Generate tab content
  for (i in seq_along(pathogens)) {
    pathogen <- pathogens[i]
    active_class <- if (i == 1) " active" else ""
    
    # Find data for this pathogen
    pathogen_ir <- NULL
    pathogen_summary <- NULL
    
    for (data in ir_data_list) {
      if (!is.null(data) && !is.null(data$pathogen) && data$pathogen[1] == pathogen) {
        pathogen_ir <- data
        break
      }
    }
    
    for (summary in summary_data_list) {
      if (!is.null(summary) && summary$pathogen == pathogen) {
        pathogen_summary <- summary
        break
      }
    }
    
    html <- paste0(html, '<div id="', pathogen, '" class="tab-content', active_class, '">')
    html <- paste0(html, '<h4>', pathogen, ' Analysis Results</h4>')
    
    # Add statistics and data table
    if (!is.null(pathogen_ir)) {
      n_obs <- nrow(pathogen_ir)
      html <- paste0(html, '<p><strong>Total Observations:</strong> ', n_obs, '</p>')
      
      # Add comprehensive data table with full scrolling
      if (n_obs > 0) {
        display_data <- pathogen_ir
        # Remove pathogen column for display
        display_data$pathogen <- NULL
        
        # Create scrollable container for full data
        html <- paste0(html, '<div style="overflow: auto; max-height: 400px; margin: 15px 0; border: 1px solid #ecf0f1; border-radius: 4px;">')
        html <- paste0(html, '<table class="data-table" style="margin: 0;">')
        if (ncol(display_data) > 0) {
          # Sticky header for better scrolling experience
          html <- paste0(html, '<thead style="position: sticky; top: 0; background: #95a5a6; z-index: 10;">')
          html <- paste0(html, '<tr>')
          for (col_name in names(display_data)) {
            html <- paste0(html, '<th style="position: sticky; top: 0;">', col_name, '</th>')
          }
          html <- paste0(html, '</tr></thead><tbody>')
          
          # Include ALL rows - no limit
          for (i in 1:nrow(display_data)) {
            html <- paste0(html, '<tr>')
            for (col_name in names(display_data)) {
              value <- display_data[i, col_name]
              if (is.numeric(value)) {
                value <- round(value, 3)
              }
              html <- paste0(html, '<td>', value, '</td>')
            }
            html <- paste0(html, '</tr>')
          }
          html <- paste0(html, '</tbody>')
        }
        html <- paste0(html, '</table></div>')
        
        # Show total count
        html <- paste0(html, '<p><em>Showing all ', n_obs, ' observations (scroll to view all data)</em></p>')
      }
    }
    
    # Add summary text if available
    if (!is.null(pathogen_summary) && length(pathogen_summary$content) > 0) {
      html <- paste0(html, '<h5>Analysis Summary:</h5>')
      summary_text <- paste(head(pathogen_summary$content, 15), collapse = "\n")
      html <- paste0(html, '<div class="summary-text">', gsub("\n", "<br>", summary_text), '</div>')
    }
    
    html <- paste0(html, '</div>')
  }
  
  html <- paste0(html, '</div>')
  return(html)
}

# ========================================================================
# Main Dashboard Generation
# ========================================================================

cat("Scanning for result files...\n")
result_data <- scan_results(args$resultDir)

cat("Processing result data...\n")
# Read IR data
ir_data_list <- list()
for (ir_file in result_data$ir_files) {
  if (args$debug) {
    cat("Processing IR file:", ir_file, "\n")
  }
  data <- read_ir_data(ir_file)
  if (!is.null(data)) {
    if (args$debug) {
      cat("  - Successfully loaded with pathogen:", data$pathogen[1], "\n")
    }
    ir_data_list[[length(ir_data_list) + 1]] <- data
  } else {
    if (args$debug) {
      cat("  - Failed to load\n")
    }
  }
}

# Read summary data  
summary_data_list <- list()
for (summary_file in result_data$summary_files) {
  if (args$debug) {
    cat("Processing summary file:", summary_file, "\n")
  }
  data <- read_summary_data(summary_file)
  if (!is.null(data)) {
    if (args$debug) {
      cat("  - Successfully loaded with pathogen:", data$pathogen, "\n")
    }
    summary_data_list[[length(summary_data_list) + 1]] <- data
  } else {
    if (args$debug) {
      cat("  - Failed to load\n")
    }
  }
}

cat("Processed data for", length(ir_data_list), "pathogens\n")
cat("Found summaries for", length(summary_data_list), "pathogens\n")

# Generate dashboard HTML
cat("Creating dashboard HTML...\n")

dashboard_html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>', args$title, '</title>
    ', generate_css(), '
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="dark-mode-toggle" onclick="toggleDarkMode()">🌙 Dark Mode</div>
            <h1>FoodNetTrends Report</h1>
            <div class="subtitle">Analysis: ', gsub("^.*?([0-9]{8}_[0-9]{6}).*$", "\\1", args$title), '</div>
        </div>
        
        ', generate_overview(result_data, ir_data_list), '
        
        ', generate_image_gallery(result_data$png_files), '
        
        ', generate_pathogen_summaries(ir_data_list, summary_data_list), '
        
        <div class="footer">
            Generated on ', Sys.time(), ' by FoodNetTrends Pipeline v1.0
        </div>
    </div>
    
    <!-- Modal for image viewing -->
    <div id="imageModal" class="modal" onclick="closeModal()">
        <span class="close" onclick="closeModal()">&times;</span>
        <img class="modal-content" id="modalImage">
    </div>
    
    <script>
    // Image modal functionality
    function openModal(imgId) {
        var modal = document.getElementById("imageModal");
        var modalImg = document.getElementById("modalImage");
        var img = event.target;
        modal.style.display = "block";
        modalImg.src = img.src;
    }
    
    function closeModal() {
        document.getElementById("imageModal").style.display = "none";
    }
    
    // Tab switching functionality
    function showTab(pathogenName) {
        // Hide all tab contents
        var contents = document.getElementsByClassName("tab-content");
        for (var i = 0; i < contents.length; i++) {
            contents[i].classList.remove("active");
        }
        
        // Remove active class from all tabs
        var tabs = document.getElementsByClassName("tab");
        for (var i = 0; i < tabs.length; i++) {
            tabs[i].classList.remove("active");
        }
        
        // Show selected tab content
        document.getElementById(pathogenName).classList.add("active");
        
        // Add active class to clicked tab
        event.target.classList.add("active");
    }
    
    // Close modal on escape key
    document.addEventListener("keydown", function(event) {
        if (event.key === "Escape") {
            closeModal();
        }
    });
    
    // Dark mode toggle functionality
    function toggleDarkMode() {
        var body = document.body;
        var toggle = document.querySelector(".dark-mode-toggle");
        
        body.classList.toggle("dark-mode");
        
        if (body.classList.contains("dark-mode")) {
            toggle.innerHTML = "☀️ Light Mode";
            localStorage.setItem("darkMode", "enabled");
        } else {
            toggle.innerHTML = "🌙 Dark Mode";
            localStorage.setItem("darkMode", "disabled");
        }
    }
    
    // Check for saved dark mode preference on page load
    document.addEventListener("DOMContentLoaded", function() {
        var darkMode = localStorage.getItem("darkMode");
        var toggle = document.querySelector(".dark-mode-toggle");
        
        if (darkMode === "enabled") {
            document.body.classList.add("dark-mode");
            toggle.innerHTML = "☀️ Light Mode";
        }
    });
    </script>
</body>
</html>')

# Write dashboard file
output_path <- file.path(args$outDir, args$outputFile)
cat("Writing dashboard to", output_path, "\n")

tryCatch({
  writeLines(dashboard_html, output_path)
  cat("Dashboard generation complete\n")
}, error = function(e) {
  cat("ERROR: Could not write dashboard file:", e$message, "\n")
  quit(status = 1)
})

cat("Dashboard generation finished at", format(Sys.time()), "\n")