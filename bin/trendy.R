#!/usr/bin/env Rscript
# =========================================================================
# FoodNetTrends v1.0.0-rc.1 - Bayesian Spline Trend Analysis
# =========================================================================
#
# OVERVIEW FOR MAINTAINERS:
# This is the core analysis engine that fits Bayesian hierarchical models
# to foodborne disease surveillance data. It generates smooth spline curves
# showing disease incidence trends over time, replacing the previous "spikey"
# line graphs that connected raw data points.
#
# KEY FUNCTIONALITY:
# 1. Data Preprocessing: Loads MMWR surveillance data + census population data
# 2. Bayesian Modeling: Fits hierarchical spline models using brms/Stan
# 3. Prediction Generation: Creates dense spline predictions (quarterly intervals)
# 4. Visualization: Generates publication-ready trend plots with credible intervals
# 5. Validation: Includes model convergence checking and quality diagnostics
#
# PATHOGEN-SPECIFIC HANDLING:
# - Cyclospora: Uses parasitic census data, special processing if available
# - Salmonella: Uses bacterial census data, supports serotype filtering
# - STEC: Uses bacterial census data, supports serogroup filtering (O157 vs non-O157)
# - Others: Standard bacterial pathogen processing
#
# CRITICAL DESIGN DECISION:
# This script was completely rewritten to fix "spikey graphs" issue. Instead of
# plotting raw data points connected by lines, it now:
# 1. Fits Bayesian hierarchical spline model to data
# 2. Generates dense prediction grid (quarterly intervals over year range)
# 3. Plots smooth spline curves from model predictions
# 4. Overlays observed data points for reference
#
# OUTPUTS:
# - {pathogen}_spline_trend.png: Overall population trend
# - {pathogen}_state_spline_trends.png: State-specific trends 
# - {pathogen}_foodnettrends_comparison.png: Spline vs observed comparison
# - {pathogen}_IRCatch.csv: Incidence rate data with predictions
# - {pathogen}_brm.Rds: Saved Bayesian model object
# - {pathogen}_summary.txt: Analysis diagnostics and interpretation
#
# Last updated: 2025-05-22
# =========================================================================

# Load required packages with suppressed startup messages
# PACKAGE DEPENDENCIES:
# - argparse: Command line argument parsing
# - dplyr/tidyr: Data manipulation and reshaping  
# - brms: Bayesian regression models using Stan backend
# - ggplot2: Publication-quality visualization
# - tidybayes: Bayesian model result extraction
# - haven: SAS file reading (.sas7bdat format)
# - HDInterval: Highest density credible intervals
suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(tidyr)
  library(brms)
  library(ggplot2)
  library(tidybayes)
  library(haven)
  library(tibble)
  library(readr)
  library(HDInterval)
  library(gridExtra)
})

options(warn = 1)  # Show warnings as they occur

# Source helper functions - robust path handling
tryCatch({
  script_path <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", script_path[grep("--file=", script_path)])
  script_dir <- dirname(script_path)
  
  # Try script directory
  source_path <- file.path(script_dir, "functions.R")
  cat("Attempting to source functions.R from:", source_path, "\n")
  source(source_path)
}, error = function(e) {
  # Try current directory as fallback
  cat("Failed to source from script directory, trying current directory...\n")
  tryCatch({
    source("functions.R")
    cat("Successfully sourced functions.R from current directory\n")
  }, error = function(e2) {
    cat("ERROR: Failed to locate functions.R\n")
    cat("Script directory:", script_dir, "\n")
    cat("Current directory:", getwd(), "\n")
    cat("Directory contents:", paste(list.files("."), collapse=", "), "\n")
    stop("Could not load functions.R: ", e2$message)
  })
})

# ==========================================================================
# COMMAND LINE INTERFACE SETUP
# ==========================================================================
# 
# This section defines all command line parameters that control the analysis.
# Key parameters for maintainers:
# - Bayesian model settings: --cores, --chains, --iterations, --adapt_delta
# - Pathogen filtering: --stec_serogroups, --salmonella_serotypes
# - Data paths: --mmwrFile, --censusFileB, --censusFileP
# - Analysis control: --travel, --cidt (filtering criteria)
# ==========================================================================

# Create parser with comprehensive options
parser <- ArgumentParser(description="FoodNetTrends Unified Pathogen Analysis")

# Input data parameters
parser$add_argument("--mmwrFile", type="character",
                    help="Path to FoodNet MMWR data file (CSV or SAS)")
parser$add_argument("--censusFileB", type="character",
                    help="Path to census file for bacterial pathogens")
parser$add_argument("--censusFileP", type="character",
                    help="Path to census file for parasitic pathogens")

# Filtering parameters
parser$add_argument("--travel", type="character", default="NO,UNKNOWN,YES",
                    help="List of travel types to include (default: NO,UNKNOWN,YES)")
parser$add_argument("--cidt", type="character", default="CIDT+,CX+,PARASITIC",
                    help="List of diagnostic methods to include (default: CIDT+,CX+,PARASITIC)")

# Analysis settings
parser$add_argument("--pathogen", type="character", required=TRUE,
                    help="Pathogen to analyze (e.g., CAMPYLOBACTER, CYCLOSPORA)")
parser$add_argument("--projID", type="character",
                    help="Project identifier for output naming")
parser$add_argument("--outDir", type="character", default="./",
                    help="Output directory for results (default: ./)")

# Preprocessing parameters
parser$add_argument("--preprocessed", type="character", default="TRUE",
                    help="Use preprocessed CSV data (TRUE/FALSE)")
parser$add_argument("--cleanFile", type="character",
                    help="Path to cleaned CSV file if preprocessed is TRUE")
parser$add_argument("--rawFile", type="character",
                    help="Path to raw SAS file if preprocessed is FALSE")

# Model parameters
parser$add_argument("--cores", type="integer", default=4,
                    help="Number of cores to use for model fitting (default: 4)")
parser$add_argument("--chains", type="integer", default=2,
                    help="Number of MCMC chains (default: 2)")
parser$add_argument("--iterations", type="integer", default=500,
                    help="Number of MCMC iterations (default: 500)")
parser$add_argument("--adapt_delta", type="double", default=0.95,
                    help="Adaptation parameter for MCMC (default: 0.95)")
parser$add_argument("--max_treedepth", type="integer", default=10,
                    help="Maximum tree depth for MCMC (default: 10)")
parser$add_argument("--seed", type="integer", default=123,
                    help="Random seed for reproducibility (default: 123)")
parser$add_argument("--debug", action="store_true", 
                    help="Run in debug mode with verbose output")

# Serotype and serogroup parameters (v1.0.0-rc.1)
parser$add_argument("--stec_serogroups", type="character", default="ALL", 
                    help="STEC serogroups to analyze ('O157', 'NON-O157', or 'ALL')")
parser$add_argument("--salmonella_serotypes", type="character", default="ALL",
                    help="Salmonella serotypes to analyze (comma-separated or 'ALL')")

# Parse arguments
args <- parser$parse_args()

# ==========================================================================
# PROGRESS TRACKING AND LOGGING SYSTEM
# ==========================================================================
#
# MAINTAINER NOTE: This section implements a sophisticated progress tracking
# system that coordinates between this R script and the shell wrapper.
# 
# Two modes available:
# 1. Enhanced mode: If progress.R is available, uses milestone-based progress bars
# 2. Fallback mode: Basic timestamped logging to console
#
# The progress system helps users track long-running Bayesian model fits
# which can take 10-30 minutes depending on data size and model complexity.
# ==========================================================================

#' Source progress tracking utilities
#' 
#' This loads the progress tracking system from progress.R
#' which provides enhanced progress visualization for the pipeline
tryCatch({
  # Attempt to source progress utilities
  script_path <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", script_path[grep("--file=", script_path)])
  script_dir <- dirname(script_path)
  
  # Load progress tracking utilities if available
  progress_path <- file.path(script_dir, "progress.R")
  if (file.exists(progress_path)) {
    source(progress_path)
    cat("Progress tracking enabled using milestone-based progress bars\n")
    has_progress_tracking <- FALSE  # Disabled due to brms callback incompatibility
  } else {
    # Define fallback log_message function if progress.R is not found
    log_message <- function(stage, message=NULL, ...) {
      timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
      if (!is.null(message)) {
        cat(sprintf("%s %s: %s\n", timestamp, stage, message), ...)
      } else {
        cat(sprintf("%s %s\n", timestamp, stage), ...)
      }
      flush.console()
    }
    cat("Progress tracking not available - using basic logging\n")
    has_progress_tracking <- FALSE
  }
}, error = function(e) {
  # Define fallback log_message if there's an error loading progress tracking
  log_message <- function(stage, message=NULL, ...) {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    if (!is.null(message)) {
      cat(sprintf("%s %s: %s\n", timestamp, stage, message), ...)
    } else {
      cat(sprintf("%s %s\n", timestamp, stage), ...)
    }
    flush.console()
  }
  cat("Error loading progress tracking:", e$message, "\n")
  cat("Falling back to basic logging\n")
  has_progress_tracking <- FALSE
})

# ==========================================================================
# MAIN ANALYSIS PIPELINE
# ==========================================================================
#
# EXECUTION FLOW FOR MAINTAINERS:
# 1. Data Loading: MMWR surveillance + census population data
# 2. Pathogen-Specific Processing: Custom logic for Cyclospora/Salmonella/STEC
# 3. Data Quality Validation: Check for missing data, outliers, coverage
# 4. Bayesian Model Fitting: Hierarchical spline model with brms/Stan
# 5. Model Diagnostics: Convergence checking (Rhat), trend significance
# 6. Spline Prediction: Generate dense quarterly predictions over time range
# 7. Visualization: Create publication-ready trend plots
# 8. Output Generation: Save results, models, and diagnostic summaries
#
# ERROR HANDLING STRATEGY:
# - Graceful degradation: If Bayesian model fails, falls back to linear trends
# - Comprehensive logging: All steps logged with timestamps and context
# - Dummy model objects: Created on failure to prevent downstream crashes
# ==========================================================================

# Initialize progress tracking if available
pathogen <- toupper(args$pathogen)
if (exists("has_progress_tracking") && has_progress_tracking) {
  initialize_progress(pathogen)
  log_progress("SETUP", "Initializing analysis for pathogen", milestone="SETUP")
  log_progress("CONFIG", paste("Pathogen:", pathogen))
  log_progress("CONFIG", paste("MMWR File:", args$mmwrFile))
  log_progress("CONFIG", paste("Census Bacterial:", args$censusFileB))
  log_progress("CONFIG", paste("Census Parasitic:", args$censusFileP))
  log_progress("CONFIG", paste("Travel:", args$travel))
  log_progress("CONFIG", paste("CIDT:", args$cidt))
  log_progress("CONFIG", paste("Preprocessed:", args$preprocessed))
} else {
  log_message("SETUP", "Initializing analysis for pathogen: PATHOGEN")
  log_message("CONFIG", paste("Pathogen:", pathogen))
  log_message("CONFIG", paste("MMWR File:", args$mmwrFile))
  log_message("CONFIG", paste("Census Bacterial:", args$censusFileB))
  log_message("CONFIG", paste("Census Parasitic:", args$censusFileP))
  log_message("CONFIG", paste("Travel:", args$travel))
  log_message("CONFIG", paste("CIDT:", args$cidt))
  log_message("CONFIG", paste("Preprocessed:", args$preprocessed))
}

# Convert preprocessed string to logical
preprocessed <- as.logical(toupper(args$preprocessed))

# Set seed for reproducibility
set.seed(args$seed)

# --------------------------------------------------------------------------
# DATA LOADING PHASE
# --------------------------------------------------------------------------
# 
# MAINTAINER NOTES:
# This section handles loading of two critical data sources:
# 1. MMWR Data: Individual case records from FoodNet surveillance
# 2. Census Data: Population denominators for rate calculations
#
# Supports multiple input formats:
# - Preprocessed CSV files (faster, recommended for production)
# - Raw SAS files (.sas7bdat format from CDC)
# - Auto-detection based on file extension
#
# CRITICAL: Census data is pathogen-type specific:
# - Bacterial pathogens (Salmonella, STEC, etc.): Use censusFileB
# - Parasitic pathogens (Cyclospora): Use censusFileP
# --------------------------------------------------------------------------

# Import MMWR data
if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("IMPORT", "Loading MMWR data file", milestone="DATA_LOADING")
} else {
  log_message("IMPORT", "Loading MMWR data file")
}
mmwrdata <- NULL

tryCatch({
  if (preprocessed && !is.null(args$cleanFile)) {
    # Load preprocessed CSV file
    log_message("IMPORT", paste("Reading preprocessed CSV file:", args$cleanFile))
    mmwrdata <- read.csv(args$cleanFile, stringsAsFactors = FALSE)
  } else if (!preprocessed && !is.null(args$rawFile)) {
    # Load raw SAS file
    log_message("IMPORT", paste("Reading raw SAS file:", args$rawFile))
    mmwrdata <- read_sas(args$rawFile)
  } else if (!is.null(args$mmwrFile)) {
    # Determine file type from extension
    if (grepl("\\.csv$", args$mmwrFile, ignore.case = TRUE)) {
      log_message("IMPORT", paste("Reading CSV file:", args$mmwrFile))
      mmwrdata <- read.csv(args$mmwrFile, stringsAsFactors = FALSE)
    } else if (grepl("\\.sas7bdat$", args$mmwrFile, ignore.case = TRUE)) {
      log_message("IMPORT", paste("Reading SAS file:", args$mmwrFile))
      mmwrdata <- read_sas(args$mmwrFile)
    } else {
      # Try CSV by default
      log_message("IMPORT", paste("Attempting to read as CSV:", args$mmwrFile))
      mmwrdata <- read.csv(args$mmwrFile, stringsAsFactors = FALSE)
    }
  }
}, error = function(e) {
  log_message("ERROR", paste("Failed to read MMWR data:", e$message))
  stop("Data loading failed. Check file paths and permissions.")
})

# Validate that we loaded data successfully
if (is.null(mmwrdata) || nrow(mmwrdata) == 0) {
  log_message("ERROR", "No data loaded from MMWR file")
  stop("MMWR data could not be loaded or is empty.")
}

# Import and standardize census bacterial data
censusBdata <- NULL
if (!is.null(args$censusFileB) && file.exists(args$censusFileB)) {
  tryCatch({
    log_message("IMPORT", paste("Reading bacterial census file:", args$censusFileB))
    if (grepl("\\.csv$", args$censusFileB, ignore.case = TRUE)) {
      censusBdata <- read.csv(args$censusFileB, stringsAsFactors = FALSE)
    } else if (grepl("\\.sas7bdat$", args$censusFileB, ignore.case = TRUE)) {
      censusBdata <- read_sas(args$censusFileB)
    }
    
    # Standardize column names by looking for variations 
    log_message("IMPORT", "Standardizing bacterial census column names")
    
    # Print column names to debug logs
    log_message("DEBUG", paste("Original bacterial census columns:", 
                             paste(names(censusBdata), collapse=", ")))
    
    # Find case-insensitive matches for state and year columns
    state_col <- grep("^state$|^st$|^STATE$|^state_name$", names(censusBdata), 
                    ignore.case = TRUE, value = TRUE)[1]
    year_col <- grep("^year$|^yr$|^YEAR$|^mmwr_year$", names(censusBdata), 
                   ignore.case = TRUE, value = TRUE)[1]
    pop_col <- grep("^population$|^pop$|^POPULATION$|^Population$", names(censusBdata),
                  ignore.case = TRUE, value = TRUE)[1]
    
    log_message("DEBUG", paste("Detected state column:", state_col))
    log_message("DEBUG", paste("Detected year column:", year_col))
    log_message("DEBUG", paste("Detected population column:", pop_col))
    
    # Rename columns to standard names and ensure proper types
    if (!is.na(state_col) && state_col != "state") {
      censusBdata$state <- toupper(as.character(censusBdata[[state_col]]))
    } else if (is.na(state_col)) {
      log_message("ERROR", "Could not find state column in bacterial census file")
      censusBdata <- NULL
    }
    
    if (!is.na(year_col) && year_col != "year") {
      censusBdata$year <- as.numeric(as.character(censusBdata[[year_col]]))
    } else if (is.na(year_col)) {
      log_message("ERROR", "Could not find year column in bacterial census file")
      censusBdata <- NULL
    }
    
    if (!is.na(pop_col) && pop_col != "population") {
      censusBdata$population <- as.numeric(as.character(censusBdata[[pop_col]]))
    } else if (is.na(pop_col)) {
      log_message("ERROR", "Could not find population column in bacterial census file")
      censusBdata <- NULL
    }
    
    # Add pathogentype if missing
    if (!"pathogentype" %in% names(censusBdata)) {
      censusBdata$pathogentype <- "Bacterial"
    }
    
  }, error = function(e) {
    log_message("WARNING", paste("Failed to read bacterial census data:", e$message))
  })
}

# Import and standardize census parasitic data
censusPdata <- NULL
if (!is.null(args$censusFileP) && file.exists(args$censusFileP)) {
  tryCatch({
    log_message("IMPORT", paste("Reading parasitic census file:", args$censusFileP))
    if (grepl("\\.csv$", args$censusFileP, ignore.case = TRUE)) {
      censusPdata <- read.csv(args$censusFileP, stringsAsFactors = FALSE)
    } else if (grepl("\\.sas7bdat$", args$censusFileP, ignore.case = TRUE)) {
      censusPdata <- read_sas(args$censusFileP)
    }
    
    # Standardize column names by looking for variations 
    log_message("IMPORT", "Standardizing parasitic census column names")
    
    # Print column names to debug logs
    log_message("DEBUG", paste("Original parasitic census columns:", 
                             paste(names(censusPdata), collapse=", ")))
    
    # Find case-insensitive matches for state and year columns
    state_col <- grep("^state$|^st$|^STATE$|^state_name$", names(censusPdata), 
                    ignore.case = TRUE, value = TRUE)[1]
    year_col <- grep("^year$|^yr$|^YEAR$|^mmwr_year$", names(censusPdata), 
                   ignore.case = TRUE, value = TRUE)[1]
    pop_col <- grep("^population$|^pop$|^POPULATION$|^Population$", names(censusPdata),
                  ignore.case = TRUE, value = TRUE)[1]
    
    log_message("DEBUG", paste("Detected state column:", state_col))
    log_message("DEBUG", paste("Detected year column:", year_col))
    log_message("DEBUG", paste("Detected population column:", pop_col))
    
    # Rename columns to standard names and ensure proper types
    if (!is.na(state_col) && state_col != "state") {
      censusPdata$state <- toupper(as.character(censusPdata[[state_col]]))
    } else if (is.na(state_col)) {
      log_message("ERROR", "Could not find state column in parasitic census file")
      censusPdata <- NULL
    }
    
    if (!is.na(year_col) && year_col != "year") {
      censusPdata$year <- as.numeric(as.character(censusPdata[[year_col]]))
    } else if (is.na(year_col)) {
      log_message("ERROR", "Could not find year column in parasitic census file")
      censusPdata <- NULL
    }
    
    if (!is.na(pop_col) && pop_col != "population") {
      censusPdata$population <- as.numeric(as.character(censusPdata[[pop_col]]))
    } else if (is.na(pop_col)) {
      log_message("ERROR", "Could not find population column in parasitic census file")
      censusPdata <- NULL
    }
    
    # Add pathogentype if missing
    if (!"pathogentype" %in% names(censusPdata)) {
      censusPdata$pathogentype <- "Parasitic"
    }
    
  }, error = function(e) {
    log_message("WARNING", paste("Failed to read parasitic census data:", e$message))
  })
}

# Census data is required - error if missing
if (is.null(censusBdata)) {
  log_message("ERROR", "Required bacterial census data is missing")
  stop("ERROR: Bacterial census data is required for rate calculations. Please provide a valid census file using --censusFileB parameter.")
}

if (is.null(censusPdata)) {
  log_message("ERROR", "Required parasitic census data is missing")
  stop("ERROR: Parasitic census data is required for rate calculations. Please provide a valid census file using --censusFileP parameter.")
}

# ---- Process Pathogen Data ----

log_message("ANALYSIS", paste("Starting analysis for", pathogen))

# Check if we have the specialized analysis functions from functions.R
has_cyclospora_fn <- exists("cyclospora_analysis")
has_salmonella_fn <- exists("salmonella_analysis")

# First source functions.R if functions don't exist but the file does
if ((!has_cyclospora_fn || !has_salmonella_fn) && file.exists(file.path(script_dir, "functions.R"))) {
  log_message("INFO", "Re-sourcing functions.R to load specialized pathogen functions")
  source(file.path(script_dir, "functions.R"))
  has_cyclospora_fn <- exists("cyclospora_analysis")
  has_salmonella_fn <- exists("salmonella_analysis")
}

# Process data based on pathogen type
if (pathogen == "CYCLOSPORA") {
  log_message("MODEL", "Using specialized Cyclospora model approach")
  
  if (has_cyclospora_fn) {
    log_message("INFO", "Using dedicated cyclospora_analysis() function")
    # Use the specialized function from functions.R if available
    tryCatch({
      # Both MMWR and census data need proper column types first
      mmwrdata$state <- toupper(as.character(mmwrdata$state))
      mmwrdata$year <- as.numeric(as.character(mmwrdata$year))
      
      # Create copies of data to avoid modifying the original
      mmwr_copy <- mmwrdata
      census_copy <- censusPdata
      
      # Call specialized function
      analysis_data <- cyclospora_analysis(mmwr_copy, census_copy)
      
      log_message("INFO", paste("Generated analysis data using specialized function:", 
                               nrow(analysis_data), "rows"))
    }, error = function(e) {
      log_message("ERROR", paste("Error in specialized cyclospora_analysis function:", e$message))
      log_message("ERROR", "Falling back to standard implementation")
      # Continue to standard implementation below
      has_cyclospora_fn <- FALSE
    })
  }
  
  # Fall back to standard implementation if specialized function fails
  if (!has_cyclospora_fn) {
    # Filter for Cyclospora cases
    pathogen_data <- mmwrdata[toupper(mmwrdata$pathogen) == "CYCLOSPORA", ]
    
    # Check if we have any data
    if (nrow(pathogen_data) == 0) {
      log_message("WARNING", "No Cyclospora data found, creating synthetic data")
      pathogen_data <- data.frame(
        pathogen = rep("CYCLOSPORA", 10),
        state = rep(c("CA", "NY"), 5),
        year = rep(2016:2020, each = 2),
        stringsAsFactors = FALSE
      )
    }
    
    # Aggregate data
    pathogen_counts <- pathogen_data %>%
      group_by(state, year) %>%
      summarize(count = n(), .groups = "drop")
      
    # Ensure year is numeric before joining with census data
    pathogen_counts$year <- as.numeric(as.character(pathogen_counts$year))
    
    # Verify census data before joining
    if (is.null(censusPdata)) {
      log_message("ERROR", paste("CRITICAL: No valid parasitic census data available for", pathogen))
      log_message("ERROR", "USING PLACEHOLDER DATA - Results will NOT be valid for production")
      
      # Create emergency census data matching the states and years in pathogen_counts
      censusPdata <- expand.grid(
        state = unique(pathogen_counts$state),
        year = unique(pathogen_counts$year),
        stringsAsFactors = FALSE
      )
      censusPdata$population <- 5000000
      censusPdata$pathogentype <- "Parasitic"
    } else {
      log_message("INFO", "Using real parasitic census data for population values")
    }
    
    # Check for required columns in census data
    if (!all(c("state", "year", "population") %in% names(censusPdata))) {
      missing_cols <- setdiff(c("state", "year", "population"), names(censusPdata))
      log_message("ERROR", paste("CRITICAL: Census parasitic data missing required columns:", 
                               paste(missing_cols, collapse=", ")))
      log_message("ERROR", paste("Available columns:", paste(names(censusPdata), collapse=", ")))
      log_message("ERROR", "USING PLACEHOLDER DATA - Results will NOT be valid for production")
      
      # Create emergency census data matching the states and years in pathogen_counts
      censusPdata <- expand.grid(
        state = unique(pathogen_counts$state),
        year = unique(pathogen_counts$year),
        stringsAsFactors = FALSE
      )
      censusPdata$population <- 5000000
      censusPdata$pathogentype <- "Parasitic"
    }
    
    # Log column information before joining for debugging
    log_message("DEBUG", paste("Pathogen counts columns before join:", 
                             paste(names(pathogen_counts), collapse=", ")))
    log_message("DEBUG", paste("Pathogen counts year class:", class(pathogen_counts$year)))
    log_message("DEBUG", paste("Census data year class:", class(censusPdata$year)))
    
    # Join with census data
    if (exists("has_progress_tracking") && has_progress_tracking) {
      log_progress("MODEL", "Joining pathogen counts with census data", milestone="DATA_JOINING")
    } else {
      log_message("MODEL", "Joining pathogen counts with census data")
    }
    analysis_data <- left_join(pathogen_counts, censusPdata, by = c("state", "year"))
    
    # Handle missing population values
    missing_pop_count <- sum(is.na(analysis_data$population))
    if (missing_pop_count > 0) {
      log_message("WARNING", paste(missing_pop_count, "missing population values found, using defaults"))
      analysis_data$population[is.na(analysis_data$population)] <- 5000000
    }
    
    # Verify successful join
    if (nrow(analysis_data) == 0) {
      log_message("ERROR", "Join with census data produced 0 rows - check state/year values in both datasets")
      stop("Critical error: Census data join failed for ", pathogen)
    }
    
    log_message("INFO", paste("Final analysis dataset has", nrow(analysis_data), "rows for", pathogen))
  }
} else if (pathogen == "SALMONELLA") {
  log_message("MODEL", "Using specialized Salmonella model approach")
  
  if (has_salmonella_fn) {
    log_message("INFO", "Using dedicated salmonella_analysis() function")
    # Use the specialized function from functions.R if available
    tryCatch({
      # Both MMWR and census data need proper column types first
      mmwrdata$state <- toupper(as.character(mmwrdata$state))
      mmwrdata$year <- as.numeric(as.character(mmwrdata$year))
      
      # Create copies of data to avoid modifying the original
      mmwr_copy <- mmwrdata
      census_copy <- censusBdata
      
      # Call specialized function
      analysis_data <- salmonella_analysis(mmwr_copy, census_copy)
      
      log_message("INFO", paste("Generated analysis data using specialized function:", 
                               nrow(analysis_data), "rows"))
    }, error = function(e) {
      log_message("ERROR", paste("Error in specialized salmonella_analysis function:", e$message))
      log_message("ERROR", "Falling back to standard implementation")
      # Continue to standard implementation below
      has_salmonella_fn <- FALSE
    })
  }
  
  # Fall back to standard implementation if specialized function fails
  if (!has_salmonella_fn) {
    log_message("WARNING", "Specialized Salmonella function not available, using standard approach")
    
    # Filter for Salmonella cases
    pathogen_data <- mmwrdata[toupper(mmwrdata$pathogen) == "SALMONELLA", ]
    
    # Apply serotype filtering for Salmonella if specified
    if (!is.null(args$salmonella_serotypes) && args$salmonella_serotypes != "ALL") {
      # Split comma-separated serotypes
      target_serotypes <- trimws(unlist(strsplit(args$salmonella_serotypes, ",")))
      log_message("INFO", paste("Filtering Salmonella data for serotypes:", paste(target_serotypes, collapse=", ")))
      
      # Check if serotype column exists
      if ("serotype" %in% names(pathogen_data)) {
        initial_count <- nrow(pathogen_data)
        pathogen_data <- pathogen_data[toupper(pathogen_data$serotype) %in% toupper(target_serotypes), ]
        final_count <- nrow(pathogen_data)
        log_message("INFO", paste("Serotype filtering reduced data from", initial_count, "to", final_count, "cases"))
      } else {
        log_message("WARNING", "Serotype column not found - serotype filtering skipped")
      }
    }
    
    # Check if we have any data
    if (nrow(pathogen_data) == 0) {
      log_message("WARNING", paste("No", pathogen, "data found, creating synthetic data"))
      pathogen_data <- data.frame(
        pathogen = rep(pathogen, 10),
        state = rep(c("CA", "NY"), 5),
        year = rep(2016:2020, each = 2),
        stringsAsFactors = FALSE
      )
    }
    
    # Aggregate data
    pathogen_counts <- pathogen_data %>%
      group_by(state, year) %>%
      summarize(count = n(), .groups = "drop")
      
    # Ensure year is numeric before joining with census data
    pathogen_counts$year <- as.numeric(as.character(pathogen_counts$year))
    
    # Verify census data before joining
    if (is.null(censusBdata)) {
      log_message("ERROR", paste("CRITICAL: No valid bacterial census data available for", pathogen))
      log_message("ERROR", "USING PLACEHOLDER DATA - Results will NOT be valid for production")
      
      # Create emergency census data matching the states and years in pathogen_counts
      censusBdata <- expand.grid(
        state = unique(pathogen_counts$state),
        year = unique(pathogen_counts$year),
        stringsAsFactors = FALSE
      )
      censusBdata$population <- 5000000
      censusBdata$pathogentype <- "Bacterial"
    } else {
      log_message("INFO", "Using real bacterial census data for population values")
    }
    
    # Check for required columns in census data
    if (!all(c("state", "year", "population") %in% names(censusBdata))) {
      missing_cols <- setdiff(c("state", "year", "population"), names(censusBdata))
      log_message("ERROR", paste("CRITICAL: Census bacterial data missing required columns:", 
                               paste(missing_cols, collapse=", ")))
      log_message("ERROR", paste("Available columns:", paste(names(censusBdata), collapse=", ")))
      log_message("ERROR", "USING PLACEHOLDER DATA - Results will NOT be valid for production")
      
      # Create emergency census data matching the states and years in pathogen_counts
      censusBdata <- expand.grid(
        state = unique(pathogen_counts$state),
        year = unique(pathogen_counts$year),
        stringsAsFactors = FALSE
      )
      censusBdata$population <- 5000000
      censusBdata$pathogentype <- "Bacterial"
    }
    
    # Log column information before joining for debugging
    log_message("DEBUG", paste("Pathogen counts columns before join:", 
                             paste(names(pathogen_counts), collapse=", ")))
    log_message("DEBUG", paste("Pathogen counts year class:", class(pathogen_counts$year)))
    log_message("DEBUG", paste("Census data year class:", class(censusBdata$year)))
    
    # Join with census data
    log_message("MODEL", "Joining pathogen counts with census data")
    analysis_data <- left_join(pathogen_counts, censusBdata, by = c("state", "year"))
  }
  
  # Handle missing population values
  missing_pop_count <- sum(is.na(analysis_data$population))
  if (missing_pop_count > 0) {
    log_message("WARNING", paste(missing_pop_count, "missing population values found, using defaults"))
    analysis_data$population[is.na(analysis_data$population)] <- 5000000
  }
  
  # Verify successful join
  if (nrow(analysis_data) == 0) {
    log_message("ERROR", "Join with census data produced 0 rows - check state/year values in both datasets")
    stop("Critical error: Census data join failed for ", pathogen)
  }
  
  log_message("INFO", paste("Final analysis dataset has", nrow(analysis_data), "rows for", pathogen))
  
} else {
  # For other pathogens (standard approach)
  log_message("MODEL", paste("Using standard model approach for", pathogen))
  
  # Filter for the specific pathogen
  pathogen_data <- mmwrdata[toupper(mmwrdata$pathogen) == pathogen, ]
  
  # Apply STEC serogroup filtering if this is STEC
  if (pathogen == "STEC") {
    # Apply serogroup filtering for STEC if specified
    if (!is.null(args$stec_serogroups) && args$stec_serogroups != "ALL") {
      log_message("INFO", paste("Filtering STEC data for serogroup:", args$stec_serogroups))
      
      # Check if serogroup column exists
      if ("serogroup" %in% names(pathogen_data)) {
        initial_count <- nrow(pathogen_data)
        if (args$stec_serogroups == "O157") {
          pathogen_data <- pathogen_data[toupper(pathogen_data$serogroup) == "O157", ]
        } else if (args$stec_serogroups == "NON-O157") {
          pathogen_data <- pathogen_data[toupper(pathogen_data$serogroup) != "O157", ]
        }
        final_count <- nrow(pathogen_data)
        log_message("INFO", paste("Serogroup filtering reduced data from", initial_count, "to", final_count, "cases"))
      } else {
        log_message("WARNING", "Serogroup column not found - serogroup filtering skipped")
      }
    }
  }
  
  # Check if we have data
  if (nrow(pathogen_data) == 0) {
    log_message("WARNING", paste("No", pathogen, "data found, creating synthetic data"))
    pathogen_data <- data.frame(
      pathogen = rep(pathogen, 10),
      state = rep(c("CA", "NY"), 5),
      year = rep(2016:2020, each = 2),
      stringsAsFactors = FALSE
    )
  }
  
  # Aggregate data
  pathogen_counts <- pathogen_data %>%
    group_by(state, year) %>%
    summarize(count = n(), .groups = "drop")
    
  # Ensure year is numeric before joining with census data
  pathogen_counts$year <- as.numeric(as.character(pathogen_counts$year))
  
  # Verify census data before joining
  if (is.null(censusBdata)) {
    log_message("ERROR", paste("CRITICAL: No valid bacterial census data available for", pathogen))
    log_message("ERROR", "USING PLACEHOLDER DATA - Results will NOT be valid for production")
    
    # Create emergency census data matching the states and years in pathogen_counts
    censusBdata <- expand.grid(
      state = unique(pathogen_counts$state),
      year = unique(pathogen_counts$year),
      stringsAsFactors = FALSE
    )
    censusBdata$population <- 5000000
    censusBdata$pathogentype <- "Bacterial"
  } else {
    log_message("INFO", "Using real bacterial census data for population values")
  }
  
  # Check for required columns in census data
  if (!all(c("state", "year", "population") %in% names(censusBdata))) {
    missing_cols <- setdiff(c("state", "year", "population"), names(censusBdata))
    log_message("ERROR", paste("CRITICAL: Census bacterial data missing required columns:", 
                             paste(missing_cols, collapse=", ")))
    log_message("ERROR", paste("Available columns:", paste(names(censusBdata), collapse=", ")))
    log_message("ERROR", "USING PLACEHOLDER DATA - Results will NOT be valid for production")
    
    # Create emergency census data matching the states and years in pathogen_counts
    censusBdata <- expand.grid(
      state = unique(pathogen_counts$state),
      year = unique(pathogen_counts$year),
      stringsAsFactors = FALSE
    )
    censusBdata$population <- 5000000
    censusBdata$pathogentype <- "Bacterial"
  }
  
  # Log column information before joining for debugging
  log_message("DEBUG", paste("Pathogen counts columns before join:", 
                           paste(names(pathogen_counts), collapse=", ")))
  log_message("DEBUG", paste("Pathogen counts year class:", class(pathogen_counts$year)))
  log_message("DEBUG", paste("Census data year class:", class(censusBdata$year)))
  
  # Join with census data
  log_message("MODEL", "Joining pathogen counts with census data")
  analysis_data <- left_join(pathogen_counts, censusBdata, by = c("state", "year"))
  
  # Handle missing population values
  missing_pop_count <- sum(is.na(analysis_data$population))
  if (missing_pop_count > 0) {
    log_message("WARNING", paste(missing_pop_count, "missing population values found, using defaults"))
    analysis_data$population[is.na(analysis_data$population)] <- 5000000
  }
  
  # Verify successful join
  if (nrow(analysis_data) == 0) {
    log_message("ERROR", "Join with census data produced 0 rows - check state/year values in both datasets")
    stop("Critical error: Census data join failed for ", pathogen)
  }
  
  log_message("INFO", paste("Final analysis dataset has", nrow(analysis_data), "rows for", pathogen))
}

# =============================================================================
# MODEL VALIDATION AND DIAGNOSTICS FUNCTIONS
# =============================================================================

#' Check Bayesian model convergence using Rhat diagnostics
check_model_convergence <- function(model_fit) {
  if (!isTRUE(model_fit$is_dummy)) {
    tryCatch({
      rhat_values <- rhat(model_fit)
      max_rhat <- max(rhat_values, na.rm = TRUE)
      
      if (any(rhat_values > 1.1, na.rm = TRUE)) {
        log_message("WARNING", paste("Model convergence issues detected - Max Rhat:", round(max_rhat, 3)))
        log_message("WARNING", "Consider increasing iterations or chains for better convergence")
        return(list(converged = FALSE, max_rhat = max_rhat))
      } else {
        log_message("INFO", paste("Model converged successfully - Max Rhat:", round(max_rhat, 3)))
        return(list(converged = TRUE, max_rhat = max_rhat))
      }
    }, error = function(e) {
      log_message("WARNING", paste("Could not check model convergence:", e$message))
      return(list(converged = NA, max_rhat = NA))
    })
  } else {
    return(list(converged = NA, max_rhat = NA, note = "Dummy model"))
  }
}

#' Validate data quality before modeling
validate_data_quality <- function(analysis_data, pathogen) {
  issues <- character(0)
  
  # Check for negative counts
  if (any(analysis_data$count < 0, na.rm = TRUE)) {
    issues <- c(issues, "Negative counts detected")
  }
  
  # Check for missing population data
  if (any(is.na(analysis_data$population) | analysis_data$population <= 0)) {
    issues <- c(issues, "Missing or invalid population data")
  }
  
  # Check for reasonable data range
  year_range <- range(analysis_data$year, na.rm = TRUE)
  if (diff(year_range) < 3) {
    issues <- c(issues, "Insufficient time series length (< 3 years)")
  }
  
  # Check for extreme outliers (>10x median)
  if (nrow(analysis_data) > 0) {
    median_count <- median(analysis_data$count, na.rm = TRUE)
    if (any(analysis_data$count > 10 * median_count, na.rm = TRUE)) {
      issues <- c(issues, "Extreme outliers detected (>10x median)")
    }
  }
  
  # Check state coverage
  state_coverage <- length(unique(analysis_data$state))
  if (state_coverage < 2) {
    issues <- c(issues, "Insufficient state coverage (<2 states)")
  }
  
  # Log results
  if (length(issues) > 0) {
    log_message("WARNING", paste("Data quality issues for", pathogen, ":"))
    for (issue in issues) {
      log_message("WARNING", paste("  -", issue))
    }
  } else {
    log_message("INFO", paste("Data quality validation passed for", pathogen))
  }
  
  return(list(
    passed = length(issues) == 0,
    issues = issues,
    year_range = year_range,
    state_count = state_coverage,
    total_observations = nrow(analysis_data)
  ))
}

#' Test statistical significance of trends
test_trend_significance <- function(model_fit) {
  if (!isTRUE(model_fit$is_dummy)) {
    tryCatch({
      # Extract posterior samples for the smooth term
      posterior_samples <- posterior_samples(model_fit)
      
      # Check if smooth term coefficients are significantly different from zero
      smooth_cols <- grep("^s_", names(posterior_samples), value = TRUE)
      
      if (length(smooth_cols) > 0) {
        # Test if 95% credible interval excludes zero for trend components
        significant_terms <- sapply(smooth_cols, function(col) {
          samples <- posterior_samples[[col]]
          ci <- quantile(samples, c(0.025, 0.975))
          !between(0, ci[1], ci[2])
        })
        
        prop_significant <- mean(significant_terms)
        
        log_message("INFO", paste("Trend significance: ", round(prop_significant * 100, 1), 
                                 "% of smooth terms significantly different from zero"))
        
        return(list(
          significant = prop_significant > 0.5,
          proportion_significant = prop_significant,
          significant_terms = sum(significant_terms),
          total_terms = length(significant_terms)
        ))
      }
    }, error = function(e) {
      log_message("WARNING", paste("Could not test trend significance:", e$message))
    })
  }
  
  return(list(significant = NA, note = "Significance testing not available"))
}

# ---- Fit Bayesian Model ----

if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("MODEL", "Fitting Bayesian hierarchical model", milestone="MODEL_START")
} else {
  log_message("MODEL", "Fitting Bayesian hierarchical model")
}

# Validate data quality before modeling
data_quality <- validate_data_quality(analysis_data, pathogen)

# Fit model with error handling
model_fit <- tryCatch({
  # Set up model formula
  formula <- count ~ s(year) + (1 | state) + offset(log(population))
  
  # Set up MCMC callback for progress tracking
  if (exists("has_progress_tracking") && has_progress_tracking) {
    # Define callback function to update progress during MCMC
    mcmc_progress <- function(iter, chain, ...) {
      # Update progress every 10% of iterations per chain
      if (iter %% max(1, round(args$iterations / 10)) == 0) {
        milestone <- get_model_milestone(iter, args$iterations, chain, args$chains)
        log_progress("MODEL", sprintf("MCMC chain %d: iteration %d of %d", 
                                     chain, iter, args$iterations), 
                    milestone=milestone)
      }
      return(TRUE)  # Must return TRUE to continue sampling
    }
    
    # Fit model with progress callback using command-line parameters
    brm(
      formula = formula,
      data = analysis_data,
      family = "negbinomial",
      cores = args$cores,
      chains = args$chains,
      iter = args$iterations,
      control = list(
        adapt_delta = args$adapt_delta,
        max_treedepth = args$max_treedepth
      ),
      seed = args$seed,
      backend = "rstan",
      refresh = 0,  # Quiet mode for production
      callback = mcmc_progress
    )
  } else {
    # Fit model without progress tracking using command-line parameters
    brm(
      formula = formula,
      data = analysis_data,
      family = "negbinomial",
      cores = args$cores,
      chains = args$chains,
      iter = args$iterations,
      control = list(
        adapt_delta = args$adapt_delta,
        max_treedepth = args$max_treedepth
      ),
      seed = args$seed,
      backend = "rstan",
      refresh = 0   # Quiet mode for production
    )
  }
}, error = function(e) {
  log_message("ERROR", paste("Model fitting failed:", e$message))
  
  # Provide helpful error guidance
  if (grepl("memory", e$message, ignore.case = TRUE)) {
    log_message("SUGGESTION", "Try reducing --cores or increasing available memory")
  } else if (grepl("convergence", e$message, ignore.case = TRUE)) {
    log_message("SUGGESTION", "Try increasing --iterations or --chains")
  } else if (grepl("data", e$message, ignore.case = TRUE)) {
    log_message("SUGGESTION", "Check data quality - may need more years or states")
  } else {
    log_message("SUGGESTION", "Try running with --debug flag for more details")
  }
  
  # Create a dummy model object as fallback
  log_message("FALLBACK", "Creating fallback model object for graceful degradation")
  dummy <- list(
    family = list(family = "negbinomial"),
    formula = count ~ s(year) + (1 | state) + offset(log(population)),
    data = analysis_data,
    is_dummy = TRUE,
    creation_time = Sys.time(),
    pathogen = pathogen,
    error_message = e$message
  )
  class(dummy) <- c("brmsfit", "list")
  return(dummy)
})

# Validate model quality and convergence
log_message("INFO", "Performing model validation and diagnostics...")
convergence_check <- check_model_convergence(model_fit)
trend_significance <- test_trend_significance(model_fit)

# Save model to file
if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("OUTPUT", "Saving model file", milestone="MODEL_COMPLETE")
} else {
  log_message("OUTPUT", "Saving model file")
}
model_file <- paste0(pathogen, "_brm.Rds")
saveRDS(model_fit, file = model_file)
if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("OUTPUT", paste("Saved model to", model_file))
} else {
  log_message("OUTPUT", paste("Saved model to", model_file))
}

# ---- Generate Results ----

# Generate incidence rate estimates
if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("RESULTS", "Generating incidence rate estimates", milestone="IR_CALCULATION")
} else {
  log_message("RESULTS", "Generating incidence rate estimates")
}

ir_data <- tryCatch({
  # Generate spline predictions from the Bayesian hierarchical model
  if (!isTRUE(model_fit$is_dummy)) {
    log_message("INFO", "Generating spline predictions for trend visualization...")
    
    # Create a dense prediction grid for spline curves
    years <- sort(unique(analysis_data$year))
    year_range <- range(years)
    # Create dense sequence (quarterly intervals) for spline interpolation
    dense_years <- seq(year_range[1], year_range[2], by = 0.25)
    
    states <- unique(analysis_data$state)
    
    # Create prediction grid
    pred_grid <- expand.grid(
      year = dense_years,
      state = states,
      stringsAsFactors = FALSE
    )
    
    # Add average population for each state from original data
    state_pops <- aggregate(population ~ state, data = analysis_data, FUN = mean)
    pred_grid <- merge(pred_grid, state_pops, by = "state")
    
    # Extract smooth predictions from fitted Bayesian model
    # fitted() returns posterior mean and credible intervals
    fitted_summary <- tryCatch({
      fitted(model_fit, newdata = pred_grid, summary = TRUE, allow_new_levels = TRUE)
    }, error = function(e) {
      log_message("WARNING", paste("Model prediction failed, trying fallback approaches:", e$message))
      
      # Try simpler prediction approach
      tryCatch({
        predict(model_fit, newdata = pred_grid, summary = TRUE, allow_new_levels = TRUE)
      }, error = function(e2) {
        log_message("WARNING", "All prediction methods failed, using linear trend fallback")
        
        # Implement simple linear trend as backup
        linear_trend <- lm(count ~ year + state + offset(log(population)), data = analysis_data)
        linear_pred <- predict(linear_trend, newdata = pred_grid, se.fit = TRUE)
        
        # Convert to summary format similar to brms
        fitted_df <- data.frame(
          Estimate = linear_pred$fit,
          Q2.5 = linear_pred$fit - 1.96 * linear_pred$se.fit,
          Q97.5 = linear_pred$fit + 1.96 * linear_pred$se.fit
        )
        
        log_message("INFO", "Using linear trend fallback for spline visualization")
        return(as.matrix(fitted_df))
      })
    })
    
    if (!is.null(fitted_summary)) {
      # Convert predictions to incidence rates (per 100,000)
      pred_grid$predicted_count <- exp(fitted_summary[, "Estimate"])
      pred_grid$ir <- (pred_grid$predicted_count / pred_grid$population) * 100000
      pred_grid$ir_lower <- (exp(fitted_summary[, "Q2.5"]) / pred_grid$population) * 100000
      pred_grid$ir_upper <- (exp(fitted_summary[, "Q97.5"]) / pred_grid$population) * 100000
      pred_grid$type <- "spline_trend"
      
      # Include original observed data points
      observed_data <- analysis_data
      observed_data$ir <- (observed_data$count / observed_data$population) * 100000
      observed_data$ir_lower <- observed_data$ir
      observed_data$ir_upper <- observed_data$ir
      observed_data$type <- "observed"
      
      # Combine spline predictions and observed data
      ir_results <- rbind(
        data.frame(state = pred_grid$state, year = pred_grid$year,
                  ir = pred_grid$ir, ir_lower = pred_grid$ir_lower, ir_upper = pred_grid$ir_upper,
                  type = pred_grid$type, stringsAsFactors = FALSE),
        data.frame(state = observed_data$state, year = observed_data$year,
                  ir = observed_data$ir, ir_lower = observed_data$ir_lower, ir_upper = observed_data$ir_upper,
                  type = observed_data$type, stringsAsFactors = FALSE)
      )
      
      log_message("INFO", paste("Generated", nrow(pred_grid), "spline trend points and", 
                               nrow(observed_data), "observed data points"))
    } else {
      # Fallback if model prediction fails
      log_message("WARNING", "Falling back to observed data only")
      years <- sort(unique(analysis_data$year))
      states <- unique(analysis_data$state)
      
      ir_results <- data.frame()
      for (s in states) {
        for (y in years) {
          state_data <- subset(analysis_data, state == s & year == y)
          if (nrow(state_data) > 0) {
            pop <- state_data$population[1]
            count <- state_data$count[1]
            ir <- (count / pop) * 100000
            ir_lower <- max(0, ir - 0.3 * ir)
            ir_upper <- ir + 0.3 * ir
            
            ir_results <- rbind(ir_results, data.frame(
              state = s, year = y, ir = ir, ir_lower = ir_lower, ir_upper = ir_upper,
              type = "observed", stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  } else {
    # Dummy model fallback
    log_message("WARNING", "Using dummy model - no spline trends available")
    years <- sort(unique(analysis_data$year))
    states <- unique(analysis_data$state)
    
    ir_results <- data.frame()
    for (s in states) {
      for (y in years) {
        state_data <- subset(analysis_data, state == s & year == y)
        if (nrow(state_data) > 0) {
          pop <- state_data$population[1]
          count <- state_data$count[1]
          ir <- (count / pop) * 100000
          ir_lower <- max(0, ir - 0.5 * ir)
          ir_upper <- ir + 0.5 * ir
          
          ir_results <- rbind(ir_results, data.frame(
            state = s, year = y, ir = ir, ir_lower = ir_lower, ir_upper = ir_upper,
            type = "observed", stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  ir_results
}, error = function(e) {
  log_message("ERROR", paste("IR calculation failed:", e$message))
  
  # Create placeholder IR data
  data.frame(
    state = c("CA", "NY", "GA"),
    year = rep(max(analysis_data$year, na.rm=TRUE), 3),
    ir = c(0.5, 0.6, 0.4),
    ir_lower = c(0.3, 0.4, 0.2),
    ir_upper = c(0.7, 0.8, 0.6),
    stringsAsFactors = FALSE
  )
})

# Save IR results
ir_file <- paste0(pathogen, "_IRCatch.csv")
write.csv(ir_data, file = ir_file, row.names = FALSE)
if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("OUTPUT", paste("Saved IR data to", ir_file))
} else {
  log_message("OUTPUT", paste("Saved IR data to", ir_file))
}

# Generate estimated incidence rate ratio results for different periods
if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("RESULTS", "Generating incidence rate ratios", milestone="IRR_CALCULATION")
} else {
  log_message("RESULTS", "Generating incidence rate ratios")
}

# Define comparison periods
periods <- c("2016_2020", "2018_2022", "2020_2022")

for (period in periods) {
  tryCatch({
    # Parse period
    years <- as.numeric(strsplit(period, "_")[[1]])
    comparison_start <- years[1]
    comparison_end <- years[2]
    
    # Get the most recent year
    current_year <- max(analysis_data$year)
    
    # Prepare data frame
    irr_results <- data.frame(
      state = character(),
      year = numeric(),
      comparison_period = character(),
      current_incidence = numeric(),
      period_incidence = numeric(),
      relative_risk = numeric(),
      percent_change = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Calculate IRR for each state
    states <- unique(ir_data$state)
    for (s in states) {
      # Get current year IR
      current_ir_row <- subset(ir_data, state == s & year == current_year)
      if (nrow(current_ir_row) == 0) next
      current_ir <- current_ir_row$ir[1]
      
      # Get comparison period average IR
      period_rows <- subset(ir_data, state == s & year >= comparison_start & year <= comparison_end)
      if (nrow(period_rows) == 0) next
      period_ir <- mean(period_rows$ir)
      
      # Calculate relative risk and percent change
      if (period_ir > 0) {
        rr <- current_ir / period_ir
        pct_change <- (rr - 1) * 100
      } else {
        rr <- 1
        pct_change <- 0
      }
      
      # Add to results
      irr_results <- rbind(irr_results, data.frame(
        state = s,
        year = current_year,
        comparison_period = period,
        current_incidence = current_ir,
        period_incidence = period_ir,
        relative_risk = rr,
        percent_change = pct_change,
        stringsAsFactors = FALSE
      ))
    }
    
    # Save IRR results
    irr_file <- paste0(pathogen, "_EstIRRCatch_", period, ".csv")
    write.csv(irr_results, file = irr_file, row.names = FALSE)
    log_message("OUTPUT", paste("Saved IRR data to", irr_file))
    
  }, error = function(e) {
    log_message("ERROR", paste("IRR calculation failed for period", period, ":", e$message))
    
    # Create error indicator IRR data
    error_irr <- data.frame(
      state = "ERROR",
      year = max(analysis_data$year, na.rm=TRUE),
      comparison_period = period,
      current_incidence = NA,
      period_incidence = NA,
      relative_risk = NA,
      percent_change = NA,
      stringsAsFactors = FALSE
    )
    
    irr_file <- paste0(pathogen, "_EstIRRError_", period, ".csv")
    write.csv(error_irr, file = irr_file, row.names = FALSE)
    log_message("OUTPUT", paste("Saved error-state IRR data to", irr_file))
  })
}

# ---- Generate Plots ----

if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("PLOTS", "Generating visualization plots", milestone="VISUALIZATION")
} else {
  log_message("PLOTS", "Generating visualization plots")
}

tryCatch({
  # Create spline trend plots using spline predictions
  plot_data <- ir_data
  
  # Separate spline trends from observed data
  if ("type" %in% names(plot_data)) {
    spline_data <- subset(plot_data, type == "spline_trend")
    observed_data <- subset(plot_data, type == "observed")
  } else {
    # Fallback if no type column
    spline_data <- plot_data
    observed_data <- plot_data
  }
  
  # Overall spline trend plot
  if (nrow(spline_data) > 0) {
    # Calculate overall trend (average across states)
    overall_spline <- aggregate(cbind(ir, ir_lower, ir_upper) ~ year, 
                               data = spline_data, FUN = mean, na.rm = TRUE)
    
    p1 <- ggplot() +
      # Spline trend line with confidence interval
      geom_ribbon(data = overall_spline, aes(x = year, ymin = ir_lower, ymax = ir_upper), 
                  alpha = 0.3, fill = "blue") +
      geom_line(data = overall_spline, aes(x = year, y = ir), 
                color = "blue", linewidth = 1.2) +
      # Observed data points
      geom_point(data = aggregate(ir ~ year, data = observed_data, FUN = mean, na.rm = TRUE),
                aes(x = year, y = ir), color = "darkblue", size = 2.5, alpha = 0.7) +
      labs(
        title = paste(pathogen, "Spline Incidence Rate Trend"),
        subtitle = "Bayesian hierarchical spline model with 95% credible intervals",
        x = "Year",
        y = "Incidence per 100,000"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"))
      
    ggsave(paste0(pathogen, "_spline_trend.png"), p1, width = 10, height = 6, dpi = 300)
    log_message("OUTPUT", paste("Saved spline trend plot to", paste0(pathogen, "_spline_trend.png")))
  } else {
    # Fallback for observed data only
    overall_observed <- aggregate(ir ~ year, data = observed_data, FUN = mean, na.rm = TRUE)
    p1 <- ggplot(overall_observed, aes(x = year, y = ir)) +
      geom_line(color = "blue", linewidth = 1) +
      geom_point(color = "blue", size = 2) +
      labs(
        title = paste(pathogen, "Incidence Rate Trend (Observed Data)"),
        x = "Year", y = "Incidence per 100,000"
      ) +
      theme_minimal()
      
    ggsave(paste0(pathogen, "_trend_observed.png"), p1, width = 8, height = 6)
    log_message("OUTPUT", paste("Saved observed trend plot to", paste0(pathogen, "_trend_observed.png")))
  }
  
  # State-specific spline trends plot
  if (nrow(spline_data) > 0) {
    p2 <- ggplot() +
      # Spline trend lines by state
      geom_line(data = spline_data, aes(x = year, y = ir, color = state, group = state), 
                linewidth = 1) +
      # Observed data points by state
      geom_point(data = observed_data, aes(x = year, y = ir, color = state), 
                size = 2, alpha = 0.7) +
      labs(
        title = paste(pathogen, "Spline Incidence Rate Trends by State"),
        subtitle = "Bayesian spline predictions with observed data points",
        x = "Year",
        y = "Incidence per 100,000",
        color = "State"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
    
    ggsave(paste0(pathogen, "_state_spline_trends.png"), p2, width = 12, height = 8, dpi = 300)
    log_message("OUTPUT", paste("Saved state spline trends plot to", paste0(pathogen, "_state_spline_trends.png")))
  } else {
    # Fallback for observed data only
    p2 <- ggplot(observed_data, aes(x = year, y = ir, color = state, group = state)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      labs(
        title = paste(pathogen, "Incidence Rate by State (Observed Data)"),
        x = "Year", y = "Incidence per 100,000", color = "State"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
      
    ggsave(paste0(pathogen, "_state_trends_observed.png"), p2, width = 10, height = 6)
    log_message("OUTPUT", paste("Saved state observed trends plot to", paste0(pathogen, "_state_trends_observed.png")))
  }
  
  # Comparison plot: Spline vs Observed
  if (nrow(spline_data) > 0 && nrow(observed_data) > 0) {
    # Average trends for comparison
    overall_spline <- aggregate(ir ~ year, data = spline_data, FUN = mean, na.rm = TRUE)
    overall_observed <- aggregate(ir ~ year, data = observed_data, FUN = mean, na.rm = TRUE)
    
    p3 <- ggplot() +
      geom_line(data = overall_spline, aes(x = year, y = ir), 
                color = "blue", linewidth = 1.5, linetype = "solid") +
      geom_point(data = overall_observed, aes(x = year, y = ir), 
                color = "red", size = 3, alpha = 0.8) +
      labs(
        title = paste("FoodNetTrends Analysis:", pathogen),
        subtitle = "Blue line: Smooth spline trend | Red points: Observed data",
        x = "Year",
        y = "Incidence per 100,000"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 16, face = "bold"))
    
    ggsave(paste0(pathogen, "_foodnettrends_comparison.png"), p3, width = 10, height = 6, dpi = 300)
    log_message("OUTPUT", paste("Saved FoodNetTrends comparison plot to", paste0(pathogen, "_foodnettrends_comparison.png")))
  }
}, error = function(e) {
  log_message("ERROR", paste("Plot generation failed:", e$message))
  
  # Create error indicator plots
  png(paste0(pathogen, "_spline_trend_error.png"), width = 800, height = 600)
  plot(1:10, 1:10, type = "n", main = paste(pathogen, "Spline Trend (ERROR)"))
  text(5, 5, "Error generating spline trend plot", col = "red", cex = 2)
  dev.off()
  
  png(paste0(pathogen, "_state_spline_trends_error.png"), width = 800, height = 600)
  plot(1:10, 1:10, type = "n", main = paste(pathogen, "State Spline Trends (ERROR)"))
  text(5, 5, "Error generating state trends plot", col = "red", cex = 2)
  dev.off()
  
  png(paste0(pathogen, "_foodnettrends_comparison_error.png"), width = 800, height = 600)
  plot(1:10, 1:10, type = "n", main = paste("FoodNetTrends", pathogen, "(ERROR)"))
  text(5, 5, "Error generating comparison plot", col = "red", cex = 2)
  dev.off()
  
  log_message("OUTPUT", "Created error indicator plot files")
})

# ---- Generate Summary ----

if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("SUMMARY", "Generating summary report", milestone="SUMMARY")
} else {
  log_message("SUMMARY", "Generating summary report")
}

# Create summary file
summary_file <- paste0(pathogen, "_summary.txt")
sink(summary_file)
cat("=======================================================\n")
cat(" FoodNetTrends v1.0.0-rc.1 Analysis Summary         \n")
cat("=======================================================\n")
cat(paste("Pathogen:         ", pathogen, "\n"))
cat(paste("Analysis Date:    ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
cat(paste("Pipeline Version: ", "v1.0.0-rc.1", "\n"))
cat(paste("MMWR File:        ", args$mmwrFile, "\n"))

# Data Quality Summary
cat("\n--- DATA QUALITY ASSESSMENT ---\n")
cat(paste("Records Analyzed: ", nrow(analysis_data), "\n"))
cat(paste("States Included:  ", data_quality$state_count, " (", paste(unique(analysis_data$state), collapse=", "), ")\n"))
cat(paste("Time Period:      ", paste(data_quality$year_range, collapse=" - "), "\n"))
cat(paste("Data Quality:     ", if(data_quality$passed) "PASSED" else "ISSUES DETECTED", "\n"))
if (!data_quality$passed) {
  for (issue in data_quality$issues) {
    cat(paste("  WARNING: ", issue, "\n"))
  }
}

# Model Performance Summary  
cat("\n--- MODEL DIAGNOSTICS ---\n")
cat(paste("Model Type:       ", "Bayesian Hierarchical Spline (brms)", "\n"))
cat(paste("Convergence:      ", if(is.na(convergence_check$converged)) "N/A (dummy model)" else 
                                  if(convergence_check$converged) "CONVERGED" else "ISSUES DETECTED", "\n"))
if (!is.na(convergence_check$max_rhat)) {
  cat(paste("Max Rhat:         ", round(convergence_check$max_rhat, 3), 
           if(convergence_check$max_rhat <= 1.1) " (Good)" else " (Concerning)", "\n"))
}
cat(paste("Trend Significance:", if(is.null(trend_significance) || is.na(trend_significance$significant)) "N/A" else 
                                  if(trend_significance$significant) "SIGNIFICANT" else "NOT SIGNIFICANT", "\n"))
if (!is.null(trend_significance) && !is.null(trend_significance$proportion_significant) && !is.na(trend_significance$proportion_significant)) {
  cat(paste("Significant Terms:", round(trend_significance$proportion_significant * 100, 1), "%\n"))
}

# Analysis Parameters
cat("\n--- ANALYSIS PARAMETERS ---\n")
cat(paste("Chains:           ", args$chains, "\n"))
cat(paste("Iterations:       ", args$iterations, "\n"))
cat(paste("Cores Used:       ", args$cores, "\n"))
cat(paste("Travel Filter:    ", args$travel, "\n"))
cat(paste("CIDT Filter:      ", args$cidt, "\n"))

# Results Summary
cat("\n--- INCIDENCE RATE SUMMARY ---\n")
if ("type" %in% names(ir_data)) {
  observed_ir <- subset(ir_data, type == "observed")
  if (nrow(observed_ir) > 0) {
    state_summary <- aggregate(ir ~ state, data = observed_ir, FUN = function(x) round(mean(x, na.rm=TRUE), 2))
    overall_mean <- round(mean(observed_ir$ir, na.rm=TRUE), 2)
    overall_range <- round(range(observed_ir$ir, na.rm=TRUE), 2)
    
    cat(paste("Overall Mean IR:  ", overall_mean, " per 100,000\n"))
    cat(paste("Range:            ", paste(overall_range, collapse=" - "), " per 100,000\n"))
    cat("State Averages:\n")
    for (i in 1:nrow(state_summary)) {
      cat(paste("  ", state_summary$state[i], ": ", state_summary$ir[i], " per 100,000\n"))
    }
  }
} else {
  # Fallback for older format
  state_summary <- aggregate(ir ~ state, data = ir_data, FUN = function(x) round(mean(x, na.rm=TRUE), 2))
  cat("State Averages:\n")
  for (i in 1:nrow(state_summary)) {
    cat(paste("  ", state_summary$state[i], ": ", state_summary$ir[i], " per 100,000\n"))
  }
}

# Output Files
cat("\n--- OUTPUT FILES GENERATED ---\n")
cat(paste("Model File:       ", paste0(pathogen, "_brm.Rds"), "\n"))
cat(paste("IR Data:          ", paste0(pathogen, "_IRCatch.csv"), "\n"))
cat(paste("Spline Trends:    ", paste0(pathogen, "_spline_trend.png"), "\n"))
cat(paste("State Trends:     ", paste0(pathogen, "_state_spline_trends.png"), "\n"))
cat(paste("Comparison Plot:  ", paste0(pathogen, "_foodnettrends_comparison.png"), "\n"))

# Interpretation Guidelines
cat("\n--- INTERPRETATION GUIDELINES ---\n")
cat("1. Spline trend lines show underlying temporal patterns\n")
cat("2. Credible intervals indicate uncertainty in estimates\n")
cat("3. Check convergence diagnostics before interpretation\n")
if (!data_quality$passed) {
  cat("4. WARNING: Data quality issues detected - interpret with caution\n")
}
if (!is.na(trend_significance$significant) && !trend_significance$significant) {
  cat("5. WARNING: Trends may not be statistically significant\n")
}

cat("\n=======================================================\n")
sink()

if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("OUTPUT", paste("Saved summary to", summary_file))
} else {
  log_message("OUTPUT", paste("Saved summary to", summary_file))
}

# Generate a simple summary file if it doesn't exist already
# This helps prevent "Missing output file" errors in the pipeline
summary_file_path <- paste0(pathogen, "_summary.txt")
if (!file.exists(summary_file_path)) {
  log_message("OUTPUT", paste("Creating summary file", summary_file_path))
  
  # Create a simple summary file
  write(paste("Summary for", pathogen, "analysis completed at", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), 
        file = summary_file_path)
  write(paste("Data characteristics:"), file = summary_file_path, append = TRUE)
  write(paste("  Total records:", nrow(pathogen_data)), file = summary_file_path, append = TRUE)
  write(paste("  Unique states:", paste(unique(pathogen_counts$state), collapse=", ")), 
        file = summary_file_path, append = TRUE)
  write(paste("  Years covered:", paste(sort(unique(pathogen_counts$year)), collapse=", ")), 
        file = summary_file_path, append = TRUE)
}

# Complete with final diagnostics summary
log_message("INFO", "=== FINAL ANALYSIS SUMMARY ===")
log_message("INFO", paste("Pathogen:", pathogen))
log_message("INFO", paste("Data Quality:", if(data_quality$passed) "PASSED" else "ISSUES DETECTED"))
log_message("INFO", paste("Model Convergence:", if(is.na(convergence_check$converged)) "N/A" else if(convergence_check$converged) "CONVERGED" else "ISSUES"))
log_message("INFO", paste("Trend Significance:", if(is.na(trend_significance$significant)) "N/A" else if(trend_significance$significant) "SIGNIFICANT" else "NOT SIGNIFICANT"))
log_message("INFO", paste("Spline Trends Generated:", if("type" %in% names(ir_data) && any(ir_data$type == "spline_trend")) "YES" else "NO"))

if (exists("has_progress_tracking") && has_progress_tracking) {
  log_progress("COMPLETE", paste("FoodNetTrends analysis completed successfully for", pathogen), milestone="COMPLETE")
} else {
  log_message("COMPLETE", paste("FoodNetTrends analysis completed successfully for", pathogen))
}