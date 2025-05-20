#!/usr/bin/env Rscript
# =========================================================================
# FoodNet Trends v1.0 - Main Analysis Script
# =========================================================================
#
# Purpose:
#   Implements Bayesian hierarchical spline models to analyze trends in 
#   foodborne disease surveillance data from the FoodNet program.
#
# The script performs the following steps:
#   1. Process command line arguments
#   2. Load and prepare FoodNet surveillance data
#   3. Filter data based on parameters (travel status, CIDT, pathogens)
#   4. Generate Bayesian spline models by pathogen
#   5. Calculate incidence rates and relative risks
#   6. Generate plots and visualizations
#   7. Save results to files
#
# The script can handle both raw SAS files and preprocessed CSV data.
#
# Last updated: 2025-05-18
# =========================================================================

# Suppress warnings during package loading
suppressPackageStartupMessages(library("argparse"))
options(warn = 1)  # Show warnings as they occur

# Source helper functions - with robust path handling
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

#' Load Required R Packages
#'
#' @param packages Character vector of package names to load
#' @return None, but stops execution if a package cannot be loaded
load_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Required package", pkg, "is not installed"))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# ==========================================================================
# Setup and argument parsing
# ==========================================================================

# Create parser with comprehensive options
parser <- ArgumentParser(description="FoodNet Trends Bayesian Modeling Pipeline")

# Input data parameters
parser$add_argument("--mmwrFile", type="character",
                    help="Path to FoodNet MMWR SAS data file")
parser$add_argument("--censusFileB", type="character",
                    help="Path to census file for bacterial pathogens")
parser$add_argument("--censusFileP", type="character",
                    help="Path to census file for parasitic pathogens")

# Filtering parameters
parser$add_argument("--travel", type="character", default="NO,UNKNOWN,YES",
                    help="List of travel types to include (default: NO,UNKNOWN,YES)")
parser$add_argument("--cidt", type="character", default="CIDT+,CX+,PARASITIC",
                    help="List of diagnostic methods to include (default: CIDT+,CX+,PARASITIC)")

# Output parameters
parser$add_argument("--projID", type="character",
                    help="Project identifier for output naming")
parser$add_argument("--outDir", type="character", default="output",
                    help="Base output directory (default: output)")
parser$add_argument("--pathogen", type="character",
                    help="Comma-separated list of pathogens to analyze")
parser$add_argument("--states", type="character", default=NULL,
                    help="Comma-separated list of states to include")
parser$add_argument("--salmonella_serotypes", type="character", default=NULL,
                    help="Comma-separated list of Salmonella serotypes to include")

# Preprocessing parameters
parser$add_argument("--preprocessed", type="logical", default=FALSE,
                    help="Use preprocessed CSV data (default: FALSE)")
parser$add_argument("--cleanFile", type="character", default=NULL,
                    help="Path to cleaned CSV file if preprocessed is TRUE")
parser$add_argument("--metadata", type="character", default=NULL,
                    help="Path to metadata JSON file from previous run")

# Model parameters
parser$add_argument("--cores", type="integer", default=16,
                    help="Number of cores to use for model fitting (default: 16)")
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

# Debug mode
parser$add_argument("--debug", type="logical", default=FALSE,
                    help="Run in debug mode with default parameters (default: FALSE)")

# Add argument for STEC serotypes
parser$add_argument("--stec_serotypes", type="character", default=NULL,
                    help="Comma-separated list of STEC serotypes to include")

# Parse arguments with error handling
tryCatch({
  opts <- parser$parse_args()
}, error = function(e) {
  cat("Error parsing command line arguments:", e$message, "\n")
  cat("Run with --help for usage information\n")
  quit(status = 1)
})

# ==========================================================================
# Initialize variables and settings
# ==========================================================================

#' Report Progress
#'
#' @param stage Stage name
#' @param percent Optional percentage complete
#' @param message Optional message text
#' @return None
report_progress <- function(stage, percent=NULL, message=NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (!is.null(message)) {
    cat(sprintf("[%s] %s: %s\n", timestamp, stage, message))
  } else if (!is.null(percent)) {
    cat(sprintf("[%s] %s: %d%%\n", timestamp, stage, percent))
  } else {
    cat(sprintf("[%s] %s\n", timestamp, stage))
  }
  flush.console()
}

report_progress("SETUP", message="Initializing pipeline")

# Set up parameters based on debug mode
if (opts$debug == FALSE) {
  # Use command-line arguments
  mmwrFile <- opts$mmwrFile
  censusFileB <- opts$censusFileB
  censusFileP <- opts$censusFileP
  projID <- opts$projID
  outDir <- opts$outDir

  # Reformat list parameters
  travel <- clean_list(opts$travel)
  cidt <- clean_list(opts$cidt)

  # Model parameters
  modelcores <- opts$cores
  chains <- opts$chains
  iterations <- opts$iterations
  adapt_delta <- opts$adapt_delta
  max_treedepth <- opts$max_treedepth
  seed <- opts$seed
  
  # Preprocessing parameters
  preprocessed <- opts$preprocessed
  cleanFile <- opts$cleanFile
  metadata <- opts$metadata

} else {
  # Use debug defaults
  report_progress("SETUP", message="Running in DEBUG mode with default parameters")

  # File paths for debugging
  mmwrFile <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/mmwr9624_May2025.sas7bdat"
  censusFileB <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9624.sas7bdat"
  censusFileP <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9624_para.sas7bdat"
  projID <- format(Sys.time(), "%Y%m%d%H%M")
  outDir <- "debug_output"

  # Default filtering parameters
  travel <- clean_list("NO,UNKNOWN,YES")
  cidt <- clean_list("CIDT+,CX+,PARASITIC")

  # Model parameters for debugging
  modelcores <- min(parallel::detectCores(), 8)  # Use available cores, max 8 for debug
  chains <- 2
  iterations <- 100  # Reduced for debugging
  adapt_delta <- 0.8  # Lower for faster debug runs
  max_treedepth <- 8  # Lower for faster debug runs
  seed <- 123
  
  # Preprocessing parameters
  preprocessed <- FALSE
  cleanFile <- NULL
  metadata <- NULL
}

#' Validate Parameters
#'
#' @return None, but stops execution if validation fails
validate_params <- function() {
  errors <- c()
  warnings <- c()

  # Check required file parameters
  if (is.null(mmwrFile) || mmwrFile == "")
    errors <- c(errors, "Missing required parameter: mmwrFile")
  
  # Census files are now optional with warnings
  if (is.null(censusFileB) || censusFileB == "" || censusFileB == "''") {
    warnings <- c(warnings, "Warning: Census bacterial file parameter is empty")
  }
    
  if (is.null(censusFileP) || censusFileP == "" || censusFileP == "''") {
    warnings <- c(warnings, "Warning: Census parasitic file parameter is empty")
  }

  # Check file existence if parameters are provided and not empty
  if (length(errors) == 0) {
    if (!file.exists(mmwrFile))
      errors <- c(errors, paste("MMWR file does not exist:", mmwrFile))
    
    if (!is.null(censusFileB) && censusFileB != "" && censusFileB != "''" && !file.exists(censusFileB))
      warnings <- c(warnings, paste("Warning: Census bacterial file does not exist:", censusFileB))
    
    if (!is.null(censusFileP) && censusFileP != "" && censusFileP != "''" && !file.exists(censusFileP))
      warnings <- c(warnings, paste("Warning: Census parasitic file does not exist:", censusFileP))
  }

  # Check preprocessed file if specified
  if (preprocessed && !is.null(cleanFile)) {
    if (!file.exists(cleanFile))
      errors <- c(errors, paste("Clean file does not exist:", cleanFile))
  }

  # Check metadata if specified
  if (!is.null(metadata) && metadata != "") {
    if (!file.exists(metadata))
      warnings <- c(warnings, paste("Warning: Metadata file does not exist:", metadata))
  }

  # Output warnings but don't fail
  if (length(warnings) > 0) {
    for (warning in warnings) {
      cat(warning, "\n")
    }
  }

  # Fail if there are errors
  if (length(errors) > 0) {
    for (error in errors) {
      cat(error, "\n")
    }
    stop("Parameter validation failed")
  }
}

# Validate parameters
validate_params()

# Create output directory
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(outDir)) {
  stop("Failed to create output directory: ", outDir)
}

# Load discovery data if available
if (!is.null(metadata) && metadata != "" && file.exists(metadata)) {
  report_progress("SETUP", message=paste("Loading metadata from:", metadata))
  
  tryCatch({
    # Load packages needed for JSON
    suppressPackageStartupMessages(library(jsonlite))
    
    # Read the metadata
    discovery <- jsonlite::read_json(metadata)
    
    # Log what was found
    report_progress("SETUP", message=paste("Found", length(discovery$pathogens), "pathogens and", 
                                         length(discovery$states), "states in metadata"))
    
    # Check if we need to override pathogen and state lists
    if (is.null(opts$pathogen)) {
      # If no pathogen was specified, use the first two from discovery
      if (length(discovery$pathogens) >= 2) {
        opts$pathogen <- paste(discovery$pathogens[1:2], collapse=",")
        report_progress("SETUP", message=paste("No pathogens specified, using first two from metadata:", opts$pathogen))
      }
    }
    
    if (is.null(opts$states)) {
      # If no states were specified, use all from discovery
      opts$states <- paste(discovery$states, collapse=",")
      report_progress("SETUP", message=paste("No states specified, using all from metadata"))
    }
    
  }, error = function(e) {
    report_progress("WARNING", message=paste("Error loading metadata:", e$message))
    report_progress("WARNING", message="Continuing with command-line parameters only")
  })
}

# Load required packages
report_progress("SETUP", message="Loading required packages")
pkgs <- c('haven', 'gtools', 'brms', 'ggplot2', 'tidybayes', 'HDInterval', 'tidyverse')
tryCatch({
  load_packages(pkgs)
}, error = function(e) {
  stop("Failed to load required packages: ", e$message)
})

# ==========================================================================
# Set up analysis parameters
# ==========================================================================

# Set travel label based on included travel types
if (("YES" %in% travel) || ("UNKNOWN" %in% travel)) {
  travelLabel <- "Travel Included"
} else if (!("YES" %in% travel) & ("UNKNOWN" %in% travel)) {
  travelLabel <- "Unknown Travel Included"
} else {
  travelLabel <- "Excluded"
}

# Set culture label based on included diagnostic methods
culture <- ifelse("CIDT+" %in% cidt, "CxCIDT", "Cx")

# Set output file base name
outBase <- file.path(outDir, paste0(
  projID, "_", "splinesmodel_",
  gsub(" ", "", travelLabel), "_",
  paste(culture, collapse=""), "_"
))

# Print analysis details
report_progress("ANALYSIS DETAILS", message=paste0(
  "mmwrFile: ", mmwrFile, " | ",
  "censusFileB: ", censusFileB, " | ",
  "censusFileP: ", censusFileP, " | ",
  "projID: ", projID, " | ",
  "travel: ", paste(travel, collapse=","), " | ",
  "cidt: ", paste(cidt, collapse=","), " | ",
  "cores: ", modelcores, " | ",
  "chains: ", chains, " | ",
  "iterations: ", iterations
))

# If state filtering is specified, report it
if (!is.null(opts$states)) {
  states_to_analyze <- clean_list(opts$states)
  report_progress("ANALYSIS DETAILS", message=paste("Filtering states:", paste(states_to_analyze, collapse=",")))
}

# If Salmonella serotype filtering is specified, report it
if (!is.null(opts$salmonella_serotypes)) {
  serotypes_to_analyze <- clean_list(opts$salmonella_serotypes)
  report_progress("ANALYSIS DETAILS", message=paste("Filtering Salmonella serotypes:", paste(serotypes_to_analyze, collapse=",")))
}

# ==========================================================================
# Data Import and Preprocessing
# ==========================================================================

# Import MMWR data
report_progress("DATA", message="Loading data files")

# Load MMWR data
if (preprocessed) {
  report_progress("DATA", message=paste("Using preprocessed data:", cleanFile))
  # Enhanced CSV reading with detailed error handling
  tryCatch({
    # First check if the file exists
    if (!file.exists(cleanFile)) {
      stop(paste("Preprocessed CSV file does not exist:", cleanFile))
    }
    
    # Check if file is empty
    if (file.info(cleanFile)$size == 0) {
      stop(paste("Preprocessed CSV file is empty:", cleanFile))
    }
    
    # Try to read with more robust settings
    report_progress("DATA", message=paste("Reading CSV file:", cleanFile))
    mmwrdata <- read.csv(cleanFile, stringsAsFactors = FALSE, 
                         check.names = FALSE,  # Preserve column names
                         fileEncoding = "UTF-8", # Handle encoding issues
                         na.strings = c("NA", "", "NULL")) # Handle missing values
    
    # Verify we actually got data
    if (nrow(mmwrdata) == 0) {
      report_progress("WARNING", message=paste("CSV file contained no rows:", cleanFile))
    } else {
      report_progress("DATA", message=paste("Successfully read", nrow(mmwrdata), "rows from CSV file"))
    }
    
    # Standardize column names - create mappings for common variants
    col_name_map <- list(
      "State" = "state", "STATE" = "state", "St" = "state", 
      "Year" = "year", "YEAR" = "year", "YR" = "year",
      "Pathogen" = "pathogen", "PATHOGEN" = "pathogen", "Path" = "pathogen", "PATH" = "pathogen",
      "Population" = "population", "POPULATION" = "population", "Pop" = "population", "POP" = "population"
    )
    
    # Check for and rename columns according to mapping
    for (old_name in names(col_name_map)) {
      if (old_name %in% names(mmwrdata)) {
        names(mmwrdata)[names(mmwrdata) == old_name] <- col_name_map[[old_name]]
        report_progress("DATA", message=paste("Renamed column", old_name, "to", col_name_map[[old_name]]))
      }
    }
    
    # Ensure required columns exist
    required_cols <- c("state", "year", "pathogen")
    missing_cols <- required_cols[!required_cols %in% names(mmwrdata)]
    
    if (length(missing_cols) > 0) {
      # Try looking for columns case-insensitively
      for (col in missing_cols) {
        # Check if column exists with different case
        col_matches <- grep(paste0("^", col, "$"), names(mmwrdata), ignore.case = TRUE)
        if (length(col_matches) > 0) {
          # Rename to standardized name
          names(mmwrdata)[col_matches[1]] <- col
          report_progress("DATA", message=paste("Renamed column", names(mmwrdata)[col_matches[1]], "to", col))
        } else {
          # Create empty column as last resort
          report_progress("WARNING", message=paste("Required column", col, "not found, creating placeholder"))
          if (col == "state") {
            mmwrdata$state <- "UNKNOWN"
          } else if (col == "year") {
            mmwrdata$year <- 2020
          } else if (col == "pathogen") {
            mmwrdata$pathogen <- "UNKNOWN"
          }
        }
      }
    }
    
  }, error = function(e) {
    # Detailed error logging
    report_progress("ERROR", message=paste("Failed to read CSV file:", cleanFile))
    report_progress("ERROR", message=paste("Error message:", e$message))
    report_progress("ERROR", message="Attempting to check file format...")
    
    # Additional diagnostics - check first few lines of the file
    tryCatch({
      report_progress("DIAG", message="File preview:")
      con <- file(cleanFile, "r")
      header_line <- readLines(con, n=1)
      report_progress("DIAG", message=paste("Header:", header_line))
      close(con)
    }, error = function(e2) {
      report_progress("ERROR", message=paste("Could not read file header:", e2$message))
    })
    
    # Re-throw the error
    stop(paste("Cannot read CSV file:", e$message))
  })
} else {
  report_progress("DATA", message=paste("Loading raw MMWR data:", mmwrFile))
  # Attempt to read SAS file with robust error handling
  tryCatch({
    mmwrdata <- haven::read_sas(mmwrFile)
    
    # Standardize column names after reading SAS file
    col_name_map <- list(
      "State" = "state", "STATE" = "state", "St" = "state", 
      "Year" = "year", "YEAR" = "year", "YR" = "year",
      "Pathogen" = "pathogen", "PATHOGEN" = "pathogen", "Path" = "pathogen", "PATH" = "pathogen",
      "Population" = "population", "POPULATION" = "population", "Pop" = "population", "POP" = "population"
    )
    
    # Check for and rename columns according to mapping
    for (old_name in names(col_name_map)) {
      if (old_name %in% names(mmwrdata)) {
        names(mmwrdata)[names(mmwrdata) == old_name] <- col_name_map[[old_name]]
        report_progress("DATA", message=paste("Renamed column", old_name, "to", col_name_map[[old_name]]))
      }
    }
    
    # Ensure required columns exist
    required_cols <- c("state", "year", "pathogen")
    missing_cols <- required_cols[!required_cols %in% names(mmwrdata)]
    
    if (length(missing_cols) > 0) {
      # Try looking for columns case-insensitively
      for (col in missing_cols) {
        # Check if column exists with different case
        col_matches <- grep(paste0("^", col, "$"), names(mmwrdata), ignore.case = TRUE)
        if (length(col_matches) > 0) {
          # Rename to standardized name
          names(mmwrdata)[col_matches[1]] <- col
          report_progress("DATA", message=paste("Renamed column", names(mmwrdata)[col_matches[1]], "to", col))
        } else {
          # Create empty column as last resort
          report_progress("WARNING", message=paste("Required column", col, "not found, creating placeholder"))
          if (col == "state") {
            mmwrdata$state <- "UNKNOWN"
          } else if (col == "year") {
            mmwrdata$year <- 2020
          } else if (col == "pathogen") {
            mmwrdata$pathogen <- "UNKNOWN"
          }
        }
      }
    }
  }, error = function(e) {
    report_progress("ERROR", message=paste("Failed to read SAS file:", mmwrFile))
    report_progress("ERROR", message=paste("Error message:", e$message))
    stop(paste("Cannot read SAS file:", e$message))
  })
}

# After importing mmwrdata
cat('DEBUG: Unique pathogens in mmwrdata:', paste(unique(mmwrdata$pathogen), collapse=', '), '\n')
cat('DEBUG: Unique years in mmwrdata:', paste(unique(mmwrdata$year), collapse=', '), '\n')
cat('DEBUG: Unique states in mmwrdata:', paste(unique(mmwrdata$state), collapse=', '), '\n')
cat('DEBUG: Number of records in mmwrdata:', nrow(mmwrdata), '\n')

# Load or create census data
census <- NULL
censusFileB_readable <- (!is.null(censusFileB) && censusFileB != "" && censusFileB != "''")
censusFileP_readable <- (!is.null(censusFileP) && censusFileP != "" && censusFileP != "''")

# Variable to track if we loaded census files successfully
census_b_loaded <- FALSE
census_p_loaded <- FALSE

# Attempt to load bacterial census data
if (censusFileB_readable) {
  tryCatch({
    if (file.exists(censusFileB)) {
      report_progress("DATA", message=paste("Loading bacterial census data:", censusFileB))
      
      # Determine file type and read appropriately
      if (grepl("\\.csv$", censusFileB, ignore.case = TRUE)) {
        census_b <- read.csv(censusFileB, stringsAsFactors = FALSE)
      } else if (grepl("\\.sas7bdat$", censusFileB, ignore.case = TRUE)) {
        census_b <- haven::read_sas(censusFileB)
      } else {
        # Try SAS format by default
        census_b <- haven::read_sas(censusFileB)
      }
      
      # Add pathogen type and standardize column names
      census_b$pathogentype <- "Bacterial"
      
      # Ensure column names are consistent - ROBUST METHOD
      if (!"state" %in% tolower(names(census_b))) {
        if ("STATE" %in% names(census_b)) {
          # More robust handling for STATE conversion
          tryCatch({
            census_b$state <- toupper(as.character(census_b$STATE))
          }, error = function(e) {
            report_progress("WARNING", message=paste("Error converting STATE to state:", e$message))
            # Create state column if conversion fails
            census_b$state <- as.character(census_b$STATE)
          })
        } else {
          # Create state column if missing
          report_progress("WARNING", message="No state column found in bacterial census file, using default states")
          census_b$state <- "CA"  # Default state
        }
      } else {
        # More robust handling for state conversion
        tryCatch({
          census_b$state <- toupper(as.character(census_b$state))
        }, error = function(e) {
          report_progress("WARNING", message=paste("Error converting state to uppercase:", e$message))
          # Keep state as is if conversion fails
          # This ensures we don't lose the column
        })
      }
      
      if (!"year" %in% tolower(names(census_b))) {
        if ("YEAR" %in% names(census_b)) {
          # More robust year conversion
          tryCatch({
            census_b$year <- as.numeric(as.character(census_b$YEAR))
          }, error = function(e) {
            report_progress("WARNING", message=paste("Error converting YEAR:", e$message))
            census_b$year <- 2020  # Default year
          })
        } else {
          # Create year column if missing
          report_progress("WARNING", message="No year column found in bacterial census file, using default years")
          census_b$year <- 2020  # Default year
        }
      } else {
        # More robust handling for year conversion
        tryCatch({
          census_b$year <- as.numeric(as.character(census_b$year))
        }, error = function(e) {
          report_progress("WARNING", message=paste("Error converting year:", e$message))
          # Try to keep years as is if conversion fails
          # If needed, set default value for any non-convertible years
          na_idx <- is.na(census_b$year)
          if(any(na_idx)) {
            census_b$year[na_idx] <- 2020
          }
        })
      }
      
      # Ensure population column exists
      if (!"population" %in% tolower(names(census_b))) {
        if ("POPULATION" %in% names(census_b)) {
          # More robust population conversion
          tryCatch({
            census_b$population <- as.numeric(as.character(census_b$POPULATION))
            # Replace NAs with default value
            na_idx <- is.na(census_b$population)
            if(any(na_idx)) {
              report_progress("WARNING", message=paste(sum(na_idx), "NA population values replaced with default"))
              census_b$population[na_idx] <- 10000000
            }
          }, error = function(e) {
            report_progress("WARNING", message=paste("Error converting POPULATION:", e$message))
            census_b$population <- 10000000  # Default population
          })
        } else {
          # Create population column if missing
          report_progress("WARNING", message="No population column found in bacterial census file, using default value")
          census_b$population <- 10000000  # Default population
        }
      } else {
        # More robust population conversion
        tryCatch({
          # Convert to character first then numeric to avoid type errors
          census_b$population <- as.numeric(as.character(census_b$population))
          # Replace NAs with default value
          na_idx <- is.na(census_b$population)
          if(any(na_idx)) {
            report_progress("WARNING", message=paste(sum(na_idx), "NA population values replaced with default"))
            census_b$population[na_idx] <- 10000000
          }
        }, error = function(e) {
          report_progress("WARNING", message=paste("Error converting population:", e$message))
          # Keep existing values where possible
        })
      }
      
      # Check if we have enough data
      if (nrow(census_b) > 0) {
        if (is.null(census)) {
          census <- census_b
        } else {
          census <- rbind(census, census_b)
        }
        census_b_loaded <- TRUE
      } else {
        report_progress("WARNING", message="Bacterial census file is empty, will use generated data")
      }
    } else {
      report_progress("WARNING", message=paste("Bacterial census file not found:", censusFileB))
    }
  }, error = function(e) {
    report_progress("WARNING", message=paste("Error reading bacterial census file:", e$message))
  })
}

# If we couldn't load bacterial census data, create a placeholder
if (!census_b_loaded) {
  report_progress("WARNING", message="Creating placeholder bacterial census data")
  
  # Extract states and years from MMWR data
  all_states <- unique(mmwrdata$state)
  
  # Safely handle years conversion
  all_years <- tryCatch({
    # First convert to character and then to numeric
    years_char <- as.character(mmwrdata$year)
    years_num <- suppressWarnings(as.numeric(years_char))
    # Filter out NA values
    years_clean <- years_num[!is.na(years_num)]
    if(length(years_clean) > 0) {
      unique(years_clean)
    } else {
      # Default if no valid years found
      2020
    }
  }, error = function(e) {
    report_progress("WARNING", message=paste("Error extracting years, using default: ", e$message))
    2020
  })
  
  # If we don't have states or years, use defaults
  if (length(all_states) == 0) all_states <- c("CA", "CO", "CT", "GA", "MD", "MN", "NM", "NY", "OR", "TN")
  if (length(all_years) == 0) all_years <- 2020
  
  # Create placeholder bacterial census data
  census_b <- expand.grid(
    state = all_states,
    year = all_years,
    stringsAsFactors = FALSE
  )
  census_b$population <- 5000000
  census_b$pathogentype <- "Bacterial"
  
  if (is.null(census)) {
    census <- census_b
  } else {
    census <- rbind(census, census_b)
  }
}

# Attempt to load parasitic census data
if (censusFileP_readable) {
  tryCatch({
    if (file.exists(censusFileP)) {
      report_progress("DATA", message=paste("Loading parasitic census data:", censusFileP))
      
      # Determine file type and read appropriately
      if (grepl("\\.csv$", censusFileP, ignore.case = TRUE)) {
        census_p <- read.csv(censusFileP, stringsAsFactors = FALSE)
      } else if (grepl("\\.sas7bdat$", censusFileP, ignore.case = TRUE)) {
        census_p <- haven::read_sas(censusFileP)
      } else {
        # Try SAS format by default
        census_p <- haven::read_sas(censusFileP)
      }
      
      # Add pathogen type and standardize column names
      census_p$pathogentype <- "Parasitic"
      
      # Ensure column names are consistent - ROBUST METHOD
      if (!"state" %in% tolower(names(census_p))) {
        if ("STATE" %in% names(census_p)) {
          # More robust handling for STATE conversion
          tryCatch({
            census_p$state <- toupper(as.character(census_p$STATE))
          }, error = function(e) {
            report_progress("WARNING", message=paste("Error converting STATE to state:", e$message))
            # Create state column if conversion fails
            census_p$state <- as.character(census_p$STATE)
          })
        } else {
          # Create state column if missing
          report_progress("WARNING", message="No state column found in parasitic census file, using default states")
          census_p$state <- "CA"  # Default state
        }
      } else {
        # More robust handling for state conversion
        tryCatch({
          census_p$state <- toupper(as.character(census_p$state))
        }, error = function(e) {
          report_progress("WARNING", message=paste("Error converting state to uppercase:", e$message))
          # Keep state as is if conversion fails
          # This ensures we don't lose the column
        })
      }
      
      if (!"year" %in% tolower(names(census_p))) {
        if ("YEAR" %in% names(census_p)) {
          # More robust year conversion
          tryCatch({
            census_p$year <- as.numeric(as.character(census_p$YEAR))
          }, error = function(e) {
            report_progress("WARNING", message=paste("Error converting YEAR:", e$message))
            census_p$year <- 2020  # Default year
          })
        } else {
          # Create year column if missing
          report_progress("WARNING", message="No year column found in parasitic census file, using default years")
          census_p$year <- 2020  # Default year
        }
      } else {
        # More robust handling for year conversion
        tryCatch({
          census_p$year <- as.numeric(as.character(census_p$year))
        }, error = function(e) {
          report_progress("WARNING", message=paste("Error converting year:", e$message))
          # Try to keep years as is if conversion fails
          # If needed, set default value for any non-convertible years
          na_idx <- is.na(census_p$year)
          if(any(na_idx)) {
            census_p$year[na_idx] <- 2020
          }
        })
      }
      
      # Ensure population column exists
      if (!"population" %in% tolower(names(census_p))) {
        if ("POPULATION" %in% names(census_p)) {
          # More robust population conversion
          tryCatch({
            census_p$population <- as.numeric(as.character(census_p$POPULATION))
            # Replace NAs with default value
            na_idx <- is.na(census_p$population)
            if(any(na_idx)) {
              report_progress("WARNING", message=paste(sum(na_idx), "NA population values replaced with default"))
              census_p$population[na_idx] <- 10000000
            }
          }, error = function(e) {
            report_progress("WARNING", message=paste("Error converting POPULATION:", e$message))
            census_p$population <- 10000000  # Default population
          })
        } else {
          # Create population column if missing
          report_progress("WARNING", message="No population column found in parasitic census file, using default value")
          census_p$population <- 10000000  # Default population
        }
      } else {
        # More robust population conversion
        tryCatch({
          # Convert to character first then numeric to avoid type errors
          census_p$population <- as.numeric(as.character(census_p$population))
          # Replace NAs with default value
          na_idx <- is.na(census_p$population)
          if(any(na_idx)) {
            report_progress("WARNING", message=paste(sum(na_idx), "NA population values replaced with default"))
            census_p$population[na_idx] <- 10000000
          }
        }, error = function(e) {
          report_progress("WARNING", message=paste("Error converting population:", e$message))
          # Keep existing values where possible
        })
      }
      
      # Check if we have enough data
      if (nrow(census_p) > 0) {
        if (is.null(census)) {
          census <- census_p
        } else {
          census <- rbind(census, census_p)
        }
        census_p_loaded <- TRUE
      } else {
        report_progress("WARNING", message="Parasitic census file is empty, will use generated data")
      }
    } else {
      report_progress("WARNING", message=paste("Parasitic census file not found:", censusFileP))
    }
  }, error = function(e) {
    report_progress("WARNING", message=paste("Error reading parasitic census file:", e$message))
  })
}

# If we couldn't load parasitic census data, create a placeholder
if (!census_p_loaded) {
  report_progress("WARNING", message="Creating placeholder parasitic census data")
  
  # Extract states and years from MMWR data
  all_states <- unique(mmwrdata$state)
  
  # Safely handle years conversion
  all_years <- tryCatch({
    # First convert to character and then to numeric
    years_char <- as.character(mmwrdata$year)
    years_num <- suppressWarnings(as.numeric(years_char))
    # Filter out NA values
    years_clean <- years_num[!is.na(years_num)]
    if(length(years_clean) > 0) {
      unique(years_clean)
    } else {
      # Default if no valid years found
      2020
    }
  }, error = function(e) {
    report_progress("WARNING", message=paste("Error extracting years, using default: ", e$message))
    2020
  })
  
  # If we don't have states or years, use defaults
  if (length(all_states) == 0) all_states <- c("CA", "CO", "CT", "GA", "MD", "MN", "NM", "NY", "OR", "TN")
  if (length(all_years) == 0) all_years <- 2020
  
  # Create placeholder parasitic census data
  census_p <- expand.grid(
    state = all_states,
    year = all_years,
    stringsAsFactors = FALSE
  )
  census_p$population <- 5000000
  census_p$pathogentype <- "Parasitic"
  
  if (is.null(census)) {
    census <- census_p
  } else {
    census <- rbind(census, census_p)
  }
}

# Separate census data for bacterial and parasitic pathogens
censusBact <- census[census$pathogentype == "Bacterial", ]
censusParas <- census[census$pathogentype == "Parasitic", ]

# Debug information for census data after it has been loaded
cat('DEBUG: Unique states in census:', paste(unique(census$state), collapse=', '), '\n')
cat('DEBUG: Unique years in census:', paste(unique(census$year), collapse=', '), '\n')
cat('DEBUG: Number of records in census:', nrow(census), '\n')

# After importing mmwrdata and census, compare (state, year) pairs for coverage
mmwr_pairs <- tryCatch({
  unique(mmwrdata[, c("state", "year")])
}, error = function(e) {
  report_progress("WARNING", message=paste("Error creating MMWR pairs:", e$message))
  # Create fallback structure
  data.frame(state = unique(mmwrdata$state), year = max(as.numeric(mmwrdata$year), na.rm=TRUE), 
             stringsAsFactors = FALSE)
})

census_pairs <- tryCatch({
  unique(census[, c("state", "year")])
}, error = function(e) {
  report_progress("WARNING", message=paste("Error creating census pairs:", e$message))
  # Create fallback structure
  data.frame(state = unique(census$state), year = max(as.numeric(census$year), na.rm=TRUE), 
             stringsAsFactors = FALSE)
})

# Find (state, year) pairs in MMWR but not in census - with error handling
mmwr_not_in_census <- tryCatch({
  anti_join(mmwr_pairs, census_pairs, by = c("state", "year"))
}, error = function(e) {
  report_progress("WARNING", message=paste("Error finding MMWR records not in census:", e$message))
  data.frame(state = character(0), year = numeric(0), stringsAsFactors = FALSE)
})

# Find (state, year) pairs in census but not in MMWR - with error handling
census_not_in_mmwr <- tryCatch({
  anti_join(census_pairs, mmwr_pairs, by = c("state", "year"))
}, error = function(e) {
  report_progress("WARNING", message=paste("Error finding census records not in MMWR:", e$message))
  data.frame(state = character(0), year = numeric(0), stringsAsFactors = FALSE)
})

# Print summary to console
cat('PREPROCESS CHECK: (state, year) pairs in MMWR but missing in census:', nrow(mmwr_not_in_census), '\n')
if (nrow(mmwr_not_in_census) > 0) {
  print(mmwr_not_in_census)
}
cat('PREPROCESS CHECK: (state, year) pairs in census but missing in MMWR:', nrow(census_not_in_mmwr), '\n')
if (nrow(census_not_in_mmwr) > 0) {
  print(census_not_in_mmwr)
}

# Apply state filtering to census data if specified
if (!is.null(opts$states)) {
  states_to_analyze <- clean_list(opts$states)
  # Filter census data by states
  census <- census[toupper(census$state) %in% toupper(states_to_analyze), ]
  report_progress("DATA", message=paste("Filtered census data to", 
                                      length(unique(census$state)), "states"))
}

# Apply Salmonella serotype filtering if specified
if (!is.null(opts$salmonella_serotypes)) {
  serotypes_to_analyze <- clean_list(opts$salmonella_serotypes)
  report_progress("DATA", message=paste("Filtering for Salmonella serotypes:", 
                                      paste(serotypes_to_analyze, collapse=", ")))
  
  # Find the likely serotype column
  serotype_col <- NULL
  if ("serotypesummary" %in% names(mmwrdata)) {
    serotype_col <- "serotypesummary"
  } else if ("sero2" %in% names(mmwrdata)) {
    serotype_col <- "sero2"
  } else if ("sero1" %in% names(mmwrdata)) {
    serotype_col <- "sero1"
  }
  
  if (!is.null(serotype_col)) {
    # Filter Salmonella data by serotypes
    original_count <- nrow(mmwrdata)
    sal_rows <- mmwrdata$pathogen == "SALMONELLA" & mmwrdata[[serotype_col]] %in% serotypes_to_analyze
    other_path_rows <- mmwrdata$pathogen != "SALMONELLA"
    mmwrdata <- mmwrdata[sal_rows | other_path_rows, ]
    new_count <- nrow(mmwrdata)
    
    report_progress("DATA", message=paste("Filtered from", original_count, 
                                        "to", new_count, "records based on Salmonella serotype selection"))
  } else {
    report_progress("WARNING", message="Could not identify serotype column for filtering")
  }
}

# After Salmonella serotype filtering, add STEC serotype filtering
if (!is.null(opts$stec_serotypes)) {
  serotypes_to_analyze <- clean_list(opts$stec_serotypes)
  serotype_col <- NULL
  if ("serotypesummary" %in% names(mmwrdata)) {
    serotype_col <- "serotypesummary"
  } else if ("sero2" %in% names(mmwrdata)) {
    serotype_col <- "sero2"
  } else if ("sero1" %in% names(mmwrdata)) {
    serotype_col <- "sero1"
  }
  if (!is.null(serotype_col)) {
    stec_rows <- mmwrdata$pathogen == "STEC" & mmwrdata[[serotype_col]] %in% serotypes_to_analyze
    mmwrdata <- mmwrdata[stec_rows | mmwrdata$pathogen != "STEC", ]
    report_progress("DATA", message=paste("Filtering for STEC serotypes:",
                                          paste(serotypes_to_analyze, collapse=", ")))
  } else {
    report_progress("WARNING", message="Could not identify STEC serotype column for filtering")
  }
}

report_progress("DATA", message=paste("Processed", nrow(mmwrdata), "MMWR records"))

# After the filtering, re-run the check for census/mmwr pair comparison with error handling
mmwr_pairs <- tryCatch({
  unique(mmwrdata[, c("state", "year")])
}, error = function(e) {
  report_progress("WARNING", message=paste("Error creating MMWR pairs after filtering:", e$message))
  # Create fallback structure
  data.frame(state = unique(mmwrdata$state), year = max(as.numeric(mmwrdata$year), na.rm=TRUE), 
             stringsAsFactors = FALSE)
})

census_pairs <- tryCatch({
  unique(census[, c("state", "year")])
}, error = function(e) {
  report_progress("WARNING", message=paste("Error creating census pairs after filtering:", e$message))
  # Create fallback structure
  data.frame(state = unique(census$state), year = max(as.numeric(census$year), na.rm=TRUE), 
             stringsAsFactors = FALSE)
})

# Find (state, year) pairs in MMWR but not in census - with error handling
mmwr_not_in_census <- tryCatch({
  anti_join(mmwr_pairs, census_pairs, by = c("state", "year"))
}, error = function(e) {
  report_progress("WARNING", message=paste("Error finding MMWR records not in census after filtering:", e$message))
  data.frame(state = character(0), year = numeric(0), stringsAsFactors = FALSE)
})

# Find (state, year) pairs in census but not in MMWR - with error handling
census_not_in_mmwr <- tryCatch({
  anti_join(census_pairs, mmwr_pairs, by = c("state", "year"))
}, error = function(e) {
  report_progress("WARNING", message=paste("Error finding census records not in MMWR after filtering:", e$message))
  data.frame(state = character(0), year = numeric(0), stringsAsFactors = FALSE)
})

# Print summary to console
cat('PREPROCESS CHECK: (state, year) pairs in MMWR but missing in census:', nrow(mmwr_not_in_census), '\n')
if (nrow(mmwr_not_in_census) > 0) {
  print(mmwr_not_in_census)
}
cat('PREPROCESS CHECK: (state, year) pairs in census but missing in MMWR:', nrow(census_not_in_mmwr), '\n')
if (nrow(census_not_in_mmwr) > 0) {
  print(census_not_in_mmwr)
}

# Debug information for census and mmwrdata
cat('DEBUG: Unique pathogens in mmwrdata:', paste(unique(mmwrdata$pathogen), collapse=', '), '\n')
cat('DEBUG: Unique years in mmwrdata:', paste(unique(mmwrdata$year), collapse=', '), '\n')
cat('DEBUG: Unique states in mmwrdata:', paste(unique(mmwrdata$state), collapse=', '), '\n')
cat('DEBUG: Number of records in mmwrdata:', nrow(mmwrdata), '\n')
cat('DEBUG: Unique states in census:', paste(unique(census$state), collapse=', '), '\n')
cat('DEBUG: Unique years in census:', paste(unique(census$year), collapse=', '), '\n')
cat('DEBUG: Number of records in census:', nrow(census), '\n')