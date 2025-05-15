#!/usr/bin/env Rscript
#
# trendy.R - Main script for FoodNet Trends Bayesian modeling
#
# This script implements a Bayesian hierarchical model with splines to analyze
# foodborne illness surveillance data from the FoodNet program. It processes
# multiple pathogens, fits models, and generates incidence rate estimates.
#

# Suppress warnings during package loading
suppressPackageStartupMessages(library("argparse"))
options(warn = 1)  # Show warnings as they occur

# Determine script directory and source helper functions
script_path <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_path[grep("--file=", script_path)])
script_dir <- dirname(script_path)

# Source helper functions with robust path handling
tryCatch({
  cat("Attempting to source functions.R from script directory:", script_dir, "\n")
  source(file.path(script_dir, "functions.R"))
}, error = function(e) {
  # Try to find functions.R in the parent directory of the script
  cat("Trying parent directory...\n")
  tryCatch({
    parent_dir <- dirname(script_dir)
    source(file.path(parent_dir, "bin", "functions.R"))
  }, error = function(e2) {
    # Try the current working directory as a last resort
    cat("Trying current working directory...\n")
    tryCatch({
      source("functions.R")
    }, error = function(e3) {
      # If all attempts fail, provide diagnostic information and stop
      cat("Failed to locate functions.R. Script directory:", script_dir, "\n")
      cat("Current working directory:", getwd(), "\n")
      cat("Files in script directory:", paste(list.files(script_dir), collapse=", "), "\n")
      cat("Files in current directory:", paste(list.files("."), collapse=", "), "\n")
      stop("Error loading functions.R: ", e3$message)
    })
  })
})

# Load required packages function
LOAD_PACKAGES <- function(packages) {
  for(pkg in packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Required package", pkg, "is not installed"))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

##############################################################
# Setup and argument parsing
##############################################################

# Create parser object with comprehensive options
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
parser$add_argument("--discovery_data", type="character", default=NULL,
                    help="Path to discovery data JSON file from previous run")

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

# Parse arguments with error handling
tryCatch({
  opts <- parser$parse_args()
}, error = function(e) {
  cat("Error parsing command line arguments:", e$message, "\n")
  cat("Run with --help for usage information\n")
  quit(status = 1)
})

##############################################################
# Initialize variables based on arguments
##############################################################

# Report progress
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
  travel <- CLEAN_LIST(opts$travel)
  cidt <- CLEAN_LIST(opts$cidt)

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
  discovery_data <- opts$discovery_data

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
  travel <- CLEAN_LIST("NO,UNKNOWN,YES")
  cidt <- CLEAN_LIST("CIDT+,CX+,PARASITIC")

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
  discovery_data <- NULL
}

# Validate required parameters
validate_params <- function() {
  errors <- c()

  # Check required file parameters
  if (is.null(mmwrFile) || mmwrFile == "")
    errors <- c(errors, "Missing required parameter: mmwrFile")
  if (is.null(censusFileB) || censusFileB == "")
    errors <- c(errors, "Missing required parameter: censusFileB")
  if (is.null(censusFileP) || censusFileP == "")
    errors <- c(errors, "Missing required parameter: censusFileP")

  # Check file existence
  if (length(errors) == 0) {
    if (!file.exists(mmwrFile))
      errors <- c(errors, paste("MMWR file does not exist:", mmwrFile))
    if (!file.exists(censusFileB))
      errors <- c(errors, paste("Census bacterial file does not exist:", censusFileB))
    if (!file.exists(censusFileP))
      errors <- c(errors, paste("Census parasitic file does not exist:", censusFileP))
  }

  # Check preprocessed file if specified
  if (preprocessed && !is.null(cleanFile)) {
    if (!file.exists(cleanFile))
      errors <- c(errors, paste("Clean file does not exist:", cleanFile))
  }

  # Check discovery data if specified
  if (!is.null(discovery_data) && discovery_data != "") {
    if (!file.exists(discovery_data))
      errors <- c(errors, paste("Discovery data file does not exist:", discovery_data))
  }

  # Check project ID
  if (is.null(projID) || projID == "") {
    projID <<- format(Sys.time(), "%Y%m%d%H%M")
    report_progress("SETUP", message=paste("No projID provided, using timestamp:", projID))
  }

  # Return errors if any
  if (length(errors) > 0) {
    for (err in errors) {
      report_progress("ERROR", message=err)
    }
    stop(paste(errors, collapse="\n"))
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
if (!is.null(discovery_data) && discovery_data != "" && file.exists(discovery_data)) {
  report_progress("SETUP", message=paste("Loading discovery data from:", discovery_data))
  
  tryCatch({
    # Load packages needed for JSON
    suppressPackageStartupMessages(library(jsonlite))
    
    # Read the discovery data
    discovery <- jsonlite::read_json(discovery_data)
    
    # Log what was found
    report_progress("SETUP", message=paste("Found", length(discovery$pathogens), "pathogens and", 
                                         length(discovery$states), "states in discovery data"))
    
    # Check if we need to override pathogen and state lists
    if (is.null(opts$pathogen)) {
      # If no pathogen was specified, use the first two from discovery
      if (length(discovery$pathogens) >= 2) {
        opts$pathogen <- paste(discovery$pathogens[1:2], collapse=",")
        report_progress("SETUP", message=paste("No pathogens specified, using first two from discovery:", opts$pathogen))
      }
    }
    
    if (is.null(opts$states)) {
      # If no states were specified, use all from discovery
      opts$states <- paste(discovery$states, collapse=",")
      report_progress("SETUP", message=paste("No states specified, using all from discovery"))
    }
    
  }, error = function(e) {
    report_progress("WARNING", message=paste("Error loading discovery data:", e$message))
    report_progress("WARNING", message="Continuing with command-line parameters only")
  })
}

# Load required packages
report_progress("SETUP", message="Loading required packages")
pkgs <- c('haven', 'gtools', 'brms', 'ggplot2', 'tidybayes', 'HDInterval', 'tidyverse')
tryCatch({
  LOAD_PACKAGES(pkgs)
}, error = function(e) {
  stop("Failed to load required packages: ", e$message)
})

##############################################################
# Set up analysis parameters
##############################################################

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
  states_to_analyze <- CLEAN_LIST(opts$states)
  report_progress("ANALYSIS DETAILS", message=paste("Filtering states:", paste(states_to_analyze, collapse=",")))
}

# If Salmonella serotype filtering is specified, report it
if (!is.null(opts$salmonella_serotypes)) {
  serotypes_to_analyze <- CLEAN_LIST(opts$salmonella_serotypes)
  report_progress("ANALYSIS DETAILS", message=paste("Filtering Salmonella serotypes:", paste(serotypes_to_analyze, collapse=",")))
}

##############################################################
# Data Import and Preprocessing
##############################################################

# Import MMWR data
report_progress("DATA", message="Importing MMWR data")
tryCatch({
  if (preprocessed && !is.null(cleanFile) && file.exists(cleanFile)) {
    report_progress("DATA", message=paste("Using preprocessed data from:", cleanFile))
    # Read the preprocessed CSV file
    mmwrdata <- readr::read_csv(cleanFile, show_col_types = FALSE)
    
    # Make column names consistent - ensure key columns are lowercase
    names(mmwrdata) <- gsub("^Pathogen$", "pathogen", names(mmwrdata), ignore.case = TRUE)
    names(mmwrdata) <- gsub("^State$", "state", names(mmwrdata), ignore.case = TRUE) 
    names(mmwrdata) <- gsub("^Year$", "year", names(mmwrdata), ignore.case = TRUE)
    
    report_progress("DATA", message="Standardized column names from preprocessed file")
  } else {
    # Read raw SAS data
    report_progress("DATA", message=paste("Reading raw SAS data from:", mmwrFile))
    mmwrdata <- haven::read_sas(mmwrFile) %>% as.data.frame()
    
    # Convert all column names to lowercase for consistency
    names(mmwrdata) <- tolower(names(mmwrdata))
    
    # Additional data cleaning as needed
    report_progress("DATA", message="Standardized column names from raw SAS file")
  }

  # Define standard list of pathogens for filtering
  pathogens <- c("CAMPYLOBACTER", "CYCLOSPORA", "SALMONELLA", "SHIGELLA", "STEC", "VIBRIO", "YERSINIA")
  report_progress("DATA", message="Defined standard pathogen list")

  # Ensure pathogen column has consistent casing for filtering
  if ("pathogen" %in% names(mmwrdata)) {
    # Standardize pathogen names to uppercase
    mmwrdata$pathogen <- toupper(mmwrdata$pathogen)
    report_progress("DATA", message="Standardized pathogen column for consistent filtering")
  }

  # Ensure required columns exist
  required_cols <- c("pathogen", "year", "state")
  missing_cols <- required_cols[!required_cols %in% names(mmwrdata)]
  if (length(missing_cols) > 0) {
    stop("Required columns missing from MMWR data: ", paste(missing_cols, collapse=", "))
  }
  
  # Ensure pathogentype column exists
  if (!"pathogentype" %in% names(mmwrdata)) {
    mmwrdata$pathogentype <- ifelse(mmwrdata$pathogen %in% c("CRYPTOSPORIDIUM", "CYCLOSPORA"), 
                                   "Parasitic", "Bacterial")
  }

  # Apply state filtering if specified
  if (!is.null(opts$states)) {
    states_to_analyze <- CLEAN_LIST(opts$states)
    report_progress("DATA", message=paste("Filtering for states:", paste(states_to_analyze, collapse=", ")))
    
    # Filter data by states
    original_count <- nrow(mmwrdata)
    mmwrdata <- mmwrdata[toupper(mmwrdata$state) %in% toupper(states_to_analyze), ]
    new_count <- nrow(mmwrdata)
    
    report_progress("DATA", message=paste("Filtered from", original_count, 
                                        "to", new_count, "records based on state selection"))
  }

  # Apply Salmonella serotype filtering if specified
  if (!is.null(opts$salmonella_serotypes)) {
    serotypes_to_analyze <- CLEAN_LIST(opts$salmonella_serotypes)
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
      # This creates a mask for the entire dataset where only the specified serotypes are TRUE
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

  report_progress("DATA", message=paste("Processed", nrow(mmwrdata), "MMWR records"))
}, error = function(e) {
  stop("Error importing MMWR data: ", e$message)
})

# Import census data
report_progress("DATA", message="Importing census data")
tryCatch({
  census <- haven::read_sas(censusFileB) %>%
    setNames(tolower(names(.))) %>%
    group_by(year, state) %>%
    dplyr::summarize(population = sum(population, na.rm=TRUE)) %>%
    mutate(pathogentype = "Bacterial") %>%
    bind_rows(
      haven::read_sas(censusFileP) %>%
        setNames(tolower(names(.))) %>%
        group_by(year, state) %>%
        dplyr::summarize(population = sum(population, na.rm=TRUE)) %>%
        mutate(pathogentype = "Parasitic")
    ) %>%
    ungroup()

  # Apply state filtering to census data if specified
  if (!is.null(opts$states)) {
    states_to_analyze <- CLEAN_LIST(opts$states)
    # Filter census data by states
    census <- census[toupper(census$state) %in% toupper(states_to_analyze), ]
    report_progress("DATA", message=paste("Filtered census data to", 
                                        length(unique(census$state)), "states"))
  }

  census <- as.data.frame(census)
  report_progress("DATA", message=paste("Processed census data with",
                                       length(unique(census$year)), "years and",
                                       length(unique(census$state)), "states"))
}, error = function(e) {
  stop("Error importing census data: ", e$message)
})

##############################################################
# Pathogen Analysis
##############################################################

# Process pathogen data
report_progress("ANALYSIS", message="Processing pathogen data")
tryCatch({
  pathDf <- PATH_ANALYSIS(mmwrdata, census)
  report_progress("ANALYSIS", message=paste("Processed",
                                          length(unique(pathDf$pathogen)),
                                          "pathogens"))

  # Process Cyclospora and Salmonella if CIDT+ is included
  if("CIDT+" %in% cidt) {
    report_progress("ANALYSIS", message="Processing Cyclospora data")
    cyloDF <- CYCLOSPORA_ANALYSIS(mmwrdata, census)

    report_progress("ANALYSIS", message="Processing Salmonella data")
    salDF <- SALMONELLA_ANALYSIS(mmwrdata, census)

    # Combine all pathogen data
    bact <- gtools::smartbind(pathDf, cyloDF) %>%
      gtools::smartbind(salDF)
  } else {
    bact <- pathDf
  }

  # Post-processing
  report_progress("ANALYSIS", message="Post-processing pathogen data")

  # Clean up memory
  remove(mmwrdata)

  # Filter and prepare data for modeling
  if (!is.null(opts$pathogen)) {
    # If specific pathogens were requested, parse and filter for them
    pathogens_to_analyze <- unlist(strsplit(opts$pathogen, ","))
    report_progress("ANALYSIS", message=paste("Filtering for requested pathogens:",
                                              paste(pathogens_to_analyze, collapse=", ")))

    # Ensure consistent case for pathogen filtering
    bact$pathogen <- toupper(bact$pathogen)
    pathogens_to_analyze <- toupper(pathogens_to_analyze)

    # Filter for requested pathogens
    bact <- subset(bact, pathogen %in% pathogens_to_analyze)

    if (nrow(bact) == 0) {
      # Instead of stopping, create a minimal dataset for the pathogen
      # This will allow the pipeline to continue but produce empty results
      report_progress("WARNING", message=paste("No data found for requested pathogens - Creating minimal dataset"))

      # Create a minimal dataset with the requested pathogen for all sites
      states <- unique(census$state)
      years <- unique(census$year)

      # Create a minimal dataset for each requested pathogen
      minimal_data_list <- list()

      for (pathogen_name in pathogens_to_analyze) {
        minimal_data <- expand.grid(
          year = years,
          state = states,
          pathogen = pathogen_name,
          stringsAsFactors = FALSE
        )

        # Add required columns
        minimal_data$count <- 0

        # Determine pathogen type
        if (pathogen_name %in% c("CRYPTOSPORIDIUM", "CYCLOSPORA")) {
          pathogen_type <- "Parasitic"
        } else {
          pathogen_type <- "Bacterial"
        }

        minimal_data$pathogentype <- pathogen_type
        minimal_data <- left_join(minimal_data,
                                census %>% filter(pathogentype == pathogen_type),
                                by = c("year", "state"))

        # Remove any NA rows that might have been created in the join
        minimal_data <- minimal_data[!is.na(minimal_data$population), ]

        minimal_data_list[[pathogen_name]] <- minimal_data
      }

      # Combine all minimal datasets
      bact <- do.call(rbind, minimal_data_list)
    }
  } else {
    # Otherwise use the default filtering from the original code
    bact <- subset(bact, pathogen == "CAMPYLOBACTER" | pathogen == "CYCLOSPORA")

    # Check if any data exists for the default pathogens
    if (nrow(bact) == 0) {
      stop("No data found for default pathogens (CAMPYLOBACTER or CYCLOSPORA)")
    }
  }

  # Prepare year variables and split by pathogen
  bact$yearn <- as.numeric(as.character(bact$year))
  bact$year <- as.factor(bact$year)
  bact_list <- split(bact, bact$pathogen)
  target_pathogens <- names(bact_list)

  report_progress("ANALYSIS", message=paste("Prepared data for modeling",
                                          length(target_pathogens),
                                          "pathogens:",
                                          paste(target_pathogens, collapse=", ")))
}, error = function(e) {
  stop("Error in pathogen analysis: ", e$message)
})

##############################################################
# Model Fitting
##############################################################

# Process each pathogen
for (pathogen_name in target_pathogens) {
  report_progress("MODEL", message=paste("Fitting model for", pathogen_name))

  # Get data for current pathogen
  current_data <- bact_list[[pathogen_name]]

  # Fit Bayesian model
  tryCatch({
    # Fit model with parameters from command line
    proposed <- PROPOSED_BM(
      current_data,
      cores = modelcores,
      chains = chains,
      iterations = iterations,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      seed = seed
    )

    # Save model
    saveFile <- paste0(outDir, "/", pathogen_name, "_brm.Rds")
    saveRDS(proposed, saveFile)
    report_progress("MODEL", message=paste("Saved model to", saveFile))

    # Save model summary
    summaryFile <- paste0(outDir, "/", pathogen_name, "_summary.txt")
    sink(summaryFile)
    print(summary(proposed))
    sink()
    report_progress("MODEL", message=paste("Saved model summary to", summaryFile))

    # Draw untransformed (link-level) predictions
    report_progress("POST-PROCESSING", message=paste("Generating predictions for", pathogen_name))
    posteriorLinpred <- LINPREAD_DRAW_FN(
      data = (current_data %>% group_by(state)),
      model = proposed
    )

    # site-level estimates
    report_progress("POST-PROCESSING", message="Calculating catchment-level draws")
    site <- LINPRED_TO_SITEIR(posteriorLinpred)

    # Catchment-level draws
    report_progress("POST-PROCESSING", message="Calculating catchment-level draws")
    catch <- CATCHMENT(posteriorLinpred)

    # Catchment-level estimates
    report_progress("POST-PROCESSING", message="Calculating catchment-level estimates")
    catchir.linpred <- LINPRED_TO_CATCHIR(catch)

    # Add metadata
    catchir.linpred$pathogen <- pathogen_name
    catchir.linpred$travel <- travelLabel
    catchir.linpred$culture <- culture

    # Save estimates
    ir_file <- paste0(outDir, "/", pathogen_name, "_IRCatch.csv")
    write.csv(catchir.linpred, ir_file, row.names = FALSE)
    report_progress("OUTPUT", message=paste("Saved incidence rate estimates to", ir_file))

    # Calculate relative risks and percent changes for different comparison periods
    report_progress("ANALYSIS", message="Calculating relative risks and percent changes")

    # Calculate for 2016-2018 (the Healthy People 2030 baseline period)
    IR_COMP_CATCH(catch, catchir.linpred, 2016, 2018,
            paste0(outDir, "/", pathogen_name, "_EstIRRCatch_2016_2018.csv"))

    # Calculate for COVID-19
    IR_COMP_CATCH(catch, catchir.linpred, 2020, 2021,
            paste0(outDir, "/", pathogen_name, "_EstIRRCatch_2020_2022.csv"))

    # Calculate for earliest years where the FoodNet catchment were stable
    IR_COMP_CATCH(catch, catchir.linpred, 2004, 2006,
            paste0(outDir, "/", pathogen_name, "_EstIRRCatch_2004_2006.csv"))

    # Calculate for 2006-2008 baseline (the Healthy People 2020 baseline)
    IR_COMP_CATCH(catch, catchir.linpred, 2006, 2008,
            paste0(outDir, "/", pathogen_name, "_EstIRRCatch_2006_2008.csv"))

    # Create visualizations if enabled
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Site-specific trends plot
      site_plot <- PLOT_SITE_TRENDS(catchir.linpred, pathogen_name, outDir)

      # Overall trend plot
      overall_plot <- PLOT_OVERALL_TREND(catchir.linpred, pathogen_name, outDir)

      # Combined visualization
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        PLOT_COMBINED(site_plot, overall_plot, pathogen_name, outDir)
      }
    }

    report_progress("COMPLETE", message=paste("Completed analysis for", pathogen_name))
  }, error = function(e) {
    report_progress("ERROR", message=paste("Error in model fitting for", pathogen_name, ":", e$message))
    # Create error file with details
    error_file <- paste0(outDir, "/", pathogen_name, "_error.txt")
    sink(error_file)
    cat(paste("Error processing", pathogen_name, "at", Sys.time(), "\n"))
    cat(paste("Error message:", e$message, "\n"))
    cat("Traceback:\n")
    cat(paste(capture.output(traceback()), collapse = "\n"))
    sink()

    # Continue with next pathogen rather than stopping the entire pipeline
    next
  })
}

report_progress("PIPELINE", message="Analysis complete for all pathogens")
report_progress("PIPELINE", message=paste("Results saved to", outDir))

# Print session info for reproducibility
report_progress("SESSION", message="Session information:")
print(sessionInfo())
