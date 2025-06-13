#!/usr/bin/env Rscript
#
# trendy.R - Main script for FoodNet Trends Bayesian modeling
#
# This script implements a Bayesian hierarchical model with splines to analyze
# foodborne illness surveillance data from the FoodNet program. It processes
# multiple pathogens, fits models, and generates incidence rate estimates.

#
# Usage:
#   ...
#
# Example:
#   ...
#
################################################################################

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
parser$add_argument("--states", type="character", default=NULL,
                    help="List of states to include (default: all states)")
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
                    help="Specific pathogen to analyze (if not processing all)")
parser$add_argument("--subgroup", type="character", default="combined",
                    help="Pathogen subgroup to analyze (e.g., O157, Enteritidis)")

# Preprocessing parameters
parser$add_argument("--preprocessed", type="logical", default=FALSE,
                    help="Use preprocessed CSV data (default: FALSE)")
parser$add_argument("--cleanFile", type="character", default=NULL,
                    help="Path to cleaned CSV file if preprocessed is TRUE")

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

# Configuration parameters
parser$add_argument("--catchment-config", type="character", default=NULL,
                    help="Path to CSV file with catchment definitions (optional)")

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
  travelLabel <- "All Cases"
} else if (!("YES" %in% travel) & ("UNKNOWN" %in% travel)) {
  travelLabel <- "Domestically-Acquired (UNK Travel Included)"
} else {
  travelLabel <- "Domestically-Acquired (UNK Travel Excluded)"
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

##############################################################
# Data Import and Preprocessing
##############################################################

# Import MMWR data
report_progress("DATA", message="Importing MMWR data")
tryCatch({
  # Check if we're using preprocessed data
  if (preprocessed) {
    # Validate that cleanFile is provided and exists
    if (is.null(cleanFile) || cleanFile == "") {
      stop("preprocessed=TRUE but no cleanFile path provided")
    }
    if (!file.exists(cleanFile)) {
      stop("Preprocessed file not found: ", cleanFile)
    }
    report_progress("DATA", message=paste("Using preprocessed data from:", cleanFile))
    # Read the preprocessed CSV file
    mmwrdata <- readr::read_csv(cleanFile, show_col_types = FALSE)
  } else {
    # This should not happen in the Nextflow pipeline context
    stop("This script is designed to work with preprocessed data. Please ensure preprocess.R has been run first.")
  }
  
  # Apply filters specific to this analysis
  mmwrdata <- mmwrdata %>%
    # Filter by detection method and travel status
    filter((cxcidt %in% cidt) & (travelint %in% travel)) %>%
    # Exclude invalid counties (if not already done in preprocessing)
    filter(!county %in% c("OUT OF STATE", "UNKNOWN", "99997"))
  
  # Apply state filter if specified
  if (!is.null(opts$states) && opts$states != "") {
    states_list <- CLEAN_LIST(opts$states)
    mmwrdata <- mmwrdata %>%
      filter(state %in% states_list)
    report_progress("DATA", message=paste("Filtered to states:", paste(states_list, collapse=", ")))
  }
  
  # Ensure required columns exist
  required_cols <- c("pathogen", "year", "state", "pathogentype")
  missing_cols <- required_cols[!required_cols %in% names(mmwrdata)]
  if (length(missing_cols) > 0) {
    stop("Required columns missing from preprocessed data: ", paste(missing_cols, collapse=", "))
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
  
  census <- as.data.frame(census)
  report_progress("DATA", message=paste("Processed census data with",
                                        length(unique(census$year)), "years and",
                                        length(unique(census$state)), "states"))
}, error = function(e) {
  stop("Error importing census data: ", e$message)
})

##############################################################
# Load Catchment Configuration
##############################################################

# Load catchment configuration if provided
catchment_config <- NULL
if (!is.null(opts$`catchment-config`)) {
  report_progress("CONFIG", message=paste("Loading catchment configuration from:", opts$`catchment-config`))
  tryCatch({
    catchment_config <- read_catchment_config(opts$`catchment-config`)
    report_progress("CONFIG", message=paste("Loaded catchment configuration with", nrow(catchment_config), "sites"))
  }, error = function(e) {
    stop("Error loading catchment configuration: ", e$message)
  })
} else {
  report_progress("CONFIG", message="Using default FoodNet catchment definitions")
}

##############################################################
# Pathogen Analysis
##############################################################

# Process pathogen data
report_progress("ANALYSIS", message="Processing pathogen data")
tryCatch({
  pathDf <- PATH_ANALYSIS(mmwrdata, census, catchment_config)%>%as.data.frame()
  report_progress("ANALYSIS", message=paste("Processed",
                                            length(unique(pathDf$pathogen)),
                                            "pathogens"))
  
  # Process Cyclospora and Salmonella if CIDT+ is included
  if("CIDT+" %in% cidt) {
    report_progress("ANALYSIS", message="Processing Cyclospora data")
    cyloDF <- CYCLOSPORA_ANALYSIS(mmwrdata, census, catchment_config)%>%as.data.frame()
    
    report_progress("ANALYSIS", message="Processing Salmonella data")
    salDF <- SALMONELLA_ANALYSIS(mmwrdata, census, catchment_config)%>%as.data.frame()
    
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
    # If a specific pathogen was requested, filter for it
    bact <- subset(bact, pathogen == opts$pathogen)
    
    # Further filter by subgroup if specified
    if (opts$subgroup != "combined") {
      # Handle different subgroup types
      if (opts$pathogen == "STEC" && opts$subgroup %in% c("O157", "nonO157")) {
        # For STEC, need to create the subgroups based on stec_class
        if ("stec_class" %in% names(bact)) {
          if (opts$subgroup == "O157") {
            bact <- subset(bact, stec_class == "STEC O157")
          } else if (opts$subgroup == "nonO157") {
            bact <- subset(bact, stec_class %in% c("STEC NONO157", "STEC O AG UNDET"))
          }
        } else {
          report_progress("WARNING", message="stec_class column not found - cannot filter by STEC subgroup")
        }
      } else if (opts$pathogen == "SALMONELLA") {
        # For Salmonella, filter by serotype
        if ("serotypesummary" %in% names(bact)) {
          bact <- subset(bact, serotypesummary == opts$subgroup)
        } else {
          report_progress("WARNING", message="serotypesummary column not found - cannot filter by serotype")
        }
      }
    }
    
    if (nrow(bact) == 0) {
      # No data found for the requested pathogen/subgroup
      pathogen_desc <- ifelse(opts$subgroup == "combined", 
                              opts$pathogen, 
                              paste(opts$pathogen, opts$subgroup, sep=":"))
      report_progress("ERROR", message=paste("No data found for:", pathogen_desc))
      
      # Create error file for this pathogen
      error_file <- paste0(outDir, "/", opts$pathogen, "_", opts$subgroup, "_error.txt")
      error_content <- c(
        paste("ERROR: No data found for:", pathogen_desc),
        paste("Date:", Sys.time()),
        paste("Project ID:", projID),
        "",
        "This pathogen/subgroup had no cases in the dataset after applying the following filters:",
        paste("- Pathogen:", opts$pathogen),
        paste("- Subgroup:", opts$subgroup),
        paste("- Travel types:", paste(travel, collapse=", ")),
        paste("- CIDT types:", paste(cidt, collapse=", ")),
        paste("- Time period: Check your input data file"),
        "",
        "Please verify:",
        "1. The pathogen name is spelled correctly",
        "2. The subgroup/serotype exists in your data",
        "3. The filters (travel, CIDT) are not excluding all cases"
      )
      writeLines(error_content, error_file)
      
      # Exit with error status to indicate failure
      stop(paste("No data found for:", pathogen_desc, "- see error file for details"))
    }
  } else {
    # No specific pathogen requested - analyze all pathogens in the data
    report_progress("ANALYSIS", message="No specific pathogen requested, analyzing all pathogens in dataset")
    
    # Check if any data exists
    if (nrow(bact) == 0) {
      stop("No data found after applying filters")
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
    
    # Construct output filename prefix including subgroup if specified
    output_prefix <- paste(pathogen_name, opts$subgroup, sep="_")
    
    # Save model
    saveFile <- paste0(outDir, "/", output_prefix, "_brm.Rds")
    saveRDS(proposed, saveFile)
    report_progress("MODEL", message=paste("Saved model to", saveFile))
    
    # Save model summary
    summaryFile <- paste0(outDir, "/", output_prefix, "_summary.txt")
    sink(summaryFile)
    print(summary(proposed))
    sink()
    report_progress("MODEL", message=paste("Saved model summary to", summaryFile))
    
    # Draw untransformed (link-level) predictions
    report_progress("POST-PROCESSING", message=paste("Generating predictions for", pathogen_name))
    posteriorLinpred <- LINPREAD_DRAW_FN(
      data = (proposed$data %>% group_by(state)),
      model = proposed
    )
    
    # site-level estimates
    report_progress("POST-PROCESSING", message="Calculating catchment-level draws")
    site <- LINPRED_TO_SITEIR(posteriorLinpred)
    # Add metadata
    site$pathogen <- pathogen_name
    site$travel <- travelLabel
    site$culture <- culture
    
    # Save site-level estimates
    siteir_file <- paste0(outDir, "/", output_prefix, "_IRSite.csv")
    write.csv(site, siteir_file, row.names = FALSE)
    report_progress("OUTPUT", message=paste("Saved site incidence rate estimates to", siteir_file))
    
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
    ir_file <- paste0(outDir, "/", output_prefix, "_IRCatch.csv")
    write.csv(catchir.linpred, ir_file, row.names = FALSE)
    report_progress("OUTPUT", message=paste("Saved incidence rate estimates to", ir_file))
    
    # Calculate relative risks and percent changes for different comparison periods
    report_progress("ANALYSIS", message="Calculating relative risks and percent changes")
    
    # Calculate for 2016-2018 (the Healthy People 2030 baseline period)
    hp30<-IR_COMP_CATCH(catch, catchir.linpred, 2016, 2018,
                  paste0(outDir, "/", output_prefix, "_EstIRRCatch_2016_2018.csv"))
    
    # Calculate for COVID-19
    IR_COMP_CATCH(catch, catchir.linpred, 2020, 2021,
                  paste0(outDir, "/", output_prefix, "_EstIRRCatch_2020_2022.csv"))
    
    # Calculate for earliest years where the FoodNet catchment were stable
    IR_COMP_CATCH(catch, catchir.linpred, 2004, 2006,
                  paste0(outDir, "/", output_prefix, "_EstIRRCatch_2004_2006.csv"))
    
    # Calculate for 2006-2008 baseline (the Healthy People 2020 baseline)
    IR_COMP_CATCH(catch, catchir.linpred, 2006, 2008,
                  paste0(outDir, "/", output_prefix, "_EstIRRCatch_2006_2008.csv"))
    
    # Create visualizations if enabled
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      tryCatch({
        # Check if plotting functions exist before calling them
        if (exists("PLOT_SITE_TRENDS", mode = "function")) {
          # Site-specific trends plot
          site_plot <- PLOT_SITE_TRENDS(site, pathogen_name, outDir)
          ggsave(paste0(outDir, "/", output_prefix, "_site_trends.png"), site_plot, 
                 width = 10, height = 8, dpi = 300)
        }
        
        if (exists("PLOT_OVERALL_TREND", mode = "function")) {
          # Overall trend plot
          overall_plot <- PLOT_OVERALL_TREND(catchir.linpred, pathogen_name, outDir)
          ggsave(paste0(outDir, "/", output_prefix, "_overall_trend.png"), overall_plot,
                 width = 10, height = 6, dpi = 300)
        }
      }, error = function(e) {
        report_progress("WARNING", message=paste("Visualization skipped:", e$message))
      })
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
