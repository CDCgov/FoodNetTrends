#!/usr/bin/env Rscript
# =========================================================================
# FoodNet Trends - Data Preprocessing and Cleaning
# =========================================================================
#
# Purpose:
#   Clean and standardize raw MMWR data for analysis by the trendy.R script.
#   Generate metadata for downstream discovery and filtering.
#
# The script performs the following operations:
#   1. Read raw SAS data file
#   2. Standardize column names and formats
#   3. Clean and recode serotype values
#   4. Standardize county names
#   5. Generate optional metadata JSON about the dataset contents
#   6. Write processed data to CSV
#
# Usage:
#   Rscript calcIR.R --mmwrFile <path_to_raw_SAS_file> --outputFile <path_to_output_csv> [--generate_metadata true/false]
#
# Output:
#   - Cleaned CSV file with standardized format
#   - Optional JSON metadata file with dataset summary
#
# Last updated: 2025-05-18
# =========================================================================

# Attempt to load functions from library
tryCatch({
  source("functions.R")
  cat("Successfully sourced functions.R\n")
}, error = function(e) {
  # If functions.R isn't available, try to find it in the script directory
  script_path <- commandArgs(trailingOnly = FALSE)
  script_path <- script_path[grep("--file=", script_path)]
  
  if (length(script_path) > 0) {
    script_path <- substring(script_path, 8)
    script_dir <- dirname(script_path)
    
    tryCatch({
      source(file.path(script_dir, "functions.R"))
      cat("Successfully sourced functions.R from script directory\n")
    }, error = function(e) {
      cat("Warning: Could not load functions.R from either location\n")
      cat("Script directory:", script_dir, "\n")
      cat("Current directory:", getwd(), "\n")
      cat("Directory contents:", paste(list.files("."), collapse=", "), "\n")
      
      # Continue without functions.R as it's not strictly required for this script
      cat("Continuing without functions.R - safe_write function will not be available\n")
    })
  }
})

# Load required packages
suppressPackageStartupMessages({
  library("argparse")
  library("dplyr")
  library("haven")
  library("jsonlite")
})

# Setup argument parser
parser <- ArgumentParser(description="Clean and standardize FoodNet MMWR data")
parser$add_argument("--mmwrFile", type = "character", help = "Path to the raw MMWR SAS file", required = TRUE)
parser$add_argument("--outputFile", type = "character", help = "Path to save the cleaned CSV file", required = TRUE)
parser$add_argument("--generate_metadata", type = "logical", default = FALSE, 
                   help = "Whether to generate a metadata JSON file (default: FALSE)")
parser$add_argument("--censusFileB", type = "character", help = "Path to census file for bacterial pathogens", required = TRUE)
parser$add_argument("--censusFileP", type = "character", help = "Path to census file for parasitic pathogens", required = TRUE)
args <- parser$parse_args()

#' Define helper functions if not available from functions.R
if(!exists("get_output_filename")) {
  #' Generate Standardized Filename
  #'
  #' @param base_name Base name for the file (e.g., dataset name)
  #' @param file_type Type of file (e.g., "metadata")
  #' @param extension File extension without dot (e.g., "json")
  #' @return A standardized filename string
  get_output_filename <- function(base_name, file_type, extension) {
    # Build filename with consistent pattern
    filename <- paste0(base_name, "_", file_type, ".", extension)
    return(filename)
  }
}

#' Define a safe_write function if it's not available from functions.R
#' 
#' @param data Data frame to write
#' @param file_path Full path to the output file (.csv or .Rds)
#' @return None
if(!exists("safe_write")) {
  safe_write <- function(data, file_path) {
    tryCatch({
      # Create directory if it doesn't exist
      dir_path <- dirname(file_path)
      if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Write data based on file extension
      if (endsWith(file_path, ".csv")) {
        if (file.exists(file_path)) {
          write.table(data, file = file_path, append = TRUE, quote = TRUE, sep = ",",
                      col.names = FALSE, row.names = FALSE)
        } else {
          write.table(data, file = file_path, append = FALSE, quote = TRUE, sep = ",",
                      col.names = TRUE, row.names = FALSE)
        }
      } else if (endsWith(file_path, ".Rds")) {
        saveRDS(data, file = file_path)
      }
    }, error = function(e) {
      message("Error writing file: ", e$message)
    })
  }
}

#' Load and validate raw MMWR data
#' 
#' @param file_path Path to the SAS file
#' @return Data frame with raw MMWR data
#' @throws Error if file cannot be read
load_mmwr_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("MMWR file does not exist: ", file_path)
  }
  
  tryCatch({
    cat("Loading raw MMWR data from:", file_path, "\n")
    data <- haven::read_sas(file_path) %>% as.data.frame()
    cat("Successfully loaded data with", nrow(data), "records and", ncol(data), "columns\n")
    return(data)
  }, error = function(e) {
    stop("Failed to read MMWR file: ", e$message)
  })
}

#' Clean and standardize MMWR data
#' 
#' @param raw_data Data frame with raw MMWR data
#' @return Cleaned data frame
clean_mmwr_data <- function(raw_data) {
  # Convert all column names to lowercase for consistency
  cat("Standardizing column names to lowercase...\n")
  data <- raw_data %>% rename_all(tolower)
  
  # Recode serotype values
  cat("Recoding serotype values...\n")
  seroList <- c("NOT SPECIATED", "UNKNOWN", "PARTIAL SERO", "NOT SERO", "")
  data$sero2 <- ifelse(data$sero1 %in% seroList, "Missing", data$sero1)
  data$sero2 <- ifelse(grepl("UNDET", data$sero2), "Missing", data$sero2)
  data$serotypesummary <- data$sero2
  
  # Standardize county names
  cat("Standardizing county names...\n")
  data <- data %>%
    mutate(
      county = if_else(county %in% c("ST. MARYS'S", "ST. MARYS"), "ST. MARY'S", county),
      county = if_else(county == "PRINCE GEORGES", "PRINCE GEORGE'S", county),
      county = if_else(county == "QUEEN ANNES", "QUEEN ANNE'S", county),
      county = if_else(county == "DE BACA", "DEBACA", county)
    )
  
  # Ensure pathogen column is uppercase for consistency
  data$pathogen <- toupper(data$pathogen)
  
  # Create pathogentype column if not present
  if(!"pathogentype" %in% names(data)) {
    cat("Creating derived pathogentype column...\n")
    data <- data %>%
      mutate(pathogentype = ifelse(pathogen %in% c("CRYPTOSPORIDIUM", "CYCLOSPORA"), 
                                  "Parasitic", "Bacterial"))
  }
  
  return(data)
}

#' Generate metadata from cleaned data
#' 
#' @param data Cleaned MMWR data
#' @param source_file Original source file path
#' @return List with metadata information
generate_metadata <- function(data, source_file) {
  cat("Generating metadata from cleaned data...\n")
  
  # Extract key information for metadata
  metadata <- list(
    pathogens = sort(unique(data$pathogen)),
    states = sort(unique(data$state)),
    years = sort(unique(as.numeric(as.character(data$year)))),
    counties = sort(unique(data$county)),
    generated_timestamp = as.character(Sys.time()),
    source_file = source_file,
    record_count = nrow(data)
  )
  
  # Add counts for basic statistics
  metadata$counts <- list(
    total_records = nrow(data),
    pathogen_counts = as.list(table(data$pathogen)),
    state_counts = as.list(table(data$state))
  )
  
  # Process Salmonella serotypes if available
  if (any(data$pathogen == "SALMONELLA") && "serotypesummary" %in% names(data)) {
    cat("Processing Salmonella serotype information...\n")
    sal_data <- data[data$pathogen == "SALMONELLA", ]
    serotype_counts <- as.data.frame(table(sal_data$serotypesummary))
    serotype_counts <- serotype_counts[order(serotype_counts$Freq, decreasing=TRUE),]
    
    # Store in metadata
    metadata$salmonella_serotypes <- as.list(serotype_counts$Freq)
    names(metadata$salmonella_serotypes) <- serotype_counts$Var1
    
    # Also store a flat list of names
    metadata$salmonella_serotype_names <- as.character(serotype_counts$Var1)
    
    cat("Found", length(metadata$salmonella_serotype_names), "Salmonella serotypes\n")
  } else {
    cat("No Salmonella serotype information available\n")
  }
  
  return(metadata)
}

# Main execution
main <- function() {
  # Load data
  mmwrdata <- load_mmwr_data(args$mmwrFile)
  
  # Clean data
  cleaned_data <- clean_mmwr_data(mmwrdata)

  # Load census files for coverage check
  cat("Loading census files for coverage check...\n")
  census_b <- haven::read_sas(args$censusFileB) %>% dplyr::mutate(state = toupper(as.character(state)), year = as.numeric(as.character(year)))
  census_p <- haven::read_sas(args$censusFileP) %>% dplyr::mutate(state = toupper(as.character(state)), year = as.numeric(as.character(year)))
  census <- dplyr::bind_rows(
    census_b %>% dplyr::mutate(pathogentype = "Bacterial"),
    census_p %>% dplyr::mutate(pathogentype = "Parasitic")
  )
  census_pairs <- unique(census[, c("state", "year")])
  cleaned_data$state <- toupper(as.character(cleaned_data$state))
  cleaned_data$year <- as.numeric(as.character(cleaned_data$year))
  mmwr_pairs <- unique(cleaned_data[, c("state", "year")])

  # Find (state, year) pairs in cleaned data but not in census
  mmwr_not_in_census <- dplyr::anti_join(mmwr_pairs, census_pairs, by = c("state", "year"))

  # Write report to preprocessed directory
  output_dir <- dirname(args$outputFile)
  coverage_report_file <- file.path(output_dir, "state_year_coverage_report.txt")
  cat("State-Year Coverage Report\n", file=coverage_report_file)
  cat("========================\n", file=coverage_report_file, append=TRUE)
  cat("(state, year) pairs in cleaned data but missing in census (n=", nrow(mmwr_not_in_census), "):\n", sep="", file=coverage_report_file, append=TRUE)
  if (nrow(mmwr_not_in_census) > 0) {
    write.table(mmwr_not_in_census, file=coverage_report_file, append=TRUE, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  } else {
    cat("(None)\n", file=coverage_report_file, append=TRUE)
  }

  # Print warning and drop records if needed
  if (nrow(mmwr_not_in_census) > 0) {
    cat("WARNING: Some (state, year) pairs in the cleaned data are missing from the census files.\n")
    cat("These records will be removed before modeling.\n")
    cat("See state_year_coverage_report.txt in the preprocessed data directory for details.\n")
    # Drop records with missing census coverage
    cleaned_data <- dplyr::semi_join(cleaned_data, census_pairs, by = c("state", "year"))
    cat("Dropped", nrow(mmwr_not_in_census), "records with missing census coverage.\n")
  }

  # Extract base name without extension for consistent naming
  output_base <- tools::file_path_sans_ext(basename(args$outputFile))
  output_dir <- dirname(args$outputFile)
  
  # Generate metadata if requested
  if(args$generate_metadata) {
    # Create standardized metadata filename
    metadata_filename <- get_output_filename(output_base, "metadata", "json")
    metadata_file <- file.path(output_dir, metadata_filename)
    
    # Generate metadata
    metadata <- generate_metadata(cleaned_data, args$mmwrFile)
    
    # Write metadata JSON
    cat("Writing metadata to:", metadata_file, "\n")
    write_json(metadata, metadata_file, pretty = TRUE)
  }
  
  # Write cleaned data to CSV using safe_write
  cat("Writing cleaned data to:", args$outputFile, "\n")
  safe_write(cleaned_data, args$outputFile)
  
  cat("Data cleaning complete. Cleaned data saved to:", args$outputFile, "\n")
  if(args$generate_metadata) {
    if(exists("metadata_file")) {
      cat("Metadata saved to:", metadata_file, "\n")
    }
  }
}

# Run the main function
main()
