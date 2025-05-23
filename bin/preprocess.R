#!/usr/bin/env Rscript
# =========================================================================
# FoodNetTrends v1.0 - Data Preprocessing and Cleaning
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
#   Rscript preprocess.R --mmwrFile <path_to_raw_SAS_file> --outputFile <path_to_output_csv> [--generate_metadata true/false]
#
# Output:
#   - Cleaned CSV file with standardized format
#   - Optional JSON metadata file with dataset summary
#
# Last updated: 2025-05-22
# =========================================================================

# Attempt to load functions from library
tryCatch({
  source("functions.R")
  cat("Successfully sourced functions.R\n")
}, error = function(e) {
  # Attempt to locate functions.R in the script directory if not found
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
  
  # Ensure state and year columns are standardized
  cat("Standardizing state and year columns...\n")
  # Ensure state column is uppercase character
  data$state <- toupper(as.character(data$state))
  # Ensure year column is numeric
  data$year <- as.numeric(as.character(data$year))
  
  return(data)
}

#' Generate metadata from cleaned data
#' 
#' @param data Cleaned MMWR data
#' @param source_file Original source file path
#' @param census_file_b Path to bacterial census file
#' @param census_file_p Path to parasitic census file
#' @return List with metadata information
generate_metadata <- function(data, source_file, census_file_b = NULL, census_file_p = NULL) {
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
  
  # Add preprocessed census file paths to metadata - these are state-level aggregated files
  if (!is.null(census_file_b)) {
    cat("Adding preprocessed bacterial census file path to metadata:", census_file_b, "\n")
    # Store relative path from output directory for portability
    metadata$census_file_bacterial_preprocessed <- basename(census_file_b)
    cat("Bacterial census path stored as:", metadata$census_file_bacterial_preprocessed, "\n")
  }
  
  if (!is.null(census_file_p)) {
    cat("Adding preprocessed parasitic census file path to metadata:", census_file_p, "\n")
    # Store relative path from output directory for portability  
    metadata$census_file_parasitic_preprocessed <- basename(census_file_p)
    cat("Parasitic census path stored as:", metadata$census_file_parasitic_preprocessed, "\n")
  }
  
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

  # Load and standardize census files for coverage check
  cat("Loading census files for coverage check...\n")
  
  # Helper function to standardize census files
  standardize_census <- function(file_path, pathogen_type) {
    cat(paste("Processing", pathogen_type, "census file:", file_path, "\n"))
    
    # Determine file type and read appropriately
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      census_data <- tryCatch({
        read.csv(file_path, stringsAsFactors = FALSE)
      }, error = function(e) {
        cat(paste("Error reading CSV:", e$message, "\n"))
        return(NULL)
      })
    } else if (grepl("\\.sas7bdat$", file_path, ignore.case = TRUE)) {
      census_data <- tryCatch({
        haven::read_sas(file_path)
      }, error = function(e) {
        cat(paste("Error reading SAS:", e$message, "\n"))
        return(NULL)
      })
    } else {
      cat("Unknown file format, attempting to read as SAS\n")
      census_data <- tryCatch({
        haven::read_sas(file_path)
      }, error = function(e) {
        cat(paste("Error reading file:", e$message, "\n"))
        return(NULL)
      })
    }
    
    if (is.null(census_data)) {
      stop(paste("ERROR: Failed to read", pathogen_type, "census file. Census data is required for analysis."))
    }
    
    # Print column names to help with debugging
    cat(paste(pathogen_type, "census file columns:", paste(colnames(census_data), collapse=", "), "\n"))
    
    # Find standard column names
    state_col <- grep("^state$|^st$|^STATE$|^state_name$", colnames(census_data), 
                    ignore.case = TRUE, value = TRUE)[1]
    year_col <- grep("^year$|^yr$|^YEAR$|^mmwr_year$", colnames(census_data), 
                   ignore.case = TRUE, value = TRUE)[1]
    pop_col <- grep("^population$|^pop$|^POPULATION$", colnames(census_data), 
                  ignore.case = TRUE, value = TRUE)[1]
    
    cat(paste("Found", pathogen_type, "columns - State:", state_col, "Year:", year_col, "Population:", pop_col, "\n"))
    
    # Standardize column names and types
    if (!is.na(state_col)) {
      census_data$state <- toupper(as.character(census_data[[state_col]]))
    } else {
      cat(paste("WARNING: Could not find state column in", pathogen_type, "census file\n"))
      # Error: state column is required
      stop(paste("ERROR: Could not find state column in", pathogen_type, "census file. Census data requires state column."))
    }
    
    if (!is.na(year_col)) {
      census_data$year <- as.numeric(as.character(census_data[[year_col]]))
    } else {
      cat(paste("WARNING: Could not find year column in", pathogen_type, "census file\n"))
      # Error: year column is required
      stop(paste("ERROR: Could not find year column in", pathogen_type, "census file. Census data requires year column."))
    }
    
    if (!is.na(pop_col)) {
      census_data$population <- as.numeric(as.character(census_data[[pop_col]]))
    } else {
      cat(paste("WARNING: Could not find population column in", pathogen_type, "census file\n"))
      # Error: population column is required
      stop(paste("ERROR: Could not find population column in", pathogen_type, "census file. Census data requires population column."))
    }
    
    # Add pathogentype
    census_data$pathogentype <- pathogen_type
    
    # Final check that all required columns exist and have proper types
    if (all(c("state", "year", "population", "pathogentype") %in% names(census_data))) {
      cat(paste(pathogen_type, "census file processed successfully\n"))
    } else {
      cat(paste("WARNING:", pathogen_type, "census file missing required columns after processing\n"))
      missing_cols <- setdiff(c("state", "year", "population", "pathogentype"), names(census_data))
      cat(paste("Missing columns:", paste(missing_cols, collapse=", "), "\n"))
    }
    
    return(census_data)
  }
  
  # Extract base name early for census file naming
  output_base <- tools::file_path_sans_ext(basename(args$outputFile))
  output_dir <- dirname(args$outputFile)
  
  # Process both census files
  census_b <- standardize_census(args$censusFileB, "Bacterial")
  census_p <- standardize_census(args$censusFileP, "Parasitic")
  
  # PIPELINE FIX: Aggregate county-level census to state-level to prevent join issues
  # Census files contain county-level data but analysis requires state-level totals
  # This aggregation prevents row multiplication during joins in downstream analysis
  cat("\nAggregating census data from county to state level...\n")
  
  # Aggregate bacterial census to state-year level
  census_b_state <- census_b %>%
    dplyr::group_by(state, year, pathogentype) %>%
    dplyr::summarise(
      population = sum(population, na.rm = TRUE),
      n_counties = dplyr::n(),  # Track aggregation for validation
      .groups = "drop"
    )
  cat(paste("Bacterial census: aggregated", nrow(census_b), "county records to", nrow(census_b_state), "state records\n"))
  
  # Aggregate parasitic census to state-year level  
  census_p_state <- census_p %>%
    dplyr::group_by(state, year, pathogentype) %>%
    dplyr::summarise(
      population = sum(population, na.rm = TRUE),
      n_counties = dplyr::n(),  # Track aggregation for validation
      .groups = "drop"
    )
  cat(paste("Parasitic census: aggregated", nrow(census_p), "county records to", nrow(census_p_state), "state records\n"))
  
  # Save preprocessed census files for downstream analysis
  census_b_filename <- paste0(output_base, "_census_bacterial.csv")
  census_b_path <- file.path(output_dir, census_b_filename)
  cat("Saving preprocessed bacterial census to:", census_b_path, "\n")
  write.csv(census_b_state, census_b_path, row.names = FALSE)
  
  # Save parasitic census
  census_p_filename <- paste0(output_base, "_census_parasitic.csv")
  census_p_path <- file.path(output_dir, census_p_filename)
  cat("Saving preprocessed parasitic census to:", census_p_path, "\n")
  write.csv(census_p_state, census_p_path, row.names = FALSE)
  
  # VALIDATION: Verify census aggregation worked correctly
  cat("\n=== Census Preprocessing Validation ===\n")
  
  # Check bacterial census
  if (file.exists(census_b_path)) {
    saved_census_b <- read.csv(census_b_path)
    cat("✓ Bacterial census file created:\n")
    cat("  - Records:", nrow(saved_census_b), "(from", nrow(census_b), "county records)\n")
    cat("  - States:", length(unique(saved_census_b$state)), "\n")
    cat("  - Years:", paste(range(saved_census_b$year), collapse="-"), "\n")
    cat("  - Has n_counties column:", "n_counties" %in% names(saved_census_b), "\n")
    
    # Verify no county-level columns remain
    county_cols <- grep("county|cofip|fips", tolower(names(saved_census_b)), value = TRUE)
    if (length(county_cols) > 0) {
      cat("  ⚠ WARNING: County-level columns still present:", paste(county_cols, collapse=", "), "\n")
    }
  } else {
    cat("✗ ERROR: Bacterial census file not created!\n")
  }
  
  # Check parasitic census  
  if (file.exists(census_p_path)) {
    saved_census_p <- read.csv(census_p_path)
    cat("\n✓ Parasitic census file created:\n")
    cat("  - Records:", nrow(saved_census_p), "(from", nrow(census_p), "county records)\n")
    cat("  - States:", length(unique(saved_census_p$state)), "\n")
    cat("  - Years:", paste(range(saved_census_p$year), collapse="-"), "\n")
    cat("  - Has n_counties column:", "n_counties" %in% names(saved_census_p), "\n")
    
    # Verify no county-level columns remain
    county_cols <- grep("county|cofip|fips", tolower(names(saved_census_p)), value = TRUE)
    if (length(county_cols) > 0) {
      cat("  ⚠ WARNING: County-level columns still present:", paste(county_cols, collapse=", "), "\n")
    }
  } else {
    cat("✗ ERROR: Parasitic census file not created!\n")
  }
  
  # Validate aggregation ratios
  if (nrow(census_b) > 0 && nrow(census_b_state) > 0) {
    aggregation_ratio_b <- round(nrow(census_b) / nrow(census_b_state), 1)
    cat("\n✓ Bacterial aggregation ratio:", aggregation_ratio_b, "counties per state-year\n")
    if (aggregation_ratio_b < 1.5) {
      cat("  ⚠ WARNING: Low aggregation ratio - data may already be state-level\n")
    }
  }
  
  if (nrow(census_p) > 0 && nrow(census_p_state) > 0) {
    aggregation_ratio_p <- round(nrow(census_p) / nrow(census_p_state), 1)
    cat("✓ Parasitic aggregation ratio:", aggregation_ratio_p, "counties per state-year\n")
    if (aggregation_ratio_p < 1.5) {
      cat("  ⚠ WARNING: Low aggregation ratio - data may already be state-level\n")
    }
  }
  
  cat("\n✓ Census preprocessing validation complete\n")
  cat("=====================================\n")
  
  # Combine for validation checks (using state-level data now)
  census <- dplyr::bind_rows(census_b_state, census_p_state)
  census_pairs <- unique(census[, c("state", "year")])
  cleaned_data$state <- toupper(as.character(cleaned_data$state))
  cleaned_data$year <- as.numeric(as.character(cleaned_data$year))
  mmwr_pairs <- unique(cleaned_data[, c("state", "year")])

  # Find (state, year) pairs in MMWR but not in census
  mmwr_not_in_census <- anti_join(mmwr_pairs, census_pairs, by = c("state", "year"))
  # Find (state, year) pairs in census but not in MMWR
  census_not_in_mmwr <- anti_join(census_pairs, mmwr_pairs, by = c("state", "year"))

  # Print summary to console
  cat('PREPROCESS CHECK: (state, year) pairs in MMWR but missing in census:', nrow(mmwr_not_in_census), '\n')
  if (nrow(mmwr_not_in_census) > 0) {
    cat('  These state-year combinations in your MMWR data have no matching census data:\n')
    print(mmwr_not_in_census)
  }

  cat('PREPROCESS CHECK: (state, year) pairs in census but not used in MMWR:', nrow(census_not_in_mmwr), '\n')
  if (nrow(census_not_in_mmwr) > 0) {
    cat('  These state-year combinations in your census data have no matching MMWR data:\n')
    print(census_not_in_mmwr)
  }

  # Write detailed report to file
  coverage_report_file <- file.path(dirname(args$outputFile), "state_year_coverage_report.txt")
  cat("Writing coverage report to:", coverage_report_file, "\n")
  sink(coverage_report_file)
  cat("FoodNet State-Year Coverage Report\n")
  cat("=================================\n\n")
  cat("Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  cat("MMWR Data File:", args$mmwrFile, "\n")
  cat("Census Files:", args$censusFileB, "and", args$censusFileP, "\n\n")

  cat("Summary:\n")
  cat("- Total state-year pairs in MMWR data:", nrow(mmwr_pairs), "\n")
  cat("- Total state-year pairs in census data:", nrow(census_pairs), "\n")
  cat("- State-year pairs in MMWR but missing from census:", nrow(mmwr_not_in_census), "\n")
  cat("- State-year pairs in census but not used in MMWR:", nrow(census_not_in_mmwr), "\n\n")

  if (nrow(mmwr_not_in_census) > 0) {
    cat("MMWR state-year pairs missing from census (will be DROPPED in analysis):\n")
    print(mmwr_not_in_census)
    cat("\n")
  }

  if (nrow(census_not_in_mmwr) > 0) {
    cat("Census state-year pairs not used in MMWR (informational only):\n")
    print(census_not_in_mmwr)
    cat("\n")
  }

  cat("IMPORTANT: Records with (state, year) pairs found in MMWR but not in census will be DROPPED\n")
  cat("from analysis because population data is required for rate calculations.\n\n")

  cat("End of report\n")
  sink()

  # Note: output_base and output_dir already extracted earlier for census files
  
  # Generate metadata if requested
  if(args$generate_metadata) {
    # Create standardized metadata filename
    metadata_filename <- get_output_filename(output_base, "metadata", "json")
    metadata_file <- file.path(output_dir, metadata_filename)
    
    # Generate metadata with preprocessed census file paths
    # Pass the paths to the aggregated state-level census files
    metadata <- generate_metadata(
      cleaned_data, 
      args$mmwrFile,
      census_b_path,  # Pass preprocessed bacterial census file path
      census_p_path   # Pass preprocessed parasitic census file path
    )
    
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
  
  # Highlight census file creation
  cat("Preprocessed census files created:\n")
  cat("- Bacterial census (state-level):", census_b_path, "\n")
  cat("- Parasitic census (state-level):", census_p_path, "\n")
  
  # Final preprocessing summary
  cat("\n========== Preprocessing Summary ==========\n")
  cat("MMWR data:", args$outputFile, "\n")
  cat("  - Records:", nrow(cleaned_data), "\n")
  cat("  - Pathogens:", length(unique(cleaned_data$pathogen)), "\n")
  cat("  - States:", length(unique(cleaned_data$state)), "\n")
  cat("  - Years:", paste(range(cleaned_data$year, na.rm=TRUE), collapse="-"), "\n")
  
  cat("\nCensus files (preprocessed):\n")
  cat("  - Bacterial:", census_b_path, "\n")
  cat("  - Parasitic:", census_p_path, "\n")
  
  if(args$generate_metadata && exists("metadata_file")) {
    cat("\nMetadata file:", metadata_file, "\n")
    cat("  - Contains preprocessed census paths: YES\n")
  }
  
  cat("\n✓ All files ready for analysis pipeline\n")
  cat("==========================================\n")
}

# Run the main function
main()
