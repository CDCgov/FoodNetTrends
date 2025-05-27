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
  library("haven")
  library("jsonlite")
  library("data.table")
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
          # For append mode, use fwrite with append=TRUE
          # Use quote="auto" for better compatibility with fread
          fwrite(data, file = file_path, append = TRUE, quote = "auto")
        } else {
          # For new files, use fwrite (much faster than write.table)
          # Use quote="auto" for better compatibility with fread
          fwrite(data, file = file_path, quote = "auto")
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
    # Read directly as data.table for efficiency
    data <- as.data.table(haven::read_sas(file_path))
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
  # Ensure we're working with data.table
  if (!inherits(raw_data, "data.table")) {
    setDT(raw_data)
  }
  
  # Convert all column names to lowercase for consistency
  cat("Standardizing column names to lowercase...\n")
  setnames(raw_data, tolower(names(raw_data)))
  data <- raw_data
  
  # Recode serotype values using data.table syntax
  cat("Recoding serotype values...\n")
  seroList <- c("NOT SPECIATED", "UNKNOWN", "PARTIAL SERO", "NOT SERO", "")
  data[, sero2 := ifelse(sero1 %in% seroList, "Missing", sero1)]
  data[grepl("UNDET", sero2), sero2 := "Missing"]
  data[, serotypesummary := sero2]
  
  # Standardize county names using data.table syntax
  cat("Standardizing county names...\n")
  data[county %in% c("ST. MARYS'S", "ST. MARYS"), county := "ST. MARY'S"]
  data[county == "PRINCE GEORGES", county := "PRINCE GEORGE'S"]
  data[county == "QUEEN ANNES", county := "QUEEN ANNE'S"]
  data[county == "DE BACA", county := "DEBACA"]
  
  # Ensure pathogen column is uppercase for consistency
  data[, pathogen := toupper(pathogen)]
  
  # Create pathogentype column if not present
  if(!"pathogentype" %in% names(data)) {
    cat("Creating derived pathogentype column...\n")
    data[, pathogentype := ifelse(pathogen %in% c("CRYPTOSPORIDIUM", "CYCLOSPORA"), 
                                  "Parasitic", "Bacterial")]
  }
  
  # Ensure state and year columns are standardized
  cat("Standardizing state and year columns...\n")
  # Ensure state column is uppercase character
  data[, state := toupper(as.character(state))]
  # Ensure year column is numeric
  data[, year := as.numeric(as.character(year))]
  
  return(data)
}

#' Safely infer geographic regions from state codes
#' 
#' @param state_codes Vector of state abbreviations
#' @return Vector of region names
infer_region <- function(state_codes) {
  # Define state to region mapping based on US Census regions
  region_map <- list(
    Northeast = c("CT", "ME", "MA", "NH", "RI", "VT", "NJ", "NY", "PA"),
    Midwest = c("IL", "IN", "MI", "OH", "WI", "IA", "KS", "MN", "MO", "NE", "ND", "SD"),
    South = c("DE", "FL", "GA", "MD", "NC", "SC", "VA", "DC", "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK", "TX"),
    West = c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY", "AK", "CA", "HI", "OR", "WA")
  )
  
  # Convert mapping to lookup vector
  state_to_region <- character()
  for (region in names(region_map)) {
    states <- region_map[[region]]
    state_to_region[states] <- region
  }
  
  # Apply mapping
  regions <- state_to_region[toupper(state_codes)]
  regions[is.na(regions)] <- "Unknown"
  return(regions)
}

#' Safely infer temporal groupings from dates
#' 
#' @param years Numeric vector of years
#' @param months Numeric vector of months (optional)
#' @return Data frame with temporal groupings
infer_temporal_groups <- function(years, months = NULL) {
  result <- data.frame(year = years)
  
  # Add quarter if months available
  if (!is.null(months) && length(months) == length(years)) {
    result$quarter <- ceiling(pmin(pmax(months, 1), 12) / 3)
    result$season <- ifelse(months %in% c(12, 1, 2), "Winter",
                           ifelse(months %in% c(3, 4, 5), "Spring",
                                  ifelse(months %in% c(6, 7, 8), "Summer", "Fall")))
  }
  
  # Add epidemiological groupings
  result$year_group <- cut(years, 
                          breaks = c(-Inf, 2000, 2005, 2010, 2015, 2020, Inf),
                          labels = c("Pre-2000", "2001-2005", "2006-2010", 
                                   "2011-2015", "2016-2020", "2021+"))
  
  return(result)
}

#' Safely infer age groups from numeric ages
#' 
#' @param ages Numeric vector of ages
#' @return Data frame with age groupings
infer_age_groups <- function(ages) {
  result <- data.frame(age = ages)
  
  # Standard age groups
  result$age_group_5yr <- cut(ages, 
                             breaks = c(-1, 4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79, 84, Inf),
                             labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
                                      "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", 
                                      "75-79", "80-84", "85+"))
  
  # Pediatric vs adult
  result$age_category <- ifelse(ages < 18, "Pediatric", "Adult")
  
  # Broader groups for analysis
  result$age_group_broad <- cut(ages,
                               breaks = c(-1, 4, 17, 49, 64, Inf),
                               labels = c("<5", "5-17", "18-49", "50-64", "65+"))
  
  return(result)
}

#' Apply safe inferences to cleaned data
#' 
#' @param data Cleaned MMWR data
#' @return Data with additional inferred fields
apply_safe_inferences <- function(data) {
  cat("Applying safe geographic and temporal inferences...\n")
  
  # Infer regions from states
  if ("state" %in% names(data)) {
    data$region <- infer_region(data$state)
    cat("  - Inferred regions for", sum(data$region != "Unknown"), "records\n")
  }
  
  # Infer temporal groupings
  if ("year" %in% names(data)) {
    temporal_data <- infer_temporal_groups(
      data$year, 
      if("month" %in% names(data)) data$month else NULL
    )
    data$quarter <- temporal_data$quarter
    data$season <- temporal_data$season
    data$year_group <- temporal_data$year_group
    cat("  - Added temporal groupings (quarter, season, year_group)\n")
  }
  
  # Infer age groups if age is available
  if ("age" %in% names(data)) {
    # Only apply to non-missing ages
    valid_ages <- !is.na(data$age) & data$age >= 0 & data$age <= 120
    if (sum(valid_ages) > 0) {
      age_data <- infer_age_groups(data$age[valid_ages])
      data$age_group_5yr[valid_ages] <- as.character(age_data$age_group_5yr)
      data$age_category[valid_ages] <- age_data$age_category
      data$age_group_broad[valid_ages] <- as.character(age_data$age_group_broad)
      cat("  - Inferred age groups for", sum(valid_ages), "records\n")
    }
  }
  
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
    record_count = nrow(data),
    output_file = basename(source_file)  # Will be updated with actual output filename
  )
  
  # Add inferred field information if available
  if ("region" %in% names(data)) {
    metadata$regions = sort(unique(data$region[data$region != "Unknown"]))
  }
  
  if ("age_category" %in% names(data)) {
    metadata$has_age_groups = TRUE
  }
  
  if ("quarter" %in% names(data)) {
    metadata$has_temporal_groups = TRUE
  }
  
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
    # Ensure we're working with data.table
    if (!inherits(data, "data.table")) {
      setDT(data)
    }
    
    # Use data.table to get serotype counts efficiently
    sal_data <- data[pathogen == "SALMONELLA"]
    serotype_counts <- sal_data[, .N, by = serotypesummary][order(-N)]
    
    # Store in metadata
    metadata$salmonella_serotypes <- as.list(serotype_counts$N)
    names(metadata$salmonella_serotypes) <- serotype_counts$serotypesummary
    
    # Also store a flat list of names
    metadata$salmonella_serotype_names <- serotype_counts$serotypesummary
    
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
  
  # Apply safe inferences
  cleaned_data <- apply_safe_inferences(cleaned_data)

  # Load and standardize census files for coverage check
  cat("Loading census files for coverage check...\n")
  
  # Helper function to standardize census files
  standardize_census <- function(file_path, pathogen_type) {
    cat(paste("Processing", pathogen_type, "census file:", file_path, "\n"))
    
    # Determine file type and read appropriately
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      census_data <- tryCatch({
        # Use fread for faster CSV reading
        fread(file_path, stringsAsFactors = FALSE)
      }, error = function(e) {
        cat(paste("Error reading CSV:", e$message, "\n"))
        return(NULL)
      })
    } else if (grepl("\\.sas7bdat$", file_path, ignore.case = TRUE)) {
      census_data <- tryCatch({
        # Convert to data.table after reading
        as.data.table(haven::read_sas(file_path))
      }, error = function(e) {
        cat(paste("Error reading SAS:", e$message, "\n"))
        return(NULL)
      })
    } else {
      cat("Unknown file format, attempting to read as SAS\n")
      census_data <- tryCatch({
        # Convert to data.table after reading
        as.data.table(haven::read_sas(file_path))
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
    
    # Ensure we're working with data.table
    if (!inherits(census_data, "data.table")) {
      setDT(census_data)
    }
    
    # Standardize column names and types using data.table syntax
    if (!is.na(state_col)) {
      census_data[, state := toupper(as.character(get(state_col)))]
    } else {
      cat(paste("WARNING: Could not find state column in", pathogen_type, "census file\n"))
      # Error: state column is required
      stop(paste("ERROR: Could not find state column in", pathogen_type, "census file. Census data requires state column."))
    }
    
    if (!is.na(year_col)) {
      census_data[, year := as.numeric(as.character(get(year_col)))]
    } else {
      cat(paste("WARNING: Could not find year column in", pathogen_type, "census file\n"))
      # Error: year column is required
      stop(paste("ERROR: Could not find year column in", pathogen_type, "census file. Census data requires year column."))
    }
    
    if (!is.na(pop_col)) {
      census_data[, population := as.numeric(as.character(get(pop_col)))]
    } else {
      cat(paste("WARNING: Could not find population column in", pathogen_type, "census file\n"))
      # Error: population column is required
      stop(paste("ERROR: Could not find population column in", pathogen_type, "census file. Census data requires population column."))
    }
    
    # Add pathogentype
    census_data[, pathogentype := pathogen_type]
    
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
  
  # Aggregate county-level census data to state-level totals for downstream analysis
  # Census files contain county-level data but surveillance analysis requires state-level totals
  # This aggregation ensures proper data structure for statistical modeling processes
  cat("\nAggregating census data from county to state level...\n")
  
  # Convert to data.table for efficient aggregation
  if (!inherits(census_b, "data.table")) {
    setDT(census_b)
  }
  if (!inherits(census_p, "data.table")) {
    setDT(census_p)
  }
  
  # Save row counts for validation before aggregation
  census_b_rows <- nrow(census_b)
  census_p_rows <- nrow(census_p)
  
  # Aggregate bacterial census to state-year level using data.table
  census_b_state <- census_b[, .(
    population = sum(population, na.rm = TRUE),
    n_counties = .N  # Track aggregation for validation
  ), by = .(state, year, pathogentype)]
  
  cat(paste("Bacterial census: aggregated", census_b_rows, "county records to", nrow(census_b_state), "state records\n"))
  
  # Clean up original census_b to free memory
  rm(census_b)
  gc()
  
  # Aggregate parasitic census to state-year level using data.table
  census_p_state <- census_p[, .(
    population = sum(population, na.rm = TRUE),
    n_counties = .N  # Track aggregation for validation
  ), by = .(state, year, pathogentype)]
  
  cat(paste("Parasitic census: aggregated", census_p_rows, "county records to", nrow(census_p_state), "state records\n"))
  
  # Clean up original census_p to free memory
  rm(census_p)
  gc()
  
  # Save preprocessed census files for downstream analysis
  census_b_filename <- paste0(output_base, "_census_bacterial.csv")
  census_b_path <- file.path(output_dir, census_b_filename)
  cat("Saving preprocessed bacterial census to:", census_b_path, "\n")
  fwrite(census_b_state, census_b_path, quote = "auto")
  
  # Save parasitic census
  census_p_filename <- paste0(output_base, "_census_parasitic.csv")
  census_p_path <- file.path(output_dir, census_p_filename)
  cat("Saving preprocessed parasitic census to:", census_p_path, "\n")
  fwrite(census_p_state, census_p_path, quote = "auto")
  
  # VALIDATION: Verify census aggregation worked correctly
  cat("\n=== Census Preprocessing Validation ===\n")
  
  # Check bacterial census
  if (file.exists(census_b_path)) {
    saved_census_b <- fread(census_b_path)
    cat("✓ Bacterial census file created:\n")
    cat("  - Records:", nrow(saved_census_b), "(from", census_b_rows, "county records)\n")
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
    saved_census_p <- fread(census_p_path)
    cat("\n✓ Parasitic census file created:\n")
    cat("  - Records:", nrow(saved_census_p), "(from", census_p_rows, "county records)\n")
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
  
  # Validate aggregation ratios using saved row counts
  if (census_b_rows > 0 && nrow(census_b_state) > 0) {
    aggregation_ratio_b <- round(census_b_rows / nrow(census_b_state), 1)
    cat("\n✓ Bacterial aggregation ratio:", aggregation_ratio_b, "counties per state-year\n")
    if (aggregation_ratio_b < 1.5) {
      cat("  ⚠ WARNING: Low aggregation ratio - data may already be state-level\n")
    }
  }
  
  if (census_p_rows > 0 && nrow(census_p_state) > 0) {
    aggregation_ratio_p <- round(census_p_rows / nrow(census_p_state), 1)
    cat("✓ Parasitic aggregation ratio:", aggregation_ratio_p, "counties per state-year\n")
    if (aggregation_ratio_p < 1.5) {
      cat("  ⚠ WARNING: Low aggregation ratio - data may already be state-level\n")
    }
  }
  
  cat("\n✓ Census preprocessing validation complete\n")
  cat("=====================================\n")
  
  # Combine for validation checks (using state-level data now)
  # Use rbindlist for efficient combination
  census <- rbindlist(list(census_b_state, census_p_state), use.names = TRUE)
  
  # Ensure cleaned_data is data.table
  if (!inherits(cleaned_data, "data.table")) {
    setDT(cleaned_data)
  }
  
  # Get unique pairs efficiently
  census_pairs <- unique(census[, .(state, year)])
  cleaned_data[, `:=`(state = toupper(as.character(state)), 
                      year = as.numeric(as.character(year)))]
  mmwr_pairs <- unique(cleaned_data[, .(state, year)])

  # Find (state, year) pairs in MMWR but not in census using data.table anti-join
  setkey(mmwr_pairs, state, year)
  setkey(census_pairs, state, year)
  mmwr_not_in_census <- mmwr_pairs[!census_pairs]
  
  # Find (state, year) pairs in census but not in MMWR
  census_not_in_mmwr <- census_pairs[!mmwr_pairs]

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

  # Check for population consistency across pathogen types
  cat('\nPREPROCESS CHECK: Verifying population consistency across pathogen types...\n')
  
  # Compare bacterial and parasitic census for same state-year combinations using data.table
  # Rename columns for clarity
  census_b_compare <- census_b_state[, .(state, year, pop_bacterial = population)]
  census_p_compare <- census_p_state[, .(state, year, pop_parasitic = population)]
  
  # Inner join using data.table
  setkey(census_b_compare, state, year)
  setkey(census_p_compare, state, year)
  pop_comparison <- census_b_compare[census_p_compare, nomatch = 0]
  
  # Find discrepancies using data.table syntax
  pop_comparison[, `:=`(
    pop_diff = abs(pop_bacterial - pop_parasitic),
    pct_diff = round(100 * abs(pop_bacterial - pop_parasitic) / pmax(pop_bacterial, pop_parasitic), 2)
  )]
  
  # Flag significant discrepancies (>1% difference)
  discrepancies <- pop_comparison[pct_diff > 1]
  
  if (nrow(discrepancies) > 0) {
    cat('  WARNING: Population inconsistencies found between bacterial and parasitic census data!\n')
    cat('  State-year combinations with >1% population difference:\n')
    print(discrepancies[order(-pct_diff), .(state, year, pop_bacterial, pop_parasitic, pct_diff)])
    cat('\n  This indicates potential data quality issues. Population should be consistent for the same state-year.\n')
    cat('  Consider reviewing the source census files for accuracy.\n')
  } else {
    cat('  ✓ Population data is consistent across pathogen types (all differences <1%)\n')
    cat('  This confirms population inference across pathogen types is valid.\n')
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
  
  cat("Population Consistency Check:\n")
  cat("==============================\n")
  if (exists("discrepancies") && nrow(discrepancies) > 0) {
    cat("WARNING: Inconsistent population values found between bacterial and parasitic census data\n")
    cat("The following state-year combinations show >1% difference:\n\n")
    print(discrepancies[order(-discrepancies$pct_diff), ])
    cat("\nRecommendation: Review source census files for data quality issues\n")
  } else {
    cat("✓ All population values are consistent between bacterial and parasitic census data\n")
    cat("  Maximum difference: <1%\n") 
    cat("  This validates that population inference across pathogen types is appropriate\n")
  }
  cat("\n")

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
    
    # Add original census file paths to metadata for full traceability
    metadata$census_file_bacterial_original <- args$censusFileB
    metadata$census_file_parasitic_original <- args$censusFileP
    
    # Update output file name to actual output file
    metadata$output_file <- basename(args$outputFile)
    
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
  
  # Clean up memory before exiting
  rm(cleaned_data, census_b_state, census_p_state, census, mmwrdata)
  gc()
}

# Run the main function
main()
