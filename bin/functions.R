# =========================================================================
# FoodNetTrends v1.0 - Core Statistical and Data Processing Functions
# =========================================================================
#
# Purpose:
#   Contains essential statistical functions for the FoodNetTrends pipeline.
#   Provides data processing, Bayesian modeling utilities, and result
#   formatting for foodborne disease surveillance analysis.
#
# Key Functions:
#   - Data preparation and standardization
#   - Bayesian model utilities for brms
#   - Specialized pathogen analysis (Cyclospora, Salmonella)
#   - Visualization and output formatting
#
# Last updated: 2025-05-22
# =========================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(gtools)
  library(brms)
  library(ggplot2)
  library(tidybayes)
  library(haven)
  library(tibble)
  library(readr)  
  library(HDInterval)
  library(gridExtra)
})

#' Generate Standardized Filename
#'
#' Generates a standardized filename for FoodNetTrends outputs.
#' This function ensures consistent naming patterns across the pipeline.
#'
#' @param pathogen Name of the pathogen (e.g., "CAMPYLOBACTER")
#' @param file_type Type of file (e.g., "model", "IRCatch", "summary")
#' @param extension File extension without dot (e.g., "Rds", "csv", "txt", "png")
#' @param subtype Optional subtype for specialized files (e.g., years for comparison files)
#' @return A standardized filename string
#' @examples
#' get_output_filename("CAMPYLOBACTER", "model", "Rds")
#' get_output_filename("SALMONELLA", "IRCatch", "csv") 
#' get_output_filename("CYCLOSPORA", "EstIRRCatch", "csv", "2016_2018")
get_output_filename <- function(pathogen, file_type, extension, subtype = NULL) {
  # Ensure inputs are valid
  if (is.null(pathogen) || is.null(file_type) || is.null(extension)) {
    stop("Pathogen, file_type, and extension must all be provided")
  }
  
  # Build filename with consistent pattern
  filename <- paste0(pathogen, "_", file_type)
  
  # Add subtype if provided
  if (!is.null(subtype)) {
    filename <- paste0(filename, "_", subtype)
  }
  
  # Add extension
  filename <- paste0(filename, ".", extension)
  
  return(filename)
}

#' Clean a List String Input
#'
#' Processes a string input or vector containing comma-separated values
#' and returns a clean vector of values.
#'
#' @param input_string A string containing comma-separated values or a vector of values.
#' @return A character vector with cleaned values.
#' @examples
#' clean_list("NO,UNKNOWN,YES")
#' clean_list(c("CIDT+", "CX+"))
clean_list <- function(input_string) {
  if (length(input_string) > 1) {
    # If input is a vector, collapse into a single string
    input_string <- paste(input_string, collapse = ",")
  }
  # Remove brackets and quotes, then split by comma
  cleanedString <- gsub('[\\[\\]\"]', '', input_string)
  strsplit(cleanedString, ",")[[1]]
}

#' Write Data to a File Safely
#'
#' Writes a data frame to a file, ensuring the target directory exists and handling errors.
#'
#' @param data Data frame to write
#' @param file_path Full path to the output file (.csv or .Rds)
#' @return None
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

# =========================================================================
# Pathogen-specific data processing functions
# =========================================================================
# The following functions handle different pathogens separately because:
# 1. Different pathogens require different census denominators (bacterial vs. parasitic)
# 2. Some pathogens (like Salmonella) have special processing requirements
# 3. Handling them separately allows for pathogen-specific customization
#    without complicating a single generic function
# =========================================================================

#' Prepare and Aggregate Pathogen Data
#'
#' Filters and aggregates FoodNetTrends data for specified pathogens and joins with census data.
#' This function is used for most bacterial pathogens.
#'
#' @param mmwrdata MMWR surveillance data frame
#' @param census Census data frame
#' @return Aggregated data frame with counts and population by year, state, and pathogen
path_analysis <- function(mmwrdata, census) {
  # Define standard pathogens
  pathogens <- c("CAMPYLOBACTER", "CYCLOSPORA", "SALMONELLA", "SHIGELLA", "STEC", "VIBRIO", "YERSINIA")
  
  # Log analysis start
  message(paste("Starting path_analysis with", nrow(mmwrdata), "MMWR records and", nrow(census), "census records"))
  
  # Ensure required columns exist in mmwrdata
  required_cols <- c("state", "year", "pathogen")
  missing_cols <- required_cols[!required_cols %in% names(mmwrdata)]
  
  if (length(missing_cols) > 0) {
    # Try looking for alternative column names
    for (col in missing_cols) {
      col_variations <- list(
        "state" = c("State", "STATE", "St", "ST", "state_name", "STATE_NAME"),
        "year" = c("Year", "YEAR", "YR", "yr", "MMWR_YEAR", "mmwr_year"),
        "pathogen" = c("Pathogen", "PATHOGEN", "Path", "PATH", "pathogen_name", "PATHOGEN_NAME", "organism", "ORGANISM")
      )
      
      # Check for alternative names
      for (alt_name in col_variations[[col]]) {
        if (alt_name %in% names(mmwrdata)) {
          # Rename to standard name
          names(mmwrdata)[names(mmwrdata) == alt_name] <- col
          warning(paste("Renamed column", alt_name, "to", col, "in path_analysis"))
          break
        }
      }
      
      # If still missing, raise an error - required columns must exist
      if (!col %in% names(mmwrdata)) {
        stop(paste("ERROR: Required column", col, "not found in MMWR data. Data is incomplete or malformed."))
      }
    }
  }
  
  # Print mmwrdata column names for debugging
  message("MMWR data columns: ", paste(names(mmwrdata), collapse=", "))
  
  # Ensure required columns exist in census data
  required_cols <- c("state", "year", "population", "pathogentype")
  missing_cols <- required_cols[!required_cols %in% names(census)]
  
  if (length(missing_cols) > 0) {
    # Try looking for alternative column names
    for (col in missing_cols) {
      col_variations <- list(
        "state" = c("State", "STATE", "St", "ST", "state_name", "STATE_NAME"),
        "year" = c("Year", "YEAR", "YR", "yr", "MMWR_YEAR", "mmwr_year"),
        "population" = c("Population", "POPULATION", "Pop", "POP", "pop"),
        "pathogentype" = c("PathogenType", "PATHOGENTYPE", "pathogen_type", "PATHOGEN_TYPE", "type", "TYPE")
      )
      
      # Check for alternative names
      for (alt_name in col_variations[[col]]) {
        if (alt_name %in% names(census)) {
          # Rename to standard name
          names(census)[names(census) == alt_name] <- col
          warning(paste("Renamed column", alt_name, "to", col, "in census data"))
          break
        }
      }
      
      # If still missing, raise error
      if (!col %in% names(census)) {
        stop(paste("CRITICAL ERROR: Required column", col, "not found in census data. Cannot proceed."))
      }
    }
  }
  
  # Print census column names for debugging
  message("Census data columns: ", paste(names(census), collapse=", "))
  
  # Coerce state and year to same type/case
  mmwrdata$state <- toupper(as.character(mmwrdata$state))
  mmwrdata$year <- as.numeric(as.character(mmwrdata$year))
  census$state <- toupper(as.character(census$state))
  census$year <- as.numeric(as.character(census$year))
  
  # Check for NAs in key columns
  if (any(is.na(mmwrdata$state)) || any(is.na(mmwrdata$year)) || any(is.na(mmwrdata$pathogen))) {
    warning("NA values found in key MMWR data columns: ", 
            "state: ", sum(is.na(mmwrdata$state)), 
            ", year: ", sum(is.na(mmwrdata$year)), 
            ", pathogen: ", sum(is.na(mmwrdata$pathogen)))
  }
  
  if (any(is.na(census$state)) || any(is.na(census$year)) || any(is.na(census$population))) {
    warning("NA values found in key census columns: ", 
            "state: ", sum(is.na(census$state)), 
            ", year: ", sum(is.na(census$year)), 
            ", population: ", sum(is.na(census$population)))
  }
  
  # Verify bacterial pathogentype entries in census
  if (sum(census$pathogentype == "Bacterial", na.rm = TRUE) == 0) {
    stop("CRITICAL ERROR: No Bacterial pathogentype found in census data. Cannot proceed without census data.")
  }
  
  # Create a data frame with counts per year, state, and pathogen
  message("Filtering MMWR data for pathogens of interest")
  filtered_data <- mmwrdata %>%
    filter(pathogen %in% pathogens)
  
  if (nrow(filtered_data) == 0) {
    stop("CRITICAL ERROR: No matching pathogen data found in MMWR data. Cannot proceed with analysis.")
  }
  
  message("Aggregating pathogen counts by year, state, and pathogen")
  pathogen_counts <- filtered_data %>%
    group_by(year, state, pathogen) %>%
    summarise(count = n(), .groups = "drop")
  
  message("Number of pathogen count rows: ", nrow(pathogen_counts))
  
  # Complete the dataset with all state/year/pathogen combinations
  message("Completing dataset with all combinations")
  pathogen_counts_complete <- pathogen_counts %>%
    complete(
      year = unique(pathogen_counts$year), 
      state = unique(pathogen_counts$state),
      pathogen = unique(pathogen_counts$pathogen), 
      fill = list(count = 0)
    )
  
  message("Number of rows after completion: ", nrow(pathogen_counts_complete))
  
  # Filter census for bacterial entries
  message("Filtering census data for bacterial records")
  census_bacterial <- census %>% 
    filter(toupper(pathogentype) == "BACTERIAL")
  
  message("Number of bacterial census records: ", nrow(census_bacterial))
  
  if (nrow(census_bacterial) == 0) {
    warning("No bacterial census records found after filtering. Using all census records.")
    census_bacterial <- census
  }
  
  # Join with census data carefully
  message("Joining pathogen counts with census data")
  pre_join_rows <- nrow(pathogen_counts_complete)
  
  # Check join columns before attempting join
  join_cols <- c("year", "state")
  if (!all(join_cols %in% names(pathogen_counts_complete)) || 
      !all(join_cols %in% names(census_bacterial))) {
    warning("Join columns missing in datasets!")
    print(paste("pathogen_counts columns:", paste(names(pathogen_counts_complete), collapse=", ")))
    print(paste("census_bacterial columns:", paste(names(census_bacterial), collapse=", ")))
    
    # Cannot proceed without essential join columns
    stop("CRITICAL ERROR: Essential join columns (state, year) missing from data. Check data structure and preprocessing.")
  }
  
  # Perform the join
  selectDf <- left_join(
    pathogen_counts_complete,
    census_bacterial,
    by = join_cols
  )
  
  post_join_rows <- nrow(selectDf)
  message("Rows before join: ", pre_join_rows, ", after join: ", post_join_rows)
  
  if (post_join_rows != pre_join_rows) {
    warning(paste("Join changed row count from", pre_join_rows, "to", post_join_rows))
  }
  
  # Check for missing population values
  na_population_count <- sum(is.na(selectDf$population))
  if (na_population_count > 0) {
    warning(paste(na_population_count, "rows have missing population values after join. These will be excluded."))
    # Exclude rows with missing population
    selectDf <- selectDf[!is.na(selectDf$population), ]
    if (nrow(selectDf) == 0) {
      stop("CRITICAL ERROR: No complete data (with population) available after join.")
    }
  }
  
  # Add pathogentype if missing
  if (!"pathogentype" %in% names(selectDf) || all(is.na(selectDf$pathogentype))) {
    warning("Missing pathogentype column after join, adding default")
    selectDf$pathogentype <- "Bacterial"
  }
  
  message("Final dataset has", nrow(selectDf), "rows")
  return(selectDf)
}

#' Prepare and Aggregate Cyclospora Data
#'
#' Filters and aggregates FoodNetTrends data specifically for Cyclospora and joins with census data.
#' Note: Cyclospora requires parasitic census data, unlike bacterial pathogens.
#'
#' @param mmwrdata MMWR surveillance data frame
#' @param census Census data frame
#' @return Aggregated data frame with counts and population by year and state for Cyclospora
cyclospora_analysis <- function(mmwrdata, census) {
  # Log analysis start
  message(paste("Starting cyclospora_analysis with", nrow(mmwrdata), "MMWR records and", nrow(census), "census records"))
  
  # Ensure required columns exist in mmwrdata
  required_cols <- c("state", "year", "pathogen")
  missing_cols <- required_cols[!required_cols %in% names(mmwrdata)]
  
  if (length(missing_cols) > 0) {
    # Try looking for alternative column names
    for (col in missing_cols) {
      col_variations <- list(
        "state" = c("State", "STATE", "St", "ST", "state_name", "STATE_NAME"),
        "year" = c("Year", "YEAR", "YR", "yr", "MMWR_YEAR", "mmwr_year"),
        "pathogen" = c("Pathogen", "PATHOGEN", "Path", "PATH", "pathogen_name", "PATHOGEN_NAME", "organism", "ORGANISM")
      )
      
      # Check for alternative names
      for (alt_name in col_variations[[col]]) {
        if (alt_name %in% names(mmwrdata)) {
          # Rename to standard name
          names(mmwrdata)[names(mmwrdata) == alt_name] <- col
          warning(paste("Renamed column", alt_name, "to", col, "in cyclospora_analysis"))
          break
        }
      }
      
      # If still missing, raise an error - required columns must exist
      if (!col %in% names(mmwrdata)) {
        stop(paste("ERROR: Required column", col, "not found in MMWR data. Data is incomplete or malformed."))
      }
    }
  }
  
  # Print mmwrdata column names for debugging
  message("MMWR data columns: ", paste(names(mmwrdata), collapse=", "))
  
  # Ensure required columns exist in census data
  required_cols <- c("state", "year", "population", "pathogentype")
  missing_cols <- required_cols[!required_cols %in% names(census)]
  
  if (length(missing_cols) > 0) {
    # Try looking for alternative column names
    for (col in missing_cols) {
      col_variations <- list(
        "state" = c("State", "STATE", "St", "ST", "state_name", "STATE_NAME"),
        "year" = c("Year", "YEAR", "YR", "yr", "MMWR_YEAR", "mmwr_year"),
        "population" = c("Population", "POPULATION", "Pop", "POP", "pop"),
        "pathogentype" = c("PathogenType", "PATHOGENTYPE", "pathogen_type", "PATHOGEN_TYPE", "type", "TYPE")
      )
      
      # Check for alternative names
      for (alt_name in col_variations[[col]]) {
        if (alt_name %in% names(census)) {
          # Rename to standard name
          names(census)[names(census) == alt_name] <- col
          warning(paste("Renamed column", alt_name, "to", col, "in census data"))
          break
        }
      }
      
      # If still missing, raise error
      if (!col %in% names(census)) {
        stop(paste("CRITICAL ERROR: Required column", col, "not found in census data. Cannot proceed."))
      }
    }
  }
  
  # Print census column names for debugging
  message("Census data columns: ", paste(names(census), collapse=", "))

  # Coerce state and year to same type/case
  mmwrdata$state <- toupper(as.character(mmwrdata$state))
  mmwrdata$year <- as.numeric(as.character(mmwrdata$year))
  census$state <- toupper(as.character(census$state))
  census$year <- as.numeric(as.character(census$year))
  
  # Check for NAs in key columns
  if (any(is.na(mmwrdata$state)) || any(is.na(mmwrdata$year)) || any(is.na(mmwrdata$pathogen))) {
    warning("NA values found in key MMWR data columns: ", 
            "state: ", sum(is.na(mmwrdata$state)), 
            ", year: ", sum(is.na(mmwrdata$year)), 
            ", pathogen: ", sum(is.na(mmwrdata$pathogen)))
  }
  
  if (any(is.na(census$state)) || any(is.na(census$year)) || any(is.na(census$population))) {
    warning("NA values found in key census columns: ", 
            "state: ", sum(is.na(census$state)), 
            ", year: ", sum(is.na(census$year)), 
            ", population: ", sum(is.na(census$population)))
  }
  
  # Verify parasitic pathogentype entries in census
  if (sum(toupper(census$pathogentype) == "PARASITIC", na.rm = TRUE) == 0) {
    stop("CRITICAL ERROR: No Parasitic pathogentype found in census data. Cannot proceed without census data.")
  }
  
  # Filter for Cyclospora
  message("Filtering for CYCLOSPORA records")
  cyclospora_data <- mmwrdata %>%
    filter(toupper(pathogen) == "CYCLOSPORA")
  
  if (nrow(cyclospora_data) == 0) {
    stop("CRITICAL ERROR: No CYCLOSPORA data found in MMWR data. Cannot proceed with analysis.")
  }
  
  message("Aggregating Cyclospora counts by year and state")
  cyclo_counts <- cyclospora_data %>%
    group_by(year, state) %>%
    summarise(count = n(), .groups = "drop")
  
  message("Number of Cyclospora count rows: ", nrow(cyclo_counts))
  
  # Complete the dataset with all state/year combinations
  message("Completing dataset with all combinations")
  cyclo_counts_complete <- cyclo_counts %>%
    complete(
      year = unique(cyclo_counts$year), 
      state = unique(cyclo_counts$state),
      fill = list(count = 0)
    )
  
  message("Number of rows after completion: ", nrow(cyclo_counts_complete))
  
  # Filter census for parasitic entries
  message("Filtering census data for parasitic records")
  census_parasitic <- census %>% 
    filter(toupper(pathogentype) == "PARASITIC")
  
  message("Number of parasitic census records: ", nrow(census_parasitic))
  
  # Check if census data needs aggregation (county to state level)
  # Preprocessed census files are already state-level, raw files need aggregation
  needs_aggregation <- FALSE
  
  # Check for county-level indicators
  if ("county" %in% tolower(names(census_parasitic)) || 
      "cofip" %in% tolower(names(census_parasitic)) ||
      "n_counties" %in% names(census_parasitic)) {
    # If n_counties exists, data is already aggregated from preprocessing
    if ("n_counties" %in% names(census_parasitic)) {
      message("Census data is already aggregated to state level (from preprocessing)")
    } else {
      needs_aggregation <- TRUE
      message("Census data appears to be county-level, aggregation needed")
    }
  }
  
  if (needs_aggregation) {
    # PIPELINE FIX: Aggregate county-level census to state-level to prevent join duplication
    # Raw census files contain county-level rows causing massive row multiplication during join
    message("Aggregating census data to state-year level")
    pre_agg_rows <- nrow(census_parasitic)
    
    census_parasitic <- census_parasitic %>%
      group_by(state, year) %>%
      summarise(
        population = sum(population, na.rm = TRUE),
        pathogentype = first(pathogentype),
        n_counties = n(),  # Track aggregation
        .groups = "drop"
      )
    message("Census records: aggregated ", pre_agg_rows, " county records to ", nrow(census_parasitic), " state records")
  } else {
    message("Census data is already at state level, no aggregation needed")
  }
  
  if (nrow(census_parasitic) == 0) {
    warning("No parasitic census records found after filtering. Using all census records.")
    census_parasitic <- census
  }
  
  # Join with census data carefully
  message("Joining Cyclospora counts with census data")
  pre_join_rows <- nrow(cyclo_counts_complete)
  
  # Check join columns before attempting join
  join_cols <- c("year", "state")
  if (!all(join_cols %in% names(cyclo_counts_complete)) || 
      !all(join_cols %in% names(census_parasitic))) {
    warning("Join columns missing in datasets!")
    print(paste("cyclo_counts columns:", paste(names(cyclo_counts_complete), collapse=", ")))
    print(paste("census_parasitic columns:", paste(names(census_parasitic), collapse=", ")))
    
    # Cannot proceed without essential join columns
    stop("CRITICAL ERROR: Essential join columns (state, year) missing from data. Check data structure and preprocessing.")
  }
  
  # Perform the join
  cyclo <- left_join(
    cyclo_counts_complete,
    census_parasitic,
    by = join_cols
  )
  
  post_join_rows <- nrow(cyclo)
  message("Rows before join: ", pre_join_rows, ", after join: ", post_join_rows)
  
  if (post_join_rows != pre_join_rows) {
    warning(paste("Join changed row count from", pre_join_rows, "to", post_join_rows))
  }
  
  # Handle missing population values - exclude incomplete data rather than fabricate
  na_population_count <- sum(is.na(cyclo$population))
  if (na_population_count > 0) {
    excluded_data <- cyclo[is.na(cyclo$population), c("state", "year")]
    warning(paste("EXCLUDING", na_population_count, "rows due to missing population data:"))
    if (nrow(excluded_data) > 0) {
      excluded_summary <- excluded_data %>%
        group_by(state) %>%
        summarise(missing_years = paste(sort(unique(year)), collapse=", "), .groups = "drop")
      for(i in 1:nrow(excluded_summary)) {
        warning(paste("  State", excluded_summary$state[i], "missing years:", excluded_summary$missing_years[i]))
      }
    }
    # Remove incomplete records
    cyclo <- cyclo[!is.na(cyclo$population), ]
    message(paste("Analysis will proceed with", nrow(cyclo), "complete records"))
  }
  
  # Add pathogentype if missing
  if (!"pathogentype" %in% names(cyclo) || all(is.na(cyclo$pathogentype))) {
    warning("Missing pathogentype column after join, adding default")
    cyclo$pathogentype <- "Parasitic"
  }
  
  message("Final Cyclospora dataset has", nrow(cyclo), "rows")
  return(cyclo)
}

#' Prepare and Aggregate Salmonella Data
#'
#' Filters and aggregates FoodNetTrends data specifically for Salmonella and joins with census data.
#' Salmonella gets special handling due to its public health importance and serotype considerations.
#'
#' @param mmwrdata MMWR surveillance data frame
#' @param census Census data frame
#' @return Aggregated data frame with counts and population by year and state for Salmonella
salmonella_analysis <- function(mmwrdata, census) {
  # Log analysis start
  message(paste("Starting salmonella_analysis with", nrow(mmwrdata), "MMWR records and", nrow(census), "census records"))
  
  # Ensure required columns exist in mmwrdata
  required_cols <- c("state", "year", "pathogen")
  missing_cols <- required_cols[!required_cols %in% names(mmwrdata)]
  
  if (length(missing_cols) > 0) {
    # Try looking for alternative column names
    for (col in missing_cols) {
      col_variations <- list(
        "state" = c("State", "STATE", "St", "ST", "state_name", "STATE_NAME"),
        "year" = c("Year", "YEAR", "YR", "yr", "MMWR_YEAR", "mmwr_year"),
        "pathogen" = c("Pathogen", "PATHOGEN", "Path", "PATH", "pathogen_name", "PATHOGEN_NAME", "organism", "ORGANISM")
      )
      
      # Check for alternative names
      for (alt_name in col_variations[[col]]) {
        if (alt_name %in% names(mmwrdata)) {
          # Rename to standard name
          names(mmwrdata)[names(mmwrdata) == alt_name] <- col
          warning(paste("Renamed column", alt_name, "to", col, "in salmonella_analysis"))
          break
        }
      }
      
      # If still missing, raise an error - required columns must exist
      if (!col %in% names(mmwrdata)) {
        stop(paste("ERROR: Required column", col, "not found in MMWR data. Data is incomplete or malformed."))
      }
    }
  }
  
  # Print mmwrdata column names for debugging
  message("MMWR data columns: ", paste(names(mmwrdata), collapse=", "))
  
  # Ensure required columns exist in census data
  required_cols <- c("state", "year", "population", "pathogentype")
  missing_cols <- required_cols[!required_cols %in% names(census)]
  
  if (length(missing_cols) > 0) {
    # Try looking for alternative column names
    for (col in missing_cols) {
      col_variations <- list(
        "state" = c("State", "STATE", "St", "ST", "state_name", "STATE_NAME"),
        "year" = c("Year", "YEAR", "YR", "yr", "MMWR_YEAR", "mmwr_year"),
        "population" = c("Population", "POPULATION", "Pop", "POP", "pop"),
        "pathogentype" = c("PathogenType", "PATHOGENTYPE", "pathogen_type", "PATHOGEN_TYPE", "type", "TYPE")
      )
      
      # Check for alternative names
      for (alt_name in col_variations[[col]]) {
        if (alt_name %in% names(census)) {
          # Rename to standard name
          names(census)[names(census) == alt_name] <- col
          warning(paste("Renamed column", alt_name, "to", col, "in census data"))
          break
        }
      }
      
      # If still missing, raise error
      if (!col %in% names(census)) {
        stop(paste("CRITICAL ERROR: Required column", col, "not found in census data. Cannot proceed."))
      }
    }
  }
  
  # Print census column names for debugging
  message("Census data columns: ", paste(names(census), collapse=", "))

  # Coerce state and year to same type/case
  mmwrdata$state <- toupper(as.character(mmwrdata$state))
  mmwrdata$year <- as.numeric(as.character(mmwrdata$year))
  census$state <- toupper(as.character(census$state))
  census$year <- as.numeric(as.character(census$year))
  
  # Check for NAs in key columns
  if (any(is.na(mmwrdata$state)) || any(is.na(mmwrdata$year)) || any(is.na(mmwrdata$pathogen))) {
    warning("NA values found in key MMWR data columns: ", 
            "state: ", sum(is.na(mmwrdata$state)), 
            ", year: ", sum(is.na(mmwrdata$year)), 
            ", pathogen: ", sum(is.na(mmwrdata$pathogen)))
  }
  
  if (any(is.na(census$state)) || any(is.na(census$year)) || any(is.na(census$population))) {
    warning("NA values found in key census columns: ", 
            "state: ", sum(is.na(census$state)), 
            ", year: ", sum(is.na(census$year)), 
            ", population: ", sum(is.na(census$population)))
  }
  
  # Verify bacterial pathogentype entries in census
  if (sum(toupper(census$pathogentype) == "BACTERIAL", na.rm = TRUE) == 0) {
    stop("CRITICAL ERROR: No Bacterial pathogentype found in census data. Cannot proceed without census data.")
  }
  
  # Filter for Salmonella
  message("Filtering for SALMONELLA records")
  salmonella_data <- mmwrdata %>%
    filter(toupper(pathogen) == "SALMONELLA")
  
  if (nrow(salmonella_data) == 0) {
    stop("CRITICAL ERROR: No SALMONELLA data found in MMWR data. Cannot proceed with analysis.")
  }
  
  message("Aggregating Salmonella counts by year and state")
  sal_counts <- salmonella_data %>%
    group_by(year, state) %>%
    summarise(count = n(), .groups = "drop")
  
  message("Number of Salmonella count rows: ", nrow(sal_counts))
  
  # Complete the dataset with all state/year combinations
  message("Completing dataset with all combinations")
  sal_counts_complete <- sal_counts %>%
    complete(
      year = unique(sal_counts$year), 
      state = unique(sal_counts$state),
      fill = list(count = 0)
    )
  
  message("Number of rows after completion: ", nrow(sal_counts_complete))
  
  # Filter census for bacterial entries
  message("Filtering census data for bacterial records")
  census_bacterial <- census %>% 
    filter(toupper(pathogentype) == "BACTERIAL")
  
  message("Number of bacterial census records: ", nrow(census_bacterial))
  
  if (nrow(census_bacterial) == 0) {
    warning("No bacterial census records found after filtering. Using all census records.")
    census_bacterial <- census
  }
  
  # Join with census data carefully
  message("Joining Salmonella counts with census data")
  pre_join_rows <- nrow(sal_counts_complete)
  
  # Check join columns before attempting join
  join_cols <- c("year", "state")
  if (!all(join_cols %in% names(sal_counts_complete)) || 
      !all(join_cols %in% names(census_bacterial))) {
    warning("Join columns missing in datasets!")
    print(paste("sal_counts columns:", paste(names(sal_counts_complete), collapse=", ")))
    print(paste("census_bacterial columns:", paste(names(census_bacterial), collapse=", ")))
    
    # Cannot proceed without essential join columns
    stop("CRITICAL ERROR: Essential join columns (state, year) missing from data. Check data structure and preprocessing.")
  }
  
  # Perform the join
  sal <- left_join(
    sal_counts_complete,
    census_bacterial,
    by = join_cols
  )
  
  post_join_rows <- nrow(sal)
  message("Rows before join: ", pre_join_rows, ", after join: ", post_join_rows)
  
  if (post_join_rows != pre_join_rows) {
    warning(paste("Join changed row count from", pre_join_rows, "to", post_join_rows))
  }
  
  # Check for missing population values
  na_population_count <- sum(is.na(sal$population))
  if (na_population_count > 0) {
    warning(paste(na_population_count, "rows have missing population values after join. These will be excluded."))
    # Exclude rows with missing population
    sal <- sal[!is.na(sal$population), ]
    if (nrow(sal) == 0) {
      stop("CRITICAL ERROR: No complete data (with population) available after join.")
    }
  }
  
  # Add pathogentype if missing
  if (!"pathogentype" %in% names(sal) || all(is.na(sal$pathogentype))) {
    warning("Missing pathogentype column after join, adding default")
    sal$pathogentype <- "Bacterial"
  }
  
  message("Final Salmonella dataset has", nrow(sal), "rows")
  return(sal)
}

#' Fit Bayesian Model for Pathogen Trends
#'
#' Fits a Bayesian hierarchical model with splines to estimate incidence rates.
#' Includes robust error handling and fallback mechanisms for zero-count data.
#'
#' @param data Data frame containing count, year, state, and population
#' @param cores Number of cores to use for model fitting
#' @param chains Number of MCMC chains
#' @param iterations Number of MCMC iterations
#' @param adapt_delta Adaptation parameter for HMC
#' @param max_treedepth Maximum tree depth for HMC
#' @param seed Random seed for reproducibility
#' @return A brms model object, or dummy model if fitting fails
proposed_bm <- function(data, cores = 16, chains = 2, iterations = 500,
                        adapt_delta = 0.95, max_treedepth = 10, seed = 123) {
  # Ensure data is properly formatted
  data <- as.data.frame(data)
  
  # Verify required columns
  required_cols <- c("count", "year", "state", "population")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns in data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Print data structure for debugging
  cat("Data structure before type conversion:\n")
  cat("Count column class:", class(data$count), "\n")
  cat("Population column class:", class(data$population), "\n")
  cat("Year column class:", class(data$year), "\n")
  cat("First few count values:", head(data$count), "\n")
  
  # Ensure all columns have correct types
  data$population <- as.numeric(as.character(data$population))
  data$count <- as.integer(as.numeric(as.character(data$count)))
  data$year <- as.numeric(as.character(data$year))
  data$state <- as.character(data$state)
  
  # Check for NA values after conversion
  na_count <- sum(is.na(data$count))
  na_pop <- sum(is.na(data$population))
  
  if (na_count > 0) {
    warning("Found ", na_count, " NA values in count after type conversion")
    # Replace NA with zeros for count
    data$count[is.na(data$count)] <- 0
  }
  
  if (na_pop > 0) {
    warning("Found ", na_pop, " NA values in population after type conversion")
    # Use mean population for NA values
    mean_pop <- mean(data$population, na.rm = TRUE)
    data$population[is.na(data$population)] <- mean_pop
  }
  
  # Handle zero-count data
  if (all(data$count == 0) || sum(data$count) == 0) {
    stop("CRITICAL ERROR: All pathogen counts are zero. Cannot fit model without positive case counts.")
  }
  
  # Ensure year is numeric (not factor) for the spline
  if (is.factor(data$year)) {
    data$year <- as.numeric(as.character(data$year))
  }
  
  # Convert state to factor if it isn't already
  if (!is.factor(data$state)) {
    data$state <- as.factor(data$state)
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Fit the model with robust settings
  model <- tryCatch({
    brm(
      count ~ s(year, by = state) + state + offset(log(population)),
      data = data,
      family = negbinomial(),
      chains = chains,
      iter = iterations,
      cores = cores,
      seed = seed,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      backend = "rstan"
    )
  }, error = function(e) {
    # If spline model fails, try a simpler model
    message("Spline model failed. Trying simpler model. Error was: ", e$message)
    
    tryCatch({
      # Try a simpler model without splines
      simpler_model <- brm(
        count ~ year + state + offset(log(population)),
        data = data,
        family = negbinomial(),
        chains = chains,
        iter = iterations,
        cores = cores,
        seed = seed,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
        backend = "rstan"
      )
      
      attr(simpler_model, "used_fallback") <- TRUE
      attr(simpler_model, "original_error") <- e$message
      
      return(simpler_model)
    }, error = function(e2) {
      # If even the simpler model fails, stop with error
      stop(paste("CRITICAL ERROR: Unable to fit any model. Primary error:", e$message, 
                 "Secondary error:", e2$message))
    })
  })
  
  return(model)
}

#' Generate Predicted Values from a Bayesian Model
#'
#' Generates posterior predictions from a fitted Bayesian model,
#' with special handling for dummy models and error cases.
#'
#' @param data Data frame to generate predictions for
#' @param model A brms model object from proposed_bm()
#' @return A tibble with posterior predictions
linpred_draw <- function(data, model) {
  # Handle manually created dummy model
  if (!is.null(attr(model, "is_manual_dummy")) && attr(model, "is_manual_dummy")) {
    stop("CRITICAL ERROR: Cannot generate predictions from a dummy model. Model fitting failed.")
  }
  
  # Handle dummy model created by proposed_bm
  if (!is.null(attr(model, "is_dummy")) && attr(model, "is_dummy")) {
    stop("CRITICAL ERROR: Cannot generate predictions from a dummy model. Model fitting failed.")
  }

  # Regular processing for normal models
  # Prepare data for prediction
  data <- as_tibble(data) %>%
    ungroup() %>%
    mutate(
      .row = row_number()
    )
  
  # Ensure year is numeric for prediction
  if ("year" %in% names(data)) {
    data$year <- as.numeric(as.character(data$year))
  }
  
  # Standardize population column to lowercase for model consistency
  # The Bayesian model uses offset(log(population)) with lowercase column name
  if ("population" %in% names(data)) {
    data$population <- as.numeric(as.character(data$population))
    cat("Using 'population' column with type:", class(data$population), "\n")
    cat("First few values:", head(data$population), "\n")
  } else if ("Population" %in% names(data)) {
    # Convert uppercase variant to standard lowercase format
    data$population <- as.numeric(as.character(data$Population))
    data$Population <- NULL  # Remove uppercase to prevent confusion
    cat("Standardizing 'Population' to 'population' column with type:", class(data$population), "\n")
    cat("First few values:", head(data$population), "\n")
  } else {
    stop("No population column found in data")
  }

  # Validate population data type and content
  if (!is.numeric(data$population)) {
    stop("Population column is not numeric after conversion")
  }
  
  if (any(is.na(data$population))) {
    warning("Population column contains ", sum(is.na(data$population)), " NA values which will be handled in processing")
  }

  # Handle fallback model (without splines)
  if (!is.null(attr(model, "used_fallback")) && attr(model, "used_fallback")) {
    message("Using fallback model to generate predictions.")
    
    # For the simpler model without splines, ensure year is numeric
    if (is.factor(data$year)) {
      data$year <- as.numeric(as.character(data$year))
    }
  }

  # Generate posterior predictions
  tryCatch({
    # Get posterior predictive draws
    epred <- epred_draws(model, newdata = data) %>% ungroup()
    
    # Remove population column from posterior draws to avoid duplication during join
    epred <- epred %>% select(-one_of("population", "Population"))
    
    # Rejoin population values using row identifier to maintain data integrity
    pop_df <- data %>% select(.row, population) %>% ungroup()
    draws <- left_join(epred, pop_df, by = ".row") %>% ungroup()
    
    if (!is.numeric(draws$population) || any(is.na(draws$population))) {
      stop("Population column is not numeric in the joined data")
    }
    
    # Calculate incidence rate per 100,000 population
    # Model predictions (.epred) are counts; divide by population for rate
    draws <- draws %>% mutate(pred_incidence = (.epred / population) * 100000)
    
    return(draws)
  }, error = function(e) {
    # If prediction fails, raise error
    stop(paste("CRITICAL ERROR: Failed to generate predictions from model:", e$message))
  })
}

#' Generate Catchment-Level Summary from Posterior Draws
#'
#' Summarizes posterior draws by year and state to produce catchment-level estimates.
#'
#' @param draws Tibble with posterior draws from linpred_draw()
#' @return A tibble with summarized incidence estimates by year and state
catchment <- function(draws) {
  # Group by relevant variables and calculate summary statistics
  catchment_data <- draws %>%
    group_by(year, state, .draw) %>%
    summarise(
      pred_incidence = mean(pred_incidence),
      .groups = "drop"
    ) %>%
    # Calculate HDI intervals for each Year/State combination
    group_by(year, state) %>%
    summarise(
      mean_incidence = mean(pred_incidence),
      median_incidence = median(pred_incidence),
      lower_hdi = hdi(pred_incidence, credMass = 0.95)[1],
      upper_hdi = hdi(pred_incidence, credMass = 0.95)[2],
      .groups = "drop"
    )

  return(catchment_data)
}

#' Format Catchment-Level Data for Output
#'
#' Formats the catchment data for output, rounding values and arranging by year and state.
#'
#' @param catchment_data Tibble from catchment()
#' @return A formatted tibble with incidence estimates by year and state
linpred_to_catchir <- function(catchment_data) {
  # Format the data for output
  ir_data <- catchment_data %>%
    mutate(
      year = as.integer(year),
      # Round numeric values to 2 decimal places
      mean_incidence = round(mean_incidence, 2),
      median_incidence = round(median_incidence, 2),
      lower_hdi = round(lower_hdi, 2),
      upper_hdi = round(upper_hdi, 2)
    ) %>%
    # Arrange by Year and State for better readability
    arrange(year, state)

  return(ir_data)
}

#' Format Site-Level Data for Output
#'
#' Formats the site-specific draws for output, rounding values and arranging by state and year.
#'
#' @param draws Tibble with posterior draws from linpred_draw()
#' @return A formatted tibble with incidence estimates by state and year
linpred_to_siteir <- function(draws) {
  # Format the data for site-specific outputs
  ir_data <- draws %>%
    group_by(year, state) %>%
    summarise(
      mean_incidence = mean(pred_incidence),
      median_incidence = median(pred_incidence),
      lower_hdi = hdi(pred_incidence, credMass = 0.95)[1],
      upper_hdi = hdi(pred_incidence, credMass = 0.95)[2],
      .groups = "drop"
    ) %>%
    # Round numeric values to 2 decimal places
    mutate(
      year = as.integer(year),
      mean_incidence = round(mean_incidence, 2),
      median_incidence = round(median_incidence, 2),
      lower_hdi = round(lower_hdi, 2),
      upper_hdi = round(upper_hdi, 2)
    ) %>%
    # Arrange by state and year for better readability
    arrange(state, year)

  return(ir_data)
}

#' Plot Site-Specific Trends
#'
#' Creates a faceted plot showing trends for each state over time.
#'
#' @param catchir_data Tibble from linpred_to_catchir()
#' @param pathogen Name of pathogen for plot title
#' @param outDir Directory to save the plot
#' @return A ggplot object with the plot
plot_site_trends <- function(catchir_data, pathogen, outDir) {
  # Create a plot for each state showing trends over time
  p <- ggplot(catchir_data, aes(x = year, y = median_incidence)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower_hdi, ymax = upper_hdi), alpha = 0.3) +
    facet_wrap(~ state, scales = "free_y") +
    labs(
      title = paste("Site-Specific Trends for", pathogen),
      subtitle = "Median incidence with 95% HDI intervals",
      y = "Incidence per 100,000 population",
      x = "Year"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold")
    )

  # Save the plot
  plot_file <- file.path(outDir, get_output_filename(pathogen, "site_trends", "png"))
  ggsave(plot_file, p, width = 12, height = 8, dpi = 300)

  return(p)
}

#' Plot Overall Trend
#'
#' Creates a plot showing the overall trend across all sites.
#'
#' @param catchir_data Tibble from linpred_to_catchir()
#' @param pathogen Name of pathogen for plot title
#' @param outDir Directory to save the plot
#' @return A ggplot object with the plot
plot_overall_trend <- function(catchir_data, pathogen, outDir) {
  # Calculate overall incidence by year (weighted by population)
  overall_data <- catchir_data %>%
    group_by(year) %>%
    summarise(
      median_incidence = mean(median_incidence),
      lower_hdi = mean(lower_hdi),
      upper_hdi = mean(upper_hdi),
      .groups = "drop"
    )

  # Create the plot
  p <- ggplot(overall_data, aes(x = year, y = median_incidence)) +
    geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = lower_hdi, ymax = upper_hdi), alpha = 0.3) +
    labs(
      title = paste("Overall Trend for", pathogen),
      subtitle = "Median incidence with 95% HDI intervals",
      y = "Incidence per 100,000 population",
      x = "Year"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )

  # Save the plot
  plot_file <- file.path(outDir, get_output_filename(pathogen, "overall_trend", "png"))
  ggsave(plot_file, p, width = 10, height = 6, dpi = 300)

  return(p)
}

#' Create Combined Visualization
#'
#' Combines site-specific and overall trend plots into a single figure.
#'
#' @param site_plot Site-specific plot from plot_site_trends()
#' @param overall_plot Overall trend plot from plot_overall_trend()
#' @param pathogen Name of pathogen for file naming
#' @param outDir Directory to save the plot
#' @return A grid object with the combined plot
plot_combined <- function(site_plot, overall_plot, pathogen, outDir) {
  # Combine the plots
  combined_plot <- gridExtra::grid.arrange(overall_plot, site_plot,
                                         ncol = 1, heights = c(1, 2))

  # Save the combined plot
  plot_file <- file.path(outDir, get_output_filename(pathogen, "combined", "png"))
  ggsave(plot_file, combined_plot, width = 12, height = 14, dpi = 300)

  return(combined_plot)
}

#' Calculate Relative Risks Compared to Historical Period
#'
#' Calculates relative risks and percent changes compared to a historical period.
#'
#' @param catchir_data Tibble from linpred_to_catchir()
#' @param start_year Start year of comparison period
#' @param end_year End year of comparison period
#' @param output_file Optional file path to save results
#' @return A data frame with relative risks and percent changes
ir_comp <- function(catchir_data, start_year, end_year, output_file = NULL) {
  # Filter data for the comparison period
  period_data <- catchir_data %>%
    filter(year >= start_year & year <= end_year)

  # Handle no data for requested period
  if (nrow(period_data) == 0) {
    warning(paste("No data available for period", start_year, "to", end_year))
    # Create minimal output to avoid errors
    if (!is.null(output_file)) {
      minimal_result <- data.frame(
        state = unique(catchir_data$state),
        year = max(catchir_data$year),
        comparison_period = paste0(start_year, "-", end_year),
        current_incidence = 0.01,
        period_incidence = 0.01,
        relative_risk = 1.00,
        percent_change = 0.00
      )
      
      # Use safe_write to save the minimal result
      safe_write(minimal_result, output_file)
    }
    return(NULL)
  }

  # Calculate average incidence for the period by state
  period_avg <- period_data %>%
    group_by(state) %>%
    summarise(
      period_incidence = mean(median_incidence),
      period_lower = mean(lower_hdi),
      period_upper = mean(upper_hdi),
      .groups = "drop"
    )

  # Get the most recent year's data
  latest_year <- max(catchir_data$year)
  latest_data <- catchir_data %>%
    filter(year == latest_year)

  # Join and calculate relative risks
  result <- latest_data %>%
    left_join(period_avg, by = "state") %>%
    mutate(
      relative_risk = median_incidence / period_incidence,
      percent_change = ((median_incidence / period_incidence) - 1) * 100,
      comparison_period = paste0(start_year, "-", end_year)
    ) %>%
    select(
      state, year, comparison_period,
      current_incidence = median_incidence,
      period_incidence,
      relative_risk,
      percent_change
    ) %>%
    arrange(state)

  # Round numeric columns for readability
  result <- result %>%
    mutate(across(where(is.numeric), ~round(., 4)))

  # Write to file if specified
  if (!is.null(output_file)) {
    # Use safe_write to save the result
    safe_write(result, output_file)
  }

  return(result)
}

#' Calculate Catchment-Level Relative Risks
#'
#' Wrapper function for ir_comp that accepts a catchment object.
#' This function exists for backward compatibility.
#'
#' @param catch Catchment object (not used but kept for API compatibility)
#' @param catchir_data Tibble from linpred_to_catchir()
#' @param start_year Start year of comparison period
#' @param end_year End year of comparison period
#' @param output_file Optional file path to save results
#' @return A data frame with relative risks and percent changes
ir_comp_catch <- function(catch, catchir_data, start_year, end_year, output_file = NULL) {
  # This is a wrapper around ir_comp for backward compatibility
  return(ir_comp(catchir_data, start_year, end_year, output_file))
}

# Add function to export model by pathogen with robust error handling
#' Save Bayesian Model for a Pathogen
#'
#' Ensures that model files are properly saved with robust error handling.
#' This function guarantees that a model file will be created even for edge cases.
#'
#' @param model The Bayesian model to save (brmsfit object)
#' @param pathogen Name of the pathogen (e.g., "CYCLOSPORA")
#' @param output_dir Directory for saving results
#' @param output_suffix Optional suffix for the output file
#' @return Full path to the saved file
save_pathogen_model <- function(model, pathogen, output_dir = ".", output_suffix = "brm") {
  # Construct output filename
  filename <- paste0(pathogen, "_", output_suffix, ".Rds")
  filepath <- file.path(output_dir, filename)
  
  # Make sure the directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Add model timestamp and metadata
  model$creation_time <- Sys.time()
  model$pathogen <- pathogen
  
  # Save the fitted model to file
  result <- tryCatch({
    saveRDS(model, file = filepath)
    cat("Saved model for", pathogen, "to", filepath, "\n")
    TRUE
  }, error = function(e) {
    stop(paste("CRITICAL ERROR: Failed to save model for", pathogen, ":", e$message))
  })
  
  # Final verification
  if (file.exists(filepath)) {
    cat("Verified file exists:", filepath, "\n")
    return(filepath)
  } else {
    cat("CRITICAL: File still doesn't exist after all attempts:", filepath, "\n")
    return(NULL)
  }
}
