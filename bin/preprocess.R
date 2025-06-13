#!/usr/bin/env Rscript
################################################################################
# preprocess.R
#
# Purpose:
#   This script cleans and aggregates raw MMWR SAS data and writes out a CSV file
#   with standardized column names required for downstream analysis.
#
#   The cleaning includes:
#     - Reading the raw SAS file.
#     - Converting all column names to lowercase.
#     - Recoding SERO values: creating a new column 'sero2' (and copying it into
#       'serotypesummary') so that values like "NOT SPECIATED", "UNKNOWN", etc.
#       are recoded to "Missing".
#     - Standardizing county names.
#     - (Any other cleaning steps can be added here as needed.)
#
#   Finally, key columns are renamed so that:
#     - 'pathogen' becomes 'Pathogen'
#     - 'state' becomes 'State'
#     - 'year' becomes 'Year'
#     - A derived 'pathogentype' column is created (if not already present) to 
#       distinguish between "Parasitic" and "Bacterial" pathogens.
#
# Usage:
#   Rscript preprocess.R --mmwrFile <path_to_raw_SAS_file> --outputFile <path_to_output_csv>
#
# Example:
#   Rscript preprocess.R --mmwrFile "/path/to/mmwr9623_Jan2024.sas7bdat" --outputFile "clean_mmwr.csv"
#
################################################################################

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("haven"))
suppressPackageStartupMessages(library("gtools"))
# stringdist is optional - used for fuzzy pathogen matching if available
if (requireNamespace("stringdist", quietly = TRUE)) {
  suppressPackageStartupMessages(library("stringdist"))
}

# Setup argument parser
parser <- ArgumentParser()
parser$add_argument("--mmwrFile", type = "character", help = "Path to the raw MMWR SAS file", required = TRUE)
parser$add_argument("--outputFile", type = "character", help = "Path to save the cleaned CSV file", required = TRUE)
parser$add_argument("--serotype-config", type = "character", help = "Path to CSV file with serotype recoding rules (optional)", required = FALSE, default = NULL)
parser$add_argument("--matching-sensitivity", type = "character", help = "Pathogen matching sensitivity: STRICT, MEDIUM, RELAXED (default: MEDIUM)", 
                    required = FALSE, default = "MEDIUM", choices = c("STRICT", "MEDIUM", "RELAXED"))
args <- parser$parse_args()

# --- Serotype Configuration Functions ---
read_serotype_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    # Return default configuration as data frame
    return(data.frame(
      serotype_value = c("NOT SPECIATED", "UNKNOWN", "PARTIAL SERO", "NOT SERO", "", "UNDET"),
      replacement = c("Missing", "Missing", "Missing", "Missing", "Missing", "Missing"),
      pathogen = rep("all", 6),
      match_type = c("exact", "exact", "exact", "exact", "exact", "contains"),
      stringsAsFactors = FALSE
    ))
  }
  # Read and validate CSV
  config <- read.csv(config_path, stringsAsFactors = FALSE)
  validate_serotype_config(config)
  return(config)
}

validate_serotype_config <- function(config) {
  required_cols <- c("serotype_value", "replacement", "pathogen", "match_type")
  missing_cols <- setdiff(required_cols, names(config))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in serotype config: ", paste(missing_cols, collapse = ", "))
  }
  # Validate match_type values
  valid_match_types <- c("exact", "contains")
  invalid_types <- setdiff(unique(config$match_type), valid_match_types)
  if (length(invalid_types) > 0) {
    stop("Invalid match_type values: ", paste(invalid_types, collapse = ", "))
  }
}

apply_serotype_config <- function(data, config) {
  # Initialize sero2 column with original values
  data$sero2 <- data$sero1
  
  # Apply exact matches first
  exact_rules <- config[config$match_type == "exact", ]
  for (i in seq_len(nrow(exact_rules))) {
    if (exact_rules$pathogen[i] == "all") {
      # Apply to all pathogens
      data$sero2[data$sero1 == exact_rules$serotype_value[i]] <- exact_rules$replacement[i]
    } else {
      # Apply to specific pathogen
      mask <- (data$sero1 == exact_rules$serotype_value[i]) & (data$pathogen == exact_rules$pathogen[i])
      data$sero2[mask] <- exact_rules$replacement[i]
    }
  }
  
  # Apply pattern matches
  pattern_rules <- config[config$match_type == "contains", ]
  for (i in seq_len(nrow(pattern_rules))) {
    if (pattern_rules$pathogen[i] == "all") {
      # Apply to all pathogens
      data$sero2[grepl(pattern_rules$serotype_value[i], data$sero2)] <- pattern_rules$replacement[i]
    } else {
      # Apply to specific pathogen
      mask <- grepl(pattern_rules$serotype_value[i], data$sero2) & (data$pathogen == pattern_rules$pathogen[i])
      data$sero2[mask] <- pattern_rules$replacement[i]
    }
  }
  
  return(data)
}

# --- Pathogen Standardization Functions ---
# Define known pathogen patterns and matching rules
# These patterns handle common variations in pathogen naming found in surveillance data:
# - Abbreviations (S. for Salmonella, C. for Campylobacter)
# - Misspellings (SHIGELA, SALMONEL)
# - Partial names (CAMP for Campylobacter)
# - Historical naming variations
get_pathogen_patterns <- function() {
  list(
    SALMONELLA = list(
      exact = c("SALMONELLA"),
      prefix = c("^SAL", "^S\\.", "^S "),  # Handles: SAL*, S.enteritidis, S typhimurium
      contains = c("SALMONELLA", "SALMONEL")
    ),
    CAMPYLOBACTER = list(
      exact = c("CAMPYLOBACTER"),
      prefix = c("^CAMP", "^C\\.", "^C "),
      contains = c("CAMPYLOBACTER", "CAMPY")
    ),
    SHIGELLA = list(
      exact = c("SHIGELLA"),
      prefix = c("^SHIG"),
      contains = c("SHIGELL", "SHIGELA")
    ),
    LISTERIA = list(
      exact = c("LISTERIA"),
      prefix = c("^LIST", "^L\\.", "^L "),
      contains = c("LISTERIA")
    ),
    STEC = list(
      exact = c("STEC", "E. COLI O157:H7"),
      prefix = c("^STEC", "^E.*COLI", "^ECOLI"),
      contains = c("E.*COLI.*O157", "STEC")
    ),
    VIBRIO = list(
      exact = c("VIBRIO"),
      prefix = c("^VIB", "^V\\.", "^V "),
      contains = c("VIBRIO")
    ),
    YERSINIA = list(
      exact = c("YERSINIA"),
      prefix = c("^YERS", "^Y\\.", "^Y "),
      contains = c("YERSINIA")
    ),
    CYCLOSPORA = list(
      exact = c("CYCLOSPORA"),
      prefix = c("^CYCL", "^CYC"),
      contains = c("CYCLOSPORA", "CYCLO")
    ),
    CRYPTOSPORIDIUM = list(
      exact = c("CRYPTOSPORIDIUM"),
      prefix = c("^CRYP"),
      contains = c("CRYPTOSPORID", "CRYPTO")
    )
  )
}

# Define matching sensitivity levels
# STRICT: Only exact matches - use when data quality is high and pathogen names are standardized
# MEDIUM: Exact + prefix + limited fuzzy - balanced approach for typical surveillance data
# RELAXED: All methods + broader fuzzy - use when data has many variations or historical inconsistencies
get_matching_settings <- function(sensitivity) {
  settings <- list(
    STRICT = list(
      use_exact = TRUE,
      use_prefix = FALSE,
      use_contains = FALSE,
      use_fuzzy = FALSE,
      fuzzy_distance = 0,
      confidence_threshold = 95  # Only accept very high confidence matches
    ),
    MEDIUM = list(
      use_exact = TRUE,
      use_prefix = TRUE,
      use_contains = FALSE,
      use_fuzzy = TRUE,
      fuzzy_distance = 1,
      confidence_threshold = 80
    ),
    RELAXED = list(
      use_exact = TRUE,
      use_prefix = TRUE,
      use_contains = TRUE,
      use_fuzzy = TRUE,
      fuzzy_distance = 2,
      confidence_threshold = 70
    )
  )
  return(settings[[sensitivity]])
}

# Standardize pathogen names based on patterns and sensitivity
standardize_pathogens <- function(data, sensitivity = "MEDIUM") {
  cat("Standardizing pathogen names with", sensitivity, "sensitivity\n")
  
  settings <- get_matching_settings(sensitivity)
  patterns <- get_pathogen_patterns()
  
  # Track original values for reporting
  data$pathogen_original <- data$pathogen
  
  # Initialize tracking for report
  standardization_log <- list()
  
  # Process each unique pathogen value
  unique_pathogens <- unique(data$pathogen)
  
  for (orig_value in unique_pathogens) {
    orig_upper <- toupper(trimws(orig_value))
    matched <- FALSE
    match_info <- list(
      original = orig_value,
      standardized = NA,
      match_type = "no_match",
      confidence = 0,
      count = sum(data$pathogen == orig_value)
    )
    
    # Try exact matches first
    if (settings$use_exact && !matched) {
      for (pathogen in names(patterns)) {
        if (orig_upper %in% patterns[[pathogen]]$exact) {
          data$pathogen[data$pathogen == orig_value] <- pathogen
          match_info$standardized <- pathogen
          match_info$match_type <- "exact"
          match_info$confidence <- 100
          matched <- TRUE
          break
        }
      }
    }
    
    # Try prefix matches
    if (settings$use_prefix && !matched) {
      for (pathogen in names(patterns)) {
        for (pattern in patterns[[pathogen]]$prefix) {
          if (grepl(pattern, orig_upper)) {
            data$pathogen[data$pathogen == orig_value] <- pathogen
            match_info$standardized <- pathogen
            match_info$match_type <- "prefix"
            match_info$confidence <- 90
            matched <- TRUE
            break
          }
        }
        if (matched) break
      }
    }
    
    # Try contains matches
    if (settings$use_contains && !matched) {
      for (pathogen in names(patterns)) {
        for (pattern in patterns[[pathogen]]$contains) {
          if (grepl(pattern, orig_upper)) {
            data$pathogen[data$pathogen == orig_value] <- pathogen
            match_info$standardized <- pathogen
            match_info$match_type <- "contains"
            match_info$confidence <- 80
            matched <- TRUE
            break
          }
        }
        if (matched) break
      }
    }
    
    # Try fuzzy matching
    if (settings$use_fuzzy && !matched) {
      # Load stringdist if available
      if (requireNamespace("stringdist", quietly = TRUE)) {
        known_pathogens <- names(patterns)
        distances <- stringdist::stringdist(orig_upper, known_pathogens, method = "lv")
        min_dist <- min(distances)
        
        if (min_dist <= settings$fuzzy_distance) {
          best_match <- known_pathogens[which.min(distances)]
          data$pathogen[data$pathogen == orig_value] <- best_match
          match_info$standardized <- best_match
          match_info$match_type <- paste0("fuzzy_d", min_dist)
          match_info$confidence <- 100 - (min_dist * 10)
          matched <- TRUE
        }
      }
    }
    
    # Record unmatched
    if (!matched) {
      match_info$standardized <- orig_value  # Keep original if no match
      match_info$match_type <- "no_match"
      match_info$confidence <- 0
    }
    
    standardization_log[[orig_value]] <- match_info
  }
  
  # Generate report
  report <- do.call(rbind, lapply(standardization_log, as.data.frame, stringsAsFactors = FALSE))
  report <- report[order(report$standardized, -report$count), ]
  
  return(list(data = data, report = report))
}

# --- Data Loading ---
cat("Loading raw MMWR data from:", args$mmwrFile, "\n")
mmwrdata <- haven::read_sas(args$mmwrFile) %>% 
  as.data.frame()

# Convert all column names to lowercase for consistency FIRST
mmwrdata <- mmwrdata %>% rename_all(tolower)

# Now filter using lowercase column name
mmwrdata <- mmwrdata %>% 
  filter(siteid != "COEX")  # Exclude Colorado expanded catchment

# Column names already converted to lowercase above

# --- Pathogen Standardization ---
# Standardize pathogen names before any other processing
standardization_result <- standardize_pathogens(mmwrdata, sensitivity = args$`matching-sensitivity`)
mmwrdata <- standardization_result$data

# Write standardization report
report_file <- sub("\\.csv$", "_preprocessing_report.csv", args$outputFile)
write.csv(standardization_result$report, file = report_file, row.names = FALSE)
cat("Pathogen standardization report saved to:", report_file, "\n")

# --- Data Cleaning: Recoding SERO Variables ---
# Load serotype configuration (use defaults if no config file provided)
serotype_config <- read_serotype_config(args$`serotype-config`)
if (!is.null(args$`serotype-config`)) {
  cat("Using serotype configuration from:", args$`serotype-config`, "\n")
} else {
  cat("Using default serotype configuration\n")
}

# Apply serotype recoding based on configuration
mmwrdata <- apply_serotype_config(mmwrdata, serotype_config)

# Copy sero2 to serotypesummary (the column we want to preserve downstream)
mmwrdata$serotypesummary <- mmwrdata$sero2

# --- Data Cleaning: Standardize County Names ---
# Correct common issues in county names.
mmwrdata <- mmwrdata %>%
  mutate(
    county = if_else(county %in% c("ST. MARYS'S", "ST. MARYS"), "ST. MARY'S", county),
    county = if_else(county == "PRINCE GEORGES", "PRINCE GEORGE'S", county),
    county = if_else(county == "QUEEN ANNES", "QUEEN ANNE'S", county),
    county = if_else(county == "DE BACA", "DEBACA", county)
  )

# Process special pathogen cases

# STEC processing
# STEC (Shiga toxin-producing E. coli) default handling:
# - O157 and non-O157 strains have different epidemiological patterns
# - CDC tracks them separately for public health surveillance
# - The stec_class field identifies the serogroup
stec_data <- NULL
if("STEC" %in% unique(mmwrdata$pathogen)) {
  # Check if stec_class column exists
  if("stec_class" %in% names(mmwrdata)) {
    # For preprocessing, we keep STEC data intact
    # The actual splitting will happen during analysis based on user preference
    cat("STEC data found with stec_class information available for grouping.\n")
  } else {
    cat("Warning: STEC found but stec_class column not present. STEC grouping options will be limited.\n")
  }
}

# Listeria CSTE filtering
# Listeria monocytogenes is filtered by CSTE (Council of State and Territorial Epidemiologists) 
# case definition to ensure only invasive listeriosis cases are included
# Non-invasive cases (e.g., gastroenteritis) are excluded for consistency with national reporting
if("LISTERIA" %in% unique(mmwrdata$pathogen)) {
  # Check if cste column exists
  if("cste" %in% names(mmwrdata)) {
    # Filter to only CSTE-reportable cases
    mmwrdata <- mmwrdata %>% 
      filter(!(pathogen == "LISTERIA" & cste != "YES"))
    cat("Filtered Listeria cases to CSTE-reportable only.\n")
  } else {
    cat("Warning: LISTERIA found but cste column not present. Including all LISTERIA cases.\n")
  }
}

# --- (Optional) Additional Cleaning Steps ---
# Insert any additional data cleaning or filtering here if needed.

# --- Standardize and Rename Key Columns ---
# IMPORTANT: Keep all column names lowercase for consistency throughout the pipeline
# The downstream analysis functions (PATH_ANALYSIS, etc.) expect lowercase column names
# The commented line below shows incorrect capitalization that would break compatibility:
# mmwrdata <- mmwrdata %>% rename(Pathogen = pathogen, State = state, Year = year)

# --- Create or Verify Derived Columns ---
# Create a derived column 'pathogentype' if not already present.
# Pathogen classification: CRYPTOSPORIDIUM and CYCLOSPORA are parasitic pathogens; all others are bacterial
if(!"pathogentype" %in% names(mmwrdata)) {
  mmwrdata <- mmwrdata %>%
    mutate(pathogentype = ifelse(pathogen %in% c("CRYPTOSPORIDIUM", "CYCLOSPORA"), "Parasitic", "Bacterial"))
}

# --- Write Cleaned Data to CSV ---
cat("Writing cleaned data to:", args$outputFile, "\n")
write.csv(mmwrdata, file = args$outputFile, row.names = FALSE)
cat("Data cleaning complete. Cleaned data saved to:", args$outputFile, "\n")

