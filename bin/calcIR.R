#!/usr/bin/env Rscript
################################################################################
# calcIR.R
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
#   NEW: Optionally generates a metadata JSON file with information about the 
#   data contents (pathogens, states, serotypes, etc.).
#
# Usage:
#   Rscript calcIR.R --mmwrFile <path_to_raw_SAS_file> --outputFile <path_to_output_csv> [--generate_metadata true/false]
#
# Example:
#   Rscript calcIR.R --mmwrFile "/path/to/mmwr9623_Jan2024.sas7bdat" --outputFile "clean_mmwr.csv"
#
################################################################################

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("haven"))
# Add jsonlite for metadata generation
suppressPackageStartupMessages(library("jsonlite"))

# Setup argument parser
parser <- ArgumentParser()
parser$add_argument("--mmwrFile", type = "character", help = "Path to the raw MMWR SAS file", required = TRUE)
parser$add_argument("--outputFile", type = "character", help = "Path to save the cleaned CSV file", required = TRUE)
# Add new parameter for metadata generation
parser$add_argument("--generate_metadata", type = "logical", default = FALSE, 
                   help = "Whether to generate a metadata JSON file (default: FALSE)")
args <- parser$parse_args()

# --- Data Loading ---
cat("Loading raw MMWR data from:", args$mmwrFile, "\n")
mmwrdata <- haven::read_sas(args$mmwrFile) %>% as.data.frame()

# --- Standardize Column Names ---
# Convert all column names to lowercase for consistency.
mmwrdata <- mmwrdata %>% rename_all(tolower)

# --- Data Cleaning: Recoding SERO Variables ---
# Define a list of SERO values to be considered non-informative.
seroList <- c("NOT SPECIATED", "UNKNOWN", "PARTIAL SERO", "NOT SERO", "")
# Create a new column 'sero2': recode values in SERO1 that are in seroList as "Missing"
mmwrdata$sero2 <- ifelse(mmwrdata$sero1 %in% seroList, "Missing", mmwrdata$sero1)
# Further, if 'sero2' contains the string "UNDET", recode it to "Missing"
mmwrdata$sero2 <- ifelse(grepl("UNDET", mmwrdata$sero2), "Missing", mmwrdata$sero2)
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

# --- (Optional) Additional Cleaning Steps ---
# Ensure pathogen column is uppercase for consistency
mmwrdata$pathogen <- toupper(mmwrdata$pathogen)

# --- Standardize and Rename Key Columns ---
# Ensure the raw data has the necessary columns and then rename them:
# If the cleaned file still has lowercase names (e.g., 'pathogen', 'state', 'year'),
# we explicitly rename them to the expected format.
#mmwrdata <- mmwrdata %>%
#  rename(
#    Pathogen = pathogen,
#    State    = state,
#    Year     = year
#  )

# --- Create or Verify Derived Columns ---
# Create a derived column 'pathogentype' if not already present.
# We assume that if Pathogen is one of "CRYPTOSPORIDIUM" or "CYCLOSPORA", it is "Parasitic"; otherwise "Bacterial".
if(!"pathogentype" %in% names(mmwrdata)) {
  mmwrdata <- mmwrdata %>%
    mutate(pathogentype = ifelse(pathogen %in% c("CRYPTOSPORIDIUM", "CYCLOSPORA"), "Parasitic", "Bacterial"))
}

# --- Generate Metadata if Requested ---
if(args$generate_metadata) {
  cat("Generating metadata from cleaned data...\n")
  
  # Determine output path for metadata JSON (same base name as CSV but with _metadata.json extension)
  metadata_file <- sub("\\.csv$", "_metadata.json", args$outputFile)
  
  # Extract key information for metadata
  metadata <- list(
    pathogens = sort(unique(mmwrdata$pathogen)),
    states = sort(unique(mmwrdata$state)),
    years = sort(unique(as.numeric(as.character(mmwrdata$year)))),
    counties = sort(unique(mmwrdata$county)),
    generated_timestamp = as.character(Sys.time()),
    source_file = args$mmwrFile,
    record_count = nrow(mmwrdata)
  )
  
  # Add counts for basic statistics
  metadata$counts <- list(
    total_records = nrow(mmwrdata),
    pathogen_counts = as.list(table(mmwrdata$pathogen)),
    state_counts = as.list(table(mmwrdata$state))
  )
  
  # Process Salmonella serotypes if available
  if (any(mmwrdata$pathogen == "SALMONELLA") && "serotypesummary" %in% names(mmwrdata)) {
    sal_data <- mmwrdata[mmwrdata$pathogen == "SALMONELLA", ]
    serotype_counts <- as.data.frame(table(sal_data$serotypesummary))
    serotype_counts <- serotype_counts[order(serotype_counts$Freq, decreasing=TRUE),]
    
    # Store in metadata
    metadata$salmonella_serotypes <- as.list(serotype_counts$Freq)
    names(metadata$salmonella_serotypes) <- serotype_counts$Var1
    
    # Also store a flat list of names
    metadata$salmonella_serotype_names <- as.character(serotype_counts$Var1)
    
    cat("Found", length(metadata$salmonella_serotype_names), "Salmonella serotypes\n")
  }
  
  # Write metadata JSON
  cat("Writing metadata to:", metadata_file, "\n")
  write_json(metadata, metadata_file, pretty = TRUE)
}

# --- Write Cleaned Data to CSV ---
cat("Writing cleaned data to:", args$outputFile, "\n")
write.csv(mmwrdata, file = args$outputFile, row.names = FALSE)
cat("Data cleaning complete. Cleaned data saved to:", args$outputFile, "\n")

if(args$generate_metadata) {
  metadata_file <- sub("\\.csv$", "_metadata.json", args$outputFile)
  cat("Metadata saved to:", metadata_file, "\n")
}
