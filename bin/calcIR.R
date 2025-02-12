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
# Usage:
#   Rscript calcIR.R --mmwrFile <path_to_raw_SAS_file> --outputFile <path_to_output_csv>
#
# Example:
#   Rscript calcIR.R --mmwrFile "/path/to/mmwr9623_Jan2024.sas7bdat" --outputFile "clean_mmwr.csv"
#
################################################################################

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("haven"))

# Setup argument parser
parser <- ArgumentParser()
parser$add_argument("--mmwrFile", type = "character", help = "Path to the raw MMWR SAS file", required = TRUE)
parser$add_argument("--outputFile", type = "character", help = "Path to save the cleaned CSV file", required = TRUE)
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
# Insert any additional data cleaning or filtering here if needed.

# --- Standardize and Rename Key Columns ---
# Ensure the raw data has the necessary columns and then rename them:
# If the cleaned file still has lowercase names (e.g., 'pathogen', 'state', 'year'),
# we explicitly rename them to the expected format.
mmwrdata <- mmwrdata %>%
  rename(
    Pathogen = pathogen,
    State    = state,
    Year     = year
  )

# --- Create or Verify Derived Columns ---
# Create a derived column 'pathogentype' if not already present.
# We assume that if Pathogen is one of "CRYPTOSPORIDIUM" or "CYCLOSPORA", it is "Parasitic"; otherwise "Bacterial".
if(!"pathogentype" %in% names(mmwrdata)) {
  mmwrdata <- mmwrdata %>%
    mutate(pathogentype = ifelse(Pathogen %in% c("CRYPTOSPORIDIUM", "CYCLOSPORA"), "Parasitic", "Bacterial"))
}

# --- Write Cleaned Data to CSV ---
cat("Writing cleaned data to:", args$outputFile, "\n")
write.csv(mmwrdata, file = args$outputFile, row.names = FALSE)
cat("Data cleaning complete. Cleaned data saved to:", args$outputFile, "\n")

