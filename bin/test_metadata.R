#!/usr/bin/env Rscript
# Unit test for generate_metadata function with census file paths

# Source the calcIR.R file
source("bin/calcIR.R")

# Create a minimal test dataset
test_data <- data.frame(
  pathogen = c("SALMONELLA", "CAMPYLOBACTER", "SALMONELLA"),
  state = c("CA", "NY", "GA"),
  year = c(2020, 2021, 2022),
  county = c("LOS ANGELES", "NEW YORK", "FULTON"),
  serotypesummary = c("HEIDELBERG", "NA", "ENTERITIDIS"),
  stringsAsFactors = FALSE
)

# Test file paths
test_source_file <- "/path/to/mmwr/data.sas7bdat"
test_census_b <- "/path/to/census/bacterial.csv"
test_census_p <- "/path/to/census/parasitic.csv"

# Run without census file paths
cat("Test 1: Metadata without census file paths\n")
metadata1 <- generate_metadata(test_data, test_source_file)
cat("Census file fields in metadata:", 
    ifelse("census_file_bacterial" %in% names(metadata1), "YES", "NO"), 
    ifelse("census_file_parasitic" %in% names(metadata1), "YES", "NO"), 
    "\n\n")

# Run with census file paths
cat("Test 2: Metadata with census file paths\n")
metadata2 <- generate_metadata(test_data, test_source_file, test_census_b, test_census_p)
cat("Census file fields in metadata:", 
    ifelse("census_file_bacterial" %in% names(metadata2), "YES", "NO"), 
    ifelse("census_file_parasitic" %in% names(metadata2), "YES", "NO"), 
    "\n")
cat("Bacterial census path:", metadata2$census_file_bacterial, "\n")
cat("Parasitic census path:", metadata2$census_file_parasitic, "\n\n")

# Run with only bacterial census file path
cat("Test 3: Metadata with only bacterial census file path\n")
metadata3 <- generate_metadata(test_data, test_source_file, test_census_b)
cat("Census file fields in metadata:", 
    ifelse("census_file_bacterial" %in% names(metadata3), "YES", "NO"), 
    ifelse("census_file_parasitic" %in% names(metadata3), "YES", "NO"), 
    "\n")
cat("Bacterial census path:", metadata3$census_file_bacterial, "\n\n")

# Run with only parasitic census file path
cat("Test 4: Metadata with only parasitic census file path\n")
metadata4 <- generate_metadata(test_data, test_source_file, NULL, test_census_p)
cat("Census file fields in metadata:", 
    ifelse("census_file_bacterial" %in% names(metadata4), "YES", "NO"), 
    ifelse("census_file_parasitic" %in% names(metadata4), "YES", "NO"), 
    "\n")
cat("Parasitic census path:", metadata4$census_file_parasitic, "\n\n")

cat("All tests completed.\n")