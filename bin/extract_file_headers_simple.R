#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library("haven"))

# File paths (hard-coded for simplicity)
mmwrFile <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/mmwr9623_Jan2024.sas7bdat"
censusFile_B <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623.sas7bdat"
censusFile_P <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623_para.sas7bdat"

# Output file
output_file <- "file_headers_output.txt"
cat("Saving output to:", output_file, "\n")

# Helper function to inspect data
inspect_data <- function(file_path, file_label) {
  cat("\n--- START OF:", file_label, "---\n")
  
  data <- tryCatch(
    {
      haven::read_sas(file_path)
    },
    error = function(e) {
      cat("ERROR reading file:", file_path, "\n")
      return(NULL)
    }
  )
  
  if (!is.null(data)) {
    cat("\nFile loaded successfully: ", file_path, "\n")
    cat("Number of Rows:", nrow(data), "\n")
    cat("Number of Columns:", ncol(data), "\n\n")
    cat("Column Names:\n")
    print(colnames(data))
    cat("\nFirst 5 Rows:\n")
    print(head(data, 5))
  }
  
  cat("\n--- END OF:", file_label, "---\n\n")
}

# Create or overwrite the output file
sink(output_file)

# Inspect each file
inspect_data(mmwrFile, "MMWR File")
inspect_data(censusFile_B, "Census Bacterial File")
inspect_data(censusFile_P, "Census Parasitic File")

sink() # Close the output file

cat("Extraction complete. Check the file:", output_file, "\n")

