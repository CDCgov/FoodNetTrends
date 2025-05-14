#!/usr/bin/env Rscript

# Script to discover data in MMWR file for FoodNet Trends
suppressPackageStartupMessages(library(haven))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(jsonlite))

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: discover_data.R <mmwr_file> [output_file]")
}

mmwr_file <- args[1]
output_file <- if (length(args) >= 2) args[2] else "data_inventory.json"

cat("Reading MMWR data file:", mmwr_file, "\n")

# Read the MMWR data file
data <- tryCatch({
    haven::read_sas(mmwr_file) %>% as.data.frame()
}, error = function(e) {
    stop(paste("Error reading MMWR file:", e$message))
})

cat("Successfully read", nrow(data), "records\n")

# Convert column names to lowercase for consistency
names(data) <- tolower(names(data))

# Extract unique pathogens, states, and years
inventory <- list(
  pathogens = sort(unique(toupper(data$pathogen))),
  states = sort(unique(toupper(data$state))),
  years = sort(unique(as.numeric(as.character(data$year))))
)

# Calculate some basic statistics
inventory$counts <- list(
  total_records = nrow(data),
  pathogen_counts = as.list(table(toupper(data$pathogen))),
  state_counts = as.list(table(toupper(data$state)))
)

# Check for Salmonella serotypes if they exist
if (any(toupper(data$pathogen) == "SALMONELLA")) {
    # Try different serotype columns that might exist
    serotype_col <- NULL
    if ("serotypesummary" %in% names(data)) {
        serotype_col <- "serotypesummary"
    } else if ("sero1" %in% names(data)) {
        serotype_col <- "sero1"
    } else if ("sero2" %in% names(data)) {
        serotype_col <- "sero2"
    }
    
    if (!is.null(serotype_col)) {
        sal_data <- data[toupper(data$pathogen) == "SALMONELLA", ]
        serotypes <- sal_data[[serotype_col]]
        
        # Count occurrences of each serotype
        serotype_counts <- as.data.frame(table(serotypes))
        serotype_counts <- serotype_counts[order(serotype_counts$Freq, decreasing=TRUE),]
        
        # Store in inventory
        inventory$salmonella_serotypes <- as.list(serotype_counts$Freq)
        names(inventory$salmonella_serotypes) <- serotype_counts$serotypes
        
        # Also store a flat list for easier usage
        inventory$salmonella_serotype_names <- as.character(serotype_counts$serotypes)
        
        cat("Found", length(inventory$salmonella_serotype_names), "Salmonella serotypes\n")
    }
}

# Add timestamp
inventory$generated <- as.character(Sys.time())
inventory$source_file <- mmwr_file

# Write the inventory
write_json(inventory, output_file, pretty = TRUE)

# Print quick summary
cat("Discovered", length(inventory$pathogens), "pathogens\n")
cat("Discovered", length(inventory$states), "states\n")
cat("Discovered", length(inventory$years), "years\n")
cat("Total records:", inventory$counts$total_records, "\n")
cat("Inventory saved to", output_file, "\n")
