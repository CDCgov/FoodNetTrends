#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("haven"))

# Set up argument parsing
parser <- ArgumentParser()
parser$add_argument("--sasFile", type = "character", help = "Path to the SAS input file", required = TRUE)
parser$add_argument("--nRows", type = "integer", default = 10, help = "Number of rows to display (default: 10)")
parser$add_argument("--outFile", type = "character", default = "input_snippet.txt", help = "Path to output text file (default: input_snippet.txt)")
args <- parser$parse_args()

# Read the SAS file
data <- read_sas(args$sasFile)

# Capture the structure, column names, and first few rows
snippet <- capture.output({
  cat("----- Structure of the Data -----\n")
  str(data)
  cat("\n----- Column Names -----\n")
  print(colnames(data))
  cat("\n----- First", args$nRows, "Rows -----\n")
  print(head(data, args$nRows))
})

# Write the captured output to the specified file
writeLines(snippet, con = args$outFile)

cat("Snippet saved to", args$outFile, "\n")

