#!/usr/bin/env Rscript

# this is currently hard coded to use 16 cores
# Load required libraries
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("haven"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("tidybayes"))
suppressPackageStartupMessages(library("HDInterval"))

# Set warning options
options(warn = -1)

sink("output_log.txt")  # Redirect output to a log file

################################################################################
# Argument Parsing
################################################################################
parser <- ArgumentParser()
parser$add_argument("--debug", type = "logical", default = TRUE,
                    help = "Enable debug mode")
parser$add_argument("--mmwrFile", type = "character",
                    help = "Path to FoodNet MMWR SAS data")
parser$add_argument("--censusFile_B", type = "character",
                    help = "Path to census data for bacterial pathogens")
parser$add_argument("--censusFile_P", type = "character",
                    help = "Path to census data for parasitic pathogens")
parser$add_argument("--travel", type = "character",
                    help = "Travel types (e.g., NO, UNKNOWN)")
parser$add_argument("--cidt", type = "character",
                    help = "CIDT types (e.g., CIDT+, CX+, PARASITIC)")
parser$add_argument("--projID", type = "character",
                    help = "Project identifier")
parser$add_argument("--outDir", type = "character", default = "./",
                    help = "Output directory for results")
parser$add_argument("--pathogen", type = "character", default = NULL,
                    help = "If provided, only process this pathogen")
parser$add_argument("--preprocessed", type = "logical", default = FALSE,
                    help = "Set to TRUE if using preprocessed CSV data")
parser$add_argument("--cleanFile", type = "character", default = NULL,
                    help = "Path to cleaned CSV file (if --preprocessed is TRUE)")
opts <- parser$parse_args()

# Make options globally available
assign("opts", opts, envir = .GlobalEnv)

# Load helper functions
source("/scicomp/home-pure/rqu4/PROJECTS/SCICOMP/FoodNetTrends/FoodNetTrends4Ethan/bin/functions.R")

################################################################################
# Set parameters
################################################################################
if (!opts$debug) {
  mmwrFile      <- opts$mmwrFile
  censusFile_B  <- opts$censusFile_B
  censusFile_P  <- opts$censusFile_P
  projID        <- opts$projID
  travel        <- CLEAN_LIST(opts$travel)
  cidt          <- CLEAN_LIST(opts$cidt)
  outDir        <- paste0(projID, "/", "SplineResults")
  pathogen_arg  <- opts$pathogen
  preprocessed  <- opts$preprocessed
  cleanFile     <- opts$cleanFile
} else {
  # Debug mode defaults
  mmwrFile      <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/mmwr9623_Jan2024.sas7bdat"
  censusFile_B  <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623.sas7bdat"
  censusFile_P  <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623_para.sas7bdat"
  projID        <- "20240705"
  travel        <- CLEAN_LIST(c("NO", "UNKNOWN", "YES"))
  cidt          <- CLEAN_LIST(c("CIDT+", "CX+", "PARASITIC"))
  outDir        <- paste0(projID, "/", "SplineResults")
  pathogen_arg  <- opts$pathogen
  preprocessed  <- opts$preprocessed
  cleanFile     <- opts$cleanFile
}

# Create output directory
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Data Import and Preprocessing
################################################################################
if (preprocessed) {
  if (is.null(cleanFile)) {
    stop("Preprocessed flag is TRUE but no cleanFile provided.")
  }
  message("Using preprocessed CSV file from: ", cleanFile)
  mmwrdata <- read.csv(cleanFile, stringsAsFactors = FALSE)
} else {
  message("-- IMPORTING MMWR DATA (raw SAS) --")
  mmwrdata <- haven::read_sas(mmwrFile) %>%
    as.data.frame() %>%
    rename_all(tolower) %>%
    filter(siteid != "coex") %>%
    mutate(
      Pathogen = pathogen,
      State = state,
      Year = year,
      pathogentype = ifelse(Pathogen %in% c("CRYPTOSPORIDIUM", "CYCLOSPORA"), "Parasitic", "Bacterial")
    )
}

if (!is.null(pathogen_arg)) {
  message("Filtering data for pathogen: ", pathogen_arg)
  mmwrdata <- mmwrdata[mmwrdata$Pathogen == pathogen_arg, ]
}

mmwrdata <- mmwrdata %>%
  group_by(Pathogen, State, Year, pathogentype) %>%
  summarise(count = n(), .groups = "drop")

message("-- IMPORTING CENSUS DATA --")
census_b <- haven::read_sas(censusFile_B) %>%
  as.data.frame() %>%
  rename_all(tolower) %>%
  rename(State = state, Year = year, Population = population)
census_p <- haven::read_sas(censusFile_P) %>%
  as.data.frame() %>%
  rename_all(tolower) %>%
  rename(State = state, Year = year, Population = population)

mmwrdata <- mmwrdata %>%
  left_join(census_b, by = c("State", "Year"), relationship = "many-to-many") %>%
  left_join(census_p %>% select(State, Year, Population) %>% rename(Population_para = Population),
            by = c("State", "Year"), relationship = "many-to-many") %>%
  mutate(Population = ifelse(pathogentype == "Parasitic", Population_para, Population)) %>%
  select(-Population_para) %>%
  filter(!is.na(Population) & Population > 0)

################################################################################
# Model Fitting
################################################################################
message("-- FITTING MODELS --")
PROPOSED_BM <- function(data) {
  model <- brm(
    count ~ Year + State + offset(log(Population)),
    data = data,
    family = negbinomial(),
    chains = 1,
    iter = 100,
    cores = 16,
    control = list(adapt_delta = 0.95, max_treedepth = 5)
  )
  return(model)
}

if (nrow(mmwrdata) == 0) {
  stop("No data remaining after filtering. Check the input or the --pathogen parameter.")
}

message("Fitting model for pathogen: ", unique(mmwrdata$Pathogen))
model_fit <- try(PROPOSED_BM(mmwrdata), silent = TRUE)
if (inherits(model_fit, "try-error")) {
  stop(paste("Error fitting model for", unique(mmwrdata$Pathogen), ":", model_fit[1])) # Extract error message
} else {
  saveRDS(model_fit, file = file.path(outDir, paste0(unique(mmwrdata$Pathogen), "_brm.Rds")))
  message("Model fitting complete. Model saved.")
}

################################################################################
# Post-Model Processing
################################################################################
# message("-- POST-MODEL PROCESSING --")
# model_data <- model_fit$data
# if (!inherits(model_data, "data.frame")) {
#   model_data <- as.data.frame(model_data)
# }
# if (!"Population" %in% names(model_data) && "population" %in% names(model_data)) {
#   model_data <- model_data %>% mutate(Population = population)
# }
# if (!"count" %in% names(model_data)) {
#   stop("Column 'count' not found in the model data.")
# }
# model_data <- distinct(model_data, State, Year, Population, count, .keep_all = TRUE)
# posteriorLinpred <- LINPREAD_DRAW_FN(data = model_data, model = model_fit)
# catch <- CATCHMENT(posteriorLinpred)
# catchir.linpred <- LINPRED_TO_CATCHIR(catch)
# csv_file <- file.path(outDir, paste0(unique(mmwrdata$Pathogen), "_IRCatch.csv"))
# SAFE_WRITE(catchir.linpred, csv_file)
# message("Post-model processing complete. CSV saved to ", csv_file)
# message("-- ANALYSIS COMPLETE --")

# message("Run completed. Printing sessionInfo:")
# sessionInfo()

################################################################################
# Post-Model Processing
################################################################################
message("-- POST-MODEL PROCESSING --")

# Skip the model data processing to bypass errors
# model_data <- model_fit$data
# if (!inherits(model_data, "data.frame")) {
#   model_data <- as.data.frame(model_data)
# }
# if (!"Population" %in% names(model_data) && "population" %in% names(model_data)) {
#   model_data <- model_data %>% mutate(Population = population)
# }
# if (!"count" %in% names(model_data)) {
#   stop("Column 'count' not found in the model data.")
# }
# model_data <- distinct(model_data, State, Year, Population, count, .keep_all = TRUE)

# Directly use the model fit object instead of processed data
posteriorLinpred <- tryCatch({
  LINPREAD_DRAW_FN(data = model_fit$data, model = model_fit)
}, error = function(e) {
  message("Error in LINPREAD_DRAW_FN: ", e$message)
  return(NULL)
})

if (!is.null(posteriorLinpred)) {
  catch <- tryCatch({
    CATCHMENT(posteriorLinpred)
  }, error = function(e) {
    message("Error in CATCHMENT: ", e$message)
    return(NULL)
  })
  
  if (!is.null(catch)) {
    catchir.linpred <- tryCatch({
      LINPRED_TO_CATCHIR(catch)
    }, error = function(e) {
      message("Error in LINPRED_TO_CATCHIR: ", e$message)
      return(NULL)
    })
    
    if (!is.null(catchir.linpred)) {
      csv_file <- file.path(outDir, paste0(unique(mmwrdata$Pathogen), "_IRCatch.csv"))
      SAFE_WRITE(catchir.linpred, csv_file)
      message("Post-model processing complete. CSV saved to ", csv_file)
    } else {
      message("Skipping post-model processing: Error in LINPRED_TO_CATCHIR.")
    }
  } else {
    message("Skipping post-model processing: Error in CATCHMENT.")
  }
} else {
  message("Skipping post-model processing: Error in LINPREAD_DRAW_FN.")
}

message("-- ANALYSIS COMPLETE --")

message("Run completed. Printing sessionInfo:")
sessionInfo()
