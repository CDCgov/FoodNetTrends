# FUNCTIONS.R

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(gtools)
  library(brms)
  library(ggplot2)
  library(tidybayes)
  library(haven)
})

# Helper: Clean up list strings and handle vector inputs
CLEAN_LIST <- function(input_string) {
  if (length(input_string) > 1) {
    # If input is a vector, collapse into a single string
    input_string <- paste(input_string, collapse = ",")
  }
  cleanedString <- gsub('[\\[\\]\"]', '', input_string)
  strsplit(cleanedString, ",")[[1]]
}

# Helper: Write data to a file safely
SAFE_WRITE <- function(data, file_path) {
  tryCatch({
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

# PATH_ANALYSIS function
PATH_ANALYSIS <- function(mmwrdata, census) {
  pathogens <- c("CAMPYLOBACTER", "CYCLOSPORA", "SALMONELLA", "SHIGELLA", "STEC", "VIBRIO", "YERSINIA")

  selectDf <- mmwrdata %>%
    filter(Pathogen %in% pathogens) %>%
    group_by(Year, State, Pathogen) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(Year, State, Pathogen = unique(Pathogen), fill = list(count = 0)) %>%
    left_join(census %>% filter(pathogentype == "Bacterial"), by = c("Year", "State")) %>%
    mutate(Year = as.numeric(as.character(Year))) %>%
    filter(State %in% c("CA", "CO", "CT", "GA", "MD", "MN", "NM", "NY", "OR", "TN"))

  return(selectDf)
}

# CYCLOSPORA_ANALYSIS function
CYCLOSPORA_ANALYSIS <- function(mmwrdata, census) {
  cyclo <- mmwrdata %>%
    filter(Pathogen == "CYCLOSPORA") %>%
    group_by(Year, State) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(Year, State, fill = list(count = 0)) %>%
    left_join(census %>% filter(pathogentype == "Parasitic"), by = c("Year", "State"))

  return(cyclo)
}

# SALMONELLA_ANALYSIS function
SALMONELLA_ANALYSIS <- function(mmwrdata, census) {
  sal <- mmwrdata %>%
    filter(Pathogen == "SALMONELLA") %>%
    group_by(Year, State) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(Year, State, fill = list(count = 0)) %>%
    left_join(census %>% filter(pathogentype == "Bacterial"), by = c("Year", "State"))

  return(sal)
}

# PROPOSED_BM function
PROPOSED_BM <- function(data) {
  model <- brm(
    count ~ s(Year, by = State) + State + offset(log(Population)),
    data = data,
    family = negbinomial(),
    chains = 4,
    iter = 2000,
    cores = 4,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  return(model)
}

LINPREAD_DRAW_FN <- function(data, model) {
  # Add the posterior predictive draws (epred) to the data
  data <- add_epred_draws(model, data, re_formula = NULL)
  
  # Ensure data is a tibble or data.frame
  if (!inherits(data, "data.frame")) {
    data <- as.data.frame(data)
  }

  # Remove any grouping structure (if it's a tibble or grouped data)
  data <- ungroup(data)
  
  # Create the predicted incidence variable
  data <- data %>%
    mutate(pred_incidence = .value / (Population / 100000))
  
  return(data)
}


# Placeholder stubs for CATCHMENT and LINPRED_TO_CATCHIR
# These must be defined or implemented by you.
CATCHMENT <- function(data) {
  # Implement your logic here
  data
}

LINPRED_TO_CATCHIR <- function(data) {
  # Implement your logic here
  data
}

