#!/usr/bin/env Rscript
################################################################################
# functions.R - Core Statistical Functions for FoodNet Trends Pipeline
#
# Purpose:
#   This script provides all statistical modeling and analysis functions used by
#   trendy.R in the FoodNet Trends pipeline. It contains the Bayesian hierarchical
#   model implementations, data processing functions, and visualization utilities.
#
# Workflow Integration:
#   1. Sourced by trendy.R at the beginning of execution
#   2. Uses cleaned data from preprocess.R (lowercase column names expected)
#   3. Provides pathogen-specific analysis functions called by trendy.R
#   4. Generates all statistical outputs (models, plots, CSV files)
#
# Key Function Groups:
#   - Utility Functions: CLEAN_LIST, SAFE_WRITE
#   - Pathogen Analysis: PATH_ANALYSIS, CYCLOSPORA_ANALYSIS, SALMONELLA_ANALYSIS
#   - Bayesian Modeling: PROPOSED_BM (main model fitting function)
#   - Post-processing: LINPREAD_DRAW_FN, CATCHMENT, LINPRED_TO_CATCHIR, LINPRED_TO_SITEIR
#   - Visualization: PLOT_SITE_TRENDS, PLOT_OVERALL_TREND
#   - Results Generation: IR_COMP_CATCH
#
# Usage:
#   This file is automatically sourced by trendy.R and should not be run directly.
#   All functions expect data with lowercase column names as produced by preprocess.R.
#
# Example workflow:
#   preprocess.R → clean_mmwr.csv → trendy.R → sources functions.R → analysis outputs
#
################################################################################

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
  library(readr)  # for parse_number()
  library(gridExtra)  # for arranging multiple plots
})

# Try to load HDInterval for hdi() function
hdi_available <- requireNamespace("HDInterval", quietly = TRUE)
if (hdi_available) {
  suppressPackageStartupMessages(library(HDInterval))
} else {
  # Define a fallback hdi function using quantiles
  hdi <- function(x, credMass = 0.95) {
    # Simple approximation using equal-tailed intervals
    alpha <- 1 - credMass
    c(quantile(x, probs = alpha/2, na.rm = TRUE),
      quantile(x, probs = 1 - alpha/2, na.rm = TRUE))
  }
  warning("HDInterval package not available. Using quantile-based approximation for HDI.")
}

# --- Catchment Configuration Functions ---
# FoodNet is an active surveillance network monitoring foodborne illnesses 
# in specific geographic areas. States joined at different times:
# - Original sites (1996): CT, GA, MN, OR, selected counties in CA
# - Later additions: MD (1998), NY (1998), TN (2000), CO (2001), NM (2004)
# This configuration ensures analyses only include years with active surveillance
read_catchment_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    # Return default FoodNet configuration
    return(data.frame(
      state = c("CA", "CO", "CT", "GA", "MD", "MN", "NM", "NY", "OR", "TN"),
      start_year = c(1996, 2001, 1996, 1996, 1998, 1996, 2004, 1998, 1996, 2000),
      end_year = rep(9999, 10),  # 9999 indicates ongoing participation
      pathogen_type = rep("both", 10),  # States monitor both bacterial and parasitic pathogens
      stringsAsFactors = FALSE
    ))
  }
  # Read and validate CSV
  config <- read.csv(config_path, stringsAsFactors = FALSE)
  validate_catchment_config(config)
  return(config)
}

validate_catchment_config <- function(config) {
  required_cols <- c("state", "start_year", "end_year")
  missing_cols <- setdiff(required_cols, names(config))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in catchment config: ", paste(missing_cols, collapse = ", "))
  }
  
  # Validate years are numeric
  if (!is.numeric(config$start_year) || !is.numeric(config$end_year)) {
    stop("start_year and end_year must be numeric")
  }
  
  # Validate start_year <= end_year
  invalid_years <- config$start_year > config$end_year
  if (any(invalid_years)) {
    stop("start_year must be less than or equal to end_year for all entries")
  }
  
  # Add pathogen_type column if missing
  if (!"pathogen_type" %in% names(config)) {
    config$pathogen_type <- "both"
  }
  
  return(config)
}

apply_catchment_filter <- function(data, catchment_config, pathogen_type = "both") {
  # Filter config for relevant pathogen type
  if (pathogen_type != "both" && "pathogen_type" %in% names(catchment_config)) {
    relevant_config <- catchment_config[
      catchment_config$pathogen_type %in% c("both", pathogen_type), 
    ]
  } else {
    relevant_config <- catchment_config
  }
  
  # Build dynamic filter using vectorized operations
  valid_rows <- rep(FALSE, nrow(data))
  
  for (i in seq_len(nrow(relevant_config))) {
    state_match <- data$state == relevant_config$state[i]
    year_match <- data$year >= relevant_config$start_year[i] & 
                  data$year <= relevant_config$end_year[i]
    valid_rows <- valid_rows | (state_match & year_match)
  }
  
  return(data[valid_rows, ])
}

# Helper: Clean up list strings and handle vector inputs
# Used by trendy.R to parse comma-separated pathogen and filter parameters
CLEAN_LIST <- function(input_string) {
  if (length(input_string) > 1) {
    # If input is a vector, collapse into a single string
    input_string <- paste(input_string, collapse = ",")
  }
  cleanedString <- gsub('[\\[\\]\"]', '', input_string)
  strsplit(cleanedString, ",")[[1]]
}

# Helper: Write data to a file safely
# Used throughout the pipeline to save outputs (CSV files, RDS objects)
# Handles directory creation and error logging
SAFE_WRITE <- function(data, file_path) {
  tryCatch({
    dir_path <- dirname(file_path)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    
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

################################################################################
# PATH_ANALYSIS - Prepare bacterial pathogen data for modeling
# 
# Called by: trendy.R for bacterial pathogens (CAMPYLOBACTER, SALMONELLA, etc.)
# Inputs: 
#   - mmwrdata: Cleaned MMWR data from preprocess.R (lowercase columns)
#   - census: Census data for population denominators
# Output: Data frame ready for PROPOSED_BM modeling
#
# Workflow position: trendy.R → PATH_ANALYSIS → PROPOSED_BM
################################################################################
PATH_ANALYSIS <- function(mmwrdata, census, catchment_config = NULL) {
  # Process all pathogens found in the data (not limited to a predefined list)
  # Get unique pathogens from the cleaned data (excluding those handled specially)
  all_pathogens <- unique(mmwrdata$pathogen)
  
  # Check if any pathogens exist
  if (length(all_pathogens) == 0) {
    stop("No pathogens found in the MMWR data")
  }
  
  # Define parasitic pathogens
  parasitic_pathogens <- c("CRYPTOSPORIDIUM", "CYCLOSPORA")
  
  selectDf <- mmwrdata %>%
    filter(pathogen %in% all_pathogens) %>%
    group_by(year, state, pathogen) %>%
    summarise(count = n(), .groups = "drop") %>%
    # complete() fills in zero counts for all year-state-pathogen combinations
    # This is intentional - FoodNet is an active surveillance system where all cases 
    # in participating sites are reported. Absence of reported cases in a 
    # participating state-year represents true zeros, not missing data
    complete(year, state, pathogen = unique(pathogen), fill = list(count = 0))
  
  # Join with appropriate census data based on pathogen type
  bacterial_data <- selectDf %>%
    filter(!pathogen %in% parasitic_pathogens) %>%
    left_join(census %>% filter(pathogentype == "Bacterial"), by = c("year", "state"))
  
  parasitic_data <- selectDf %>%
    filter(pathogen %in% parasitic_pathogens) %>%
    left_join(census %>% filter(pathogentype == "Parasitic"), by = c("year", "state"))
  
  # Combine the data
  selectDf <- bind_rows(bacterial_data, parasitic_data) %>%
    mutate(year = as.numeric(as.character(year)))
  
  # Check if the join produced any data with population
  if (all(is.na(selectDf$population))) {
    stop("No population data found after joining with census data. Check that census data contains 'Bacterial' pathogentype.")
  }
  
  # Remove rows with missing population data
  rows_before <- nrow(selectDf)
  selectDf <- selectDf %>% filter(!is.na(population))
  rows_after <- nrow(selectDf)
  
  if (rows_before > rows_after) {
    warning(paste("Removed", rows_before - rows_after, "rows with missing population data"))
  }
  
  if (nrow(selectDf) == 0) {
    stop("No data remaining after removing rows with missing population")
  }
    
  # Drop year-state combinations from the dataset for years before the given state entered the FoodNet catchment
  # Configurable via the catchment_config parameter to support different surveillance periods
  if (is.null(catchment_config)) {
    catchment_config <- read_catchment_config()
  }
  
  # Apply catchment filter  
  selectDf <- apply_catchment_filter(selectDf, catchment_config, "bacterial")
  
  return(selectDf)
}

################################################################################
# CYCLOSPORA_ANALYSIS - Prepare Cyclospora (parasitic) data for modeling
# 
# Called by: trendy.R specifically for CYCLOSPORA pathogen
# Inputs: 
#   - mmwrdata: Cleaned MMWR data from preprocess.R (lowercase columns)
#   - census: Census data with parasitic population denominators
# Output: Data frame ready for PROPOSED_BM modeling
#
# Note: Uses different census data (Parasitic) than bacterial pathogens
# Workflow position: trendy.R → CYCLOSPORA_ANALYSIS → PROPOSED_BM
################################################################################
CYCLOSPORA_ANALYSIS <- function(mmwrdata, census, catchment_config = NULL) {
  cyclo <- mmwrdata %>%
    filter(pathogen == "CYCLOSPORA") %>%
    group_by(year, state) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(year, state, fill = list(count = 0)) %>%
    left_join(census %>% filter(pathogentype == "Parasitic"), by = c("year", "state"))
  
  # Remove rows with missing population data
  cyclo <- cyclo %>% filter(!is.na(population))
    
  # Drop year-state combinations from the dataset for years before the given state entered the FoodNet catchment
  # Configurable via the catchment_config parameter to support different surveillance periods
  if (is.null(catchment_config)) {
    catchment_config <- read_catchment_config()
  }
  
  # Apply catchment filter for parasitic pathogens
  cyclo <- apply_catchment_filter(cyclo, catchment_config, "parasitic")
  
  return(cyclo)
}

################################################################################
# SALMONELLA_ANALYSIS - Prepare Salmonella-specific data for modeling
# 
# Called by: trendy.R specifically for SALMONELLA pathogen
# Inputs: 
#   - mmwrdata: Cleaned MMWR data from preprocess.R (lowercase columns)
#   - census: Census data with bacterial population denominators
# Output: Data frame ready for PROPOSED_BM modeling
#
# Note: Similar to PATH_ANALYSIS but Salmonella-specific
# Workflow position: trendy.R → SALMONELLA_ANALYSIS → PROPOSED_BM
################################################################################
SALMONELLA_ANALYSIS <- function(mmwrdata, census, catchment_config = NULL) {
  sal <- mmwrdata %>%
    filter(pathogen == "SALMONELLA") %>%
    group_by(year, state) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(year, state, fill = list(count = 0)) %>%
    left_join(census %>% filter(pathogentype == "Bacterial"), by = c("year", "state"))
  
  # Remove rows with missing population data
  sal <- sal %>% filter(!is.na(population))
    
  # Drop year-state combinations from the dataset for years before the given state entered the FoodNet catchment
  # Configurable via the catchment_config parameter to support different surveillance periods
  if (is.null(catchment_config)) {
    catchment_config <- read_catchment_config()
  }
  
  # Apply catchment filter for bacterial pathogens
  sal <- apply_catchment_filter(sal, catchment_config, "bacterial")
  
  return(sal)
}

################################################################################
# PROPOSED_BM - Main Bayesian hierarchical modeling function
# 
# Called by: trendy.R after data preparation by PATH/CYCLOSPORA/SALMONELLA_ANALYSIS
# Purpose: Fits a Bayesian hierarchical model with splines to estimate incidence rates
# 
# Inputs:
#   - data: Prepared data frame with count, year, state, population columns
#   - cores: Number of CPU cores for parallel processing
#   - chains: Number of MCMC chains (default: 2, paper uses 6)
#   - iterations: Iterations per chain (default: 500, paper uses 10,001)
#   - adapt_delta: HMC adaptation parameter (default: 0.95)
#   - max_treedepth: Maximum tree depth for HMC (default: 10)
#   - seed: Random seed for reproducibility
#
# Output: brms model object with posterior samples
#
# Workflow position: Data preparation functions → PROPOSED_BM → Post-processing
# 
# Memory scaling (approximate):
# - 2 chains: 24GB RAM
# - 4 chains: 48GB RAM
# - 6 chains: 56GB RAM (publication quality)
# - 8 chains: 72GB RAM
#
# Parameter guidelines:
# - adapt_delta: Increase (0.95-0.99) if divergent transitions occur
# - max_treedepth: Increase if hitting max treedepth warnings
# - iterations: Publication quality typically requires 5000-10000
################################################################################
PROPOSED_BM <- function(data, cores = 16, chains = 2, iterations = 500,
                        adapt_delta = 0.95, max_treedepth = 10, seed = 123) {
  # Ensure data is properly formatted
  data <- as.data.frame(data)
  
  # Verify required columns
  required_cols <- c("count", "year", "state", "population")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns in data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure population is numeric and handle NA values
  data$population <- as.numeric(data$population)
  if (any(is.na(data$population))) {
    stop("Population column contains NA values after conversion")
  }
  if (any(data$population <= 0)) {
    stop("Population column contains zero or negative values")
  }
  
  # Ensure count is integer
  data$count <- as.integer(as.numeric(data$count))
  if (any(is.na(data$count))) {
    stop("Count column contains NA values after conversion")
  }
  
  # Check if all counts are zero - this will cause model fitting issues
  if (all(data$count == 0) || sum(data$count) == 0) {
    stop("All counts zero. Cannot fit model")
  }
  
  # Ensure year is numeric (not factor) for the spline
  if (is.factor(data$year)) {
    data$year <- as.numeric(as.character(data$year))
  }
  
  # Convert state to factor if it isn't already
  if (!is.factor(data$state)) {
    data$state <- as.factor(data$state)
  }
  
  # Set a reasonable seed for reproducibility
  set.seed(seed)
  
  # Fit the model with more robust settings
  model <- tryCatch({
    brm(
      # Model formula explained:
      # count ~ s(year, by = state) + state + offset(log(population))
      # - count: observed case counts (response variable)
      # - s(year, by = state): state-specific smoothing splines over time
      # - state: state fixed effects (baseline differences)
      # - offset(log(population)): log population offset for rate modeling
      count ~ s(year, by = state) + state + offset(log(population)),
      data = data,
      family = negbinomial(),  # Handles overdispersion common in count data
      chains = chains,
      iter = iterations,
      cores = cores,
      seed = seed,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      backend = "rstan"  # Explicitly use rstan backend for stability
    )
  }, error = function(e) {
    stop("Model did not converge. May need to run a simpler version, use more iterations, or more robust adapt_delta/max_treedepth values.")
  })
  
  return(model)
}

################################################################################
# LINPREAD_DRAW_FN - Extract posterior predictions from fitted model
# 
# Called by: trendy.R after PROPOSED_BM completes successfully
# Purpose: Generate posterior draws of incidence rates for each state-year
# 
# Inputs:
#   - data: Original data used in model fitting
#   - model: Fitted brms model object from PROPOSED_BM
#
# Output: Data frame with posterior draws and calculated incidence rates
# Note: Uses add_linpred_draws to get untransformed predictions
#
# Workflow position: PROPOSED_BM → LINPREAD_DRAW_FN → CATCHMENT
################################################################################
LINPREAD_DRAW_FN <- function(data, model) {
  # Prepare data: convert to tibble, ungroup, add a row identifier
  data <- as_tibble(data) %>%
    ungroup() %>%
    mutate(.row = row_number())
  
  # Get posterior predictive draws (using tidybayes's epred_draws).
  tryCatch({
    draws <- epred_draws(model, newdata = data) %>% 
              ungroup()
    
    if (!is.numeric(draws$population) || any(is.na(draws$population))) {
      stop("Population column is not numeric in the joined data")
    }
    
    return(draws)
  }, error = function(e) {
    # If prediction fails, create synthetic draws
    stop("Error generating predictions")
  })
}

################################################################################
# CATCHMENT - Aggregate state-level draws to catchment-level draws. Calculate 
# catchment-level incidence
# 
# Called by: trendy.R after LINPREAD_DRAW_FN
# Purpose: Combine state-level posterior draws into catchment-level draws
# 
# Input: Posterior draws from LINPREAD_DRAW_FN (state-level)
# Output: Catchment-level aggregated draws
#
# Workflow position: LINPREAD_DRAW_FN → CATCHMENT → LINPRED_TO_CATCHIR
################################################################################
CATCHMENT <- function(draws) {
  # Group by relevant variables and calculate summary statistics
  catchment_data <- draws %>%
    group_by(year, .draw) %>%
    summarise(
      count = sum(count),
      population = sum(population),
      .epred = sum(.epred),
      .groups = "drop"
    )
  
  return(catchment_data)
}
################################################################################
# LINPRED_TO_CATCHIR - Generate catchment-level point estimates by summarizing 
# incidence across draws 
# 
# Called by: trendy.R after CATCHMENT
# Purpose: Convert posterior draws to point estimates and credible intervals
# 
# Input: Catchment-level draws from CATCHMENT function
# Output: Summary statistics (median, mean, CI) for catchment incidence rates
# Note: Uses median as primary estimate (more robust for skewed distributions)
#
# Workflow position: CATCHMENT → LINPRED_TO_CATCHIR → CSV output
################################################################################
LINPRED_TO_CATCHIR <- function(catchment_data) {
  ir_data<-catchment_data %>% 
    mutate(ir=.epred/(population/100000))%>%
    group_by(year)%>%
    summarise(
      # Population: This is a check, SD should be 0
      population=round(median(population),6),
      population_check=round(sd(population),6),
      # Raw or Reported Count: This is a check, SD should be 0
      raw_count=round(median(count),6),
      raw_check=round(sd(count),6),
      # Raw or Reported Incidence: This is a check, SD should be 0
      raw_ir=round(raw_count/(population/100000),6),
      # Estimated Count
      median=round(median(.epred),6),
      mean=round(mean(.epred),6),
      lower_equitailed=round(quantile(.epred, probs = 0.025, na.rm=TRUE),6),
      upper_equitailed=round(quantile(.epred, probs = 0.975, na.rm=TRUE),6),
      lower_hdi = round(hdi(.epred, credMass = 0.95)[1],6),
      upper_hdi = round(hdi(.epred, credMass = 0.95)[2],6),
      # Estimated Incidence
      median_ir= round(median(ir),6),
      mean_ir= round(mean(ir),6),
      lower_equitailed_ir=round(quantile(ir, probs = 0.025, na.rm=TRUE),6),
      upper_equitailed_ir=round(quantile(ir, probs = 0.975, na.rm=TRUE),6),
      lower_hdi_ir = round(hdi(ir, credMass = 0.95)[1],6),
      upper_hdi_ir = round(hdi(ir, credMass = 0.95)[2],6))%>%
    # Arrange by Year and State for better readability
    arrange(year)
  return(ir_data)
}

################################################################################
# LINPRED_TO_SITEIR - Calculate site-level point estimates by summarizing 
# incidence across draws 
#
# Called by: trendy.R after LINPREAD_DRAW_FN
# Purpose: Convert posterior draws to point estimates and credible intervals by state
# 
# Input: State-level draws from LINPREAD_DRAW_FN
# Output: Summary statistics (median, mean, CI) for state-specific incidence rates
# Note: Parallel to LINPRED_TO_CATCHIR but maintains state-level granularity
#
# Workflow position: LINPREAD_DRAW_FN → LINPRED_TO_SITEIR → PLOT_SITE_TRENDS
################################################################################
LINPRED_TO_SITEIR <- function(site_data) {
  ir_data<-site_data %>% 
    mutate(ir=.epred/(population/100000))%>%
    group_by(year, state)%>%
    summarise(
      # Population: This is a check, SD should be 0
      population=round(median(population),6),
      population_check=round(sd(population),6),
      # Raw or Reported Count: This is a check, SD should be 0
      raw_count=round(median(count),6),
      raw_check=round(sd(count),6),
      # Raw or Reported Incidence: This is a check, SD should be 0
      raw_ir=round(raw_count/(population/100000),6),
      # Estimated Count
      median=round(median(.epred),6),
      mean=round(mean(.epred),6),
      lower_equitailed=round(quantile(.epred, probs = 0.025, na.rm=TRUE),6),
      upper_equitailed=round(quantile(.epred, probs = 0.975, na.rm=TRUE),6),
      lower_hdi = round(hdi(.epred, credMass = 0.95)[1],6),
      upper_hdi = round(hdi(.epred, credMass = 0.95)[2],6),
      # Estimated Incidence
      median_ir= round(median(ir),6),
      mean_ir= round(mean(ir),6),
      lower_equitailed_ir=round(quantile(ir, probs = 0.025, na.rm=TRUE),6),
      upper_equitailed_ir=round(quantile(ir, probs = 0.975, na.rm=TRUE),6),
      lower_hdi_ir = round(hdi(ir, credMass = 0.95)[1],6),
      upper_hdi_ir = round(hdi(ir, credMass = 0.95)[2],6))%>%
    # Arrange by Year and State for better readability
    arrange(year)
  return(ir_data)
}

################################################################################
# PLOT_SITE_TRENDS - Generate state-specific trend plots
# 
# Called by: trendy.R after LINPRED_TO_SITEIR
# Purpose: Create faceted plots showing incidence trends by state
# 
# Inputs:
#   - site: State-level summary data from LINPRED_TO_SITEIR
#   - pathogen: Name of pathogen for plot title
#   - outDir: Output directory for saving plot
#
# Output: PNG file with state-specific trend plots
#
# Workflow position: LINPRED_TO_SITEIR → PLOT_SITE_TRENDS → PNG output
################################################################################
PLOT_SITE_TRENDS <- function(site, pathogen, outDir) {
  # Create a plot for each state showing trends over time
  p <- ggplot(site, aes(x = year, y = median_ir)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower_hdi_ir, ymax = upper_hdi_ir), alpha = 0.3) +
    facet_wrap(~ state, scales = "free_y") +
    labs(
      title = paste("Site-Specific Trends for", pathogen),
      subtitle = "Median incidence with 95% HDI intervals",
      y = "Incidence per 100,000 population",
      x = "Year"
    ) +
    theme_minimal() +
    geom_vline(aes(xintercept = 2004), linetype="dashed", color="red")+
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold")
    )
  
  # Save the plot
  plot_file <- file.path(outDir, paste0(pathogen, "_site_trends.png"))
  ggsave(plot_file, p, width = 12, height = 8, dpi = 300)
  
  return(p)
}

################################################################################
# PLOT_OVERALL_TREND - Generate catchment-wide trend plot
# 
# Called by: trendy.R after LINPRED_TO_CATCHIR
# Purpose: Create single plot showing overall catchment incidence trends
# 
# Inputs:
#   - catchir_data: Catchment-level summary data from LINPRED_TO_CATCHIR
#   - pathogen: Name of pathogen for plot title
#   - outDir: Output directory for saving plot
#
# Output: PNG file with catchment-wide trend plot
#
# Workflow position: LINPRED_TO_CATCHIR → PLOT_OVERALL_TREND → PNG output
################################################################################
PLOT_OVERALL_TREND <- function(catchir_data, pathogen, outDir) {
# Create the plot
  p <- ggplot(catchir_data, aes(x = year, y = median_ir)) +
    geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = lower_hdi_ir, ymax = upper_hdi_ir), alpha = 0.3) +
    geom_vline(aes(xintercept = 2004), linetype="dashed", color="red")+
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
  plot_file <- file.path(outDir, paste0(pathogen, "_overall_trend.png"))
  ggsave(plot_file, p, width = 10, height = 6, dpi = 300)
  
  return(p)
}

################################################################################
# IR_COMP_CATCH - Calculate incidence rate ratios between time periods
# 
# Called by: trendy.R to generate comparative statistics
# Purpose: Compare incidence rates between specified years (e.g., 2023 vs 1996-1998)
# 
# Inputs:
#   - catch: Raw posterior draws from CATCHMENT
#   - start_year, end_year: Years to compare
#   - output_file: Optional CSV output path
#
# Output: Data frame with IRR estimates and confidence intervals
# Note: Uses median() for robust estimation with skewed posteriors
#
# Workflow position: CATCHMENT → IR_COMP_CATCH → CSV output
################################################################################
IR_COMP_CATCH <- function(catch, start_year, end_year, output_file = NULL) {
  
  # Filter data for the comparison period
  period_data <- catch %>%
    filter(year >= start_year & year <= end_year)%>% group_by(.draw)%>%
    mutate(ir=.epred/(population/100000))%>%
    summarise(baseline_ir=median(ir),
              baseline_count=median(count))
  colnames(period_data)<-c(".draw", "baseline_ir", "baseline_count")
  
  # Check if we have data for the requested period
  if (nrow(period_data) == 0) {
    stop(paste("No data available for period", start_year, "to", end_year))
  }
 # Join the data for the baseline period with the data for all other years. This is because you do all calcualtions on the draws THEN average
 comb<-left_join(catch, period_data, by=c(".draw"))%>% 
   mutate(ir=.epred/(population/100000),
          est_ir=.epred/(population/100000),
          relative_risk=  est_ir/baseline_ir,
          percent_change= ((est_ir-baseline_ir)/baseline_ir)*100)%>%
   group_by(year)%>%
 # extract estimates from the draws
   summarise(
     # Population: This is a check, SD should be 0
     population=round(median(population),6),
     population_check=round(sd(population),6),
     # Raw or Reported Count: This is a check, SD should be 0
     raw_count=round(median(count),6),
     raw_check=round(sd(count),6),
     # Raw or Reported Incidence: This is a check, SD should be 0
     raw_ir=round(raw_count/(population/100000),6),
     # Estimated Count
     median=round(median(.epred),6),
     mean=round(mean(.epred),6),
     lower_equitailed=round(quantile(.epred, probs = 0.025, na.rm=TRUE),6),
     upper_equitailed=round(quantile(.epred, probs = 0.975, na.rm=TRUE),6),
     lower_hdi = round(hdi(.epred, credMass = 0.95)[1],6),
     upper_hdi = round(hdi(.epred, credMass = 0.95)[2],6),
     # Estimated Incidence
     median_ir= round(median(ir),6),
     mean_ir= round(mean(ir),6),
     lower_equitailed_ir=round(quantile(ir, probs = 0.025, na.rm=TRUE),6),
     upper_equitailed_ir=round(quantile(ir, probs = 0.975, na.rm=TRUE),6),
     lower_hdi_ir = round(hdi(ir, credMass = 0.95)[1],6),
     upper_hdi_ir = round(hdi(ir, credMass = 0.95)[2],6),
     # Relative Risk and Percent Change
     relative_risk_lower_hdi = round(hdi(relative_risk, credMass = 0.95)[1],6),
     relative_risk_upper_hdi = round(hdi(relative_risk, credMass = 0.95)[2],6),
     relative_risk_est=round(median(relative_risk),6),
     percent_change_lower_hdi = round(hdi(percent_change, credMass = 0.95)[1],6),
     percent_change_upper_hdi = round(hdi(percent_change, credMass = 0.95)[2],6),
     percent_change_est=round(median(percent_change),6))%>% # Using median as it's more robust for potentially skewed posterior distributions
     mutate(comparison_period = paste0(start_year, "-", end_year))
  # Calculate relative risks for each year in the dataset relative to the baseline period
  # Alternative approach: To calculate IRR for only the most recent year vs baseline:
  # latest_year <- max(catchir_data$year)
  # Then filter results to show only year == latest_year
  # Get the most recent year's data
  # latest_data <- catchir_data %>% filter(year == latest_year)
  
  # Round numeric columns for readability
  result <- comb %>%
    mutate(across(where(is.numeric), ~round(., 6)))
  
  # Write to file if specified
  if (!is.null(output_file)) {
    # Create directory if it doesn't exist
    dir_path <- dirname(output_file)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    
    write.csv(result, output_file, row.names = FALSE)
  }
  
  return(result)
}


# New function: Plot percent change trend
PLOT_PCTCHange_TREND <- function(hp30) {
  # Create the plot
  p <- ggplot(hp30, aes(x = year, y = relative_risk_est)) +
    geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = relative_risk_lower_hdi, ymax = relative_risk_upper_hdi), alpha = 0.3) +
    geom_vline(aes(xintercept = 2004), type="dashed", color="red")+
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
  plot_file <- file.path(outDir, paste0(pathogen, "_overall_trend.png"))
  ggsave(plot_file, p, width = 10, height = 6, dpi = 300)
  
  return(p)
}

# This allows you to pull in our outputs for the same time period and combine into a single file. Helpful for making multi-pathogen tables and graphs
combine_files<-function(file_path, pattern){
  dir_path <- dirname(file_path)
  setwd(dir_path) 
  df = list.files(all.files = T,  pattern = pattern, full.names = F, recursive = TRUE) 
  df %>%
    set_names(.) %>%
    map_df(~mutate_all(read.csv(.x), as.character), .id = 'grp') %>%
    mutate(grp = str_remove(basename(grp), ".xlsx")) %>%
    separate(grp, c('pathogen', 'drop'), sep = '_', extra = 'merge')%>%
    select(-c(drop)) -> datas
  return(datas)
}