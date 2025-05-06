# FUNCTIONS.R
################################################################################
# can we add text similar to calcIR.R? maybe for each function that explains how/where it fits into the workflow? I am having trouble linking these functions to trendy.R and calcIR.R
# functions.R
#
# Purpose:
#   This script ...
#
#   This includes:
#     - ...
#
# Usage:
#   ...
#
# Example:
#   ...
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
  library(HDInterval)  # for hdi() function
  library(gridExtra)  # for arranging multiple plots
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

# PATH_ANALYSIS function
PATH_ANALYSIS <- function(mmwrdata, census) {
  pathogens <- c("CAMPYLOBACTER", "CYCLOSPORA", "SALMONELLA", "SHIGELLA", "STEC", "VIBRIO", "YERSINIA", "LISTERIA", "STEC", "STEC NONO157","STEC O157")
  
  selectDf <- mmwrdata %>%
    filter(pathogen %in% pathogens) %>%
    group_by(year, state, pathogen) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(year, state, pathogen = unique(pathogen), fill = list(count = 0)) %>%
    left_join(census %>% filter(pathogentype == "Bacterial"), by = c("year", "state")) %>%
    mutate(year = as.numeric(as.character(year))) %>%
    
    # Drop year-state combinations from the dataset for years before the given state entered the FoodNet catchment
    ## To make the function more flexible for non-FoodNet datasets, is there a way we can have users upload a file with a column for year and a column for state
    ## that can be used in this step to drop states in years where they weren't part of a catchment? This is a want not a need. I'd be interested in learning how to do this too
    ## maybe we could do it on a Teams session together?
    subset((state=="CA") | (state=="CO" & year>=2001) | (state=="CT") | (state=="GA") | (state=="MD" & year>=1998) | (state=="MN") | (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | (state=="TN" & year>=2000))
  
  return(selectDf)
}

# CYCLOSPORA_ANALYSIS function
CYCLOSPORA_ANALYSIS <- function(mmwrdata, census) {
  cyclo <- mmwrdata %>%
    filter(pathogen == "CYCLOSPORA") %>%
    group_by(year, state) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(year, state, fill = list(count = 0)) %>%
    left_join(census %>% filter(pathogentype == "Parasitic"), by = c("year", "state"))%>%
    
    # Drop year-state combinations from the dataset for years before the given state entered the FoodNet catchment
    ## To make the function more flexible for non-FoodNet datasets, is there a way we can have users upload a file with a column for year and a column for state
    ## that can be used in this step to drop states in years where they weren't part of a catchment? This is a want not a need. I'd be interested in learning how to do this too
    ## maybe we could do it on a Teams session together?
    subset((state=="CA") | (state=="CO" & year>=2001) | (state=="CT") | (state=="GA") | (state=="MD" & year>=1998) | (state=="MN") | (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | (state=="TN" & year>=2000))
  
  return(cyclo)
}

# SALMONELLA_ANALYSIS function
SALMONELLA_ANALYSIS <- function(mmwrdata, census) {
  sal <- mmwrdata %>%
    filter(pathogen == "SALMONELLA") %>%
    group_by(year, state) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(year, state, fill = list(count = 0)) %>%
    left_join(census %>% filter(pathogentype == "Bacterial"), by = c("year", "state"))%>%
    
    # Drop year-state combinations from the dataset for years before the given state entered the FoodNet catchment
    ## To make the function more flexible for non-FoodNet datasets, is there a way we can have users upload a file with a column for year and a column for state
    ## that can be used in this step to drop states in years where they weren't part of a catchment? This is a want not a need. I'd be interested in learning how to do this too
    ## maybe we could do it on a Teams session together?
    subset((state=="CA") | (state=="CO" & year>=2001) | (state=="CT") | (state=="GA") | (state=="MD" & year>=1998) | (state=="MN") | (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | (state=="TN" & year>=2000))
  
  return(sal)
}

# PROPOSED_BM function - Updated to handle zero-count data
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
  
  # Ensure population is numeric
  data$population <- as.numeric(data$population)
  
  # Ensure count is integer
  data$count <- as.integer(as.numeric(data$count))
  
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
      count ~ s(year, by = state) + state + offset(log(population)),
      data = data,
      family = negbinomial(),
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

# LINPREAD_DRAW_FN function
## Draw untransformed (link-level) predictions for a new (or the original) data using add_linpred (which is an alternate spelling of add_fitted_draws) and transform them
## This generates a distribution of estimates for each site
LINPREAD_DRAW_FN <- function(data, model) {
  # Prepare data: convert to tibble, ungroup, add a row identifier, and force the Population column to be numeric.
  data <- as_tibble(data) %>%
    ungroup() %>%
    mutate(
      .row = row_number(),
      Population = if ("Population" %in% names(.)) {
        parse_number(as.character(Population))
      } else if ("population" %in% names(.)) {
        parse_number(as.character(population))
      } else {
        stop("No population column found")
      }
    )
  
  # Ensure that Population is numeric and no NA values were introduced.
  if (!is.numeric(data$Population) || any(is.na(data$Population))) {
    stop("Population column is not numeric after conversion")
  }
  
  # Get posterior predictive draws (using tidybayes's epred_draws).
  tryCatch({
    epred <- epred_draws(model, newdata = data) %>% ungroup()
    
    # Remove any Population column in the posterior draws to avoid conflict.
    epred <- epred %>% select(-one_of("Population"))
    
    # Join the Population values back by the unique row identifier.
    pop_df <- data %>% select(.row, Population) %>% ungroup()
    draws <- left_join(epred, pop_df, by = ".row") %>% ungroup()
    
    if (!is.numeric(draws$Population) || any(is.na(draws$Population))) {
      stop("Population column is not numeric in the joined data")
    }
    
    # Now compute predicted incidence (per 100,000)
    draws <- draws %>% mutate(pred_incidence = .epred / (Population / 100000))
    
    return(draws)
  }, error = function(e) {
    # If prediction fails, create synthetic draws
    stop("Error generating predictions")
  })
}

# Implementation of CATCHMENT function
## Convert draws from site-level to catchment-level estimates
## This uses the output from LINPREAD_DRAW_FN
CATCHMENT <- function(draws) {
  # Group by relevant variables and calculate summary statistics
  catchment_data <- draws %>%
    group_by(year, .draw) %>%
    summarise(
      count = sum(count),
      population = sum(population),
      .epred = sum(.epred),
      Population = sum(Population),
      .groups = "drop"
    )
  
  return(catchment_data)
}

# Implementation of LINPRED_TO_CATCHIR function
## Convert catchment-level draws to catchment-level estimates, including equal-tailed credibility interval
LINPRED_TO_CATCHIR <- function(catchment_data) {
  ir_data<-catchment_data %>% 
    group_by(year) %>% 
    summarise(
      raw_count=median(count),
      raw_check=sd(count),
      median=median(.epred),
      mean=mean(.epred),
      lower_equitailed=quantile(.epred, probs = 0.025, na.rm=TRUE),
      upper_equitailed=quantile(.epred, probs = 0.975, na.rm=TRUE),
      lower_hdi = (hdi(.epred, credMass = 0.95)[1]),
      upper_hdi = (hdi(.epred, credMass = 0.95)[2]),
      population=mean(population),
      population_check=sd(population))%>%
    mutate(
      median_ir= round(median/(population/100000),2),
      mean_ir= round(mean/(population/100000),2),
      lower_equitailed_ir=round(lower_equitailed/(population/100000),2),
      upper_equitailed_ir=round(upper_equitailed/(population/100000),2),
      lower_hdi_ir = round(lower_hdi/(population/100000),2),
      upper_hdi_ir = round(upper_hdi/(population/100000),2)) %>%
    # Arrange by Year and State for better readability
    arrange(year)
  return(ir_data)
}

# Implementation of LINPRED_TO_SITEIR function
## Convert catchment-level draws to catchment-level estimates, including equal-tailed credibility interval
LINPRED_TO_SITEIR <- function(site_data) {
  ir_data<-site_data %>% 
    group_by(year, state) %>% 
    summarise(
      raw_count=median(count),
      raw_check=sd(count),
      median=median(.epred),
      mean=mean(.epred),
      lower_equitailed=quantile(.epred, probs = 0.025, na.rm=TRUE),
      upper_equitailed=quantile(.epred, probs = 0.975, na.rm=TRUE),
      lower_hdi = (hdi(.epred, credMass = 0.95)[1]),
      upper_hdi = (hdi(.epred, credMass = 0.95)[2]),
      population=mean(population),
      population_check=sd(population))%>%
    mutate(
      median_ir= round(median/(population/100000),2),
      mean_ir= round(mean/(population/100000),2),
      lower_equitailed_ir=round(lower_equitailed/(population/100000),2),
      upper_equitailed_ir=round(upper_equitailed/(population/100000),2),
      lower_hdi_ir = round(lower_hdi/(population/100000),2),
      upper_hdi_ir = round(upper_hdi/(population/100000),2)
    ) %>%
    # Arrange by Year and State for better readability
    arrange(year, state)
  return(ir_data)
}

##### stopped work here.

# New function: Plot site-specific trends
PLOT_SITE_TRENDS <- function(catchir_data, pathogen, outDir) {
  # Create a plot for each state showing trends over time
  p <- ggplot(catchir_data, aes(x = year, y = median_ir)) +
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
    geom_vline(xintercept = 2004)+
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

# New function: Plot overall trend
PLOT_OVERALL_TREND <- function(catchir_data, pathogen, outDir) {
  # Calculate overall incidence by year (weighted by population)
  overall_data <- catchir_data %>%
    group_by(year) %>%
    summarise(
      median_incidence = mean(median_incidence),
      lower_hdi = mean(lower_hdi),
      upper_hdi = mean(upper_hdi),
      .groups = "drop"
    )
  
  # Create the plot
  p <- ggplot(overall_data, aes(x = year, y = median_ir)) +
    geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = lower_hdi_ir, ymax = upper_hdi_ir), alpha = 0.3) +
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

# New function: Create a combined visualization
PLOT_COMBINED <- function(site_plot, overall_plot, pathogen, outDir) {
  # Combine the plots
  combined_plot <- gridExtra::grid.arrange(overall_plot, site_plot,
                                           ncol = 1, heights = c(1, 2))
  
  # Save the combined plot
  plot_file <- file.path(outDir, paste0(pathogen, "_combined.png"))
  ggsave(plot_file, combined_plot, width = 12, height = 14, dpi = 300)
  
  return(combined_plot)
}

# Implementation of IR_COMP function for calculating relative risks
IR_COMP_CATCH <- function(catch, catchir_data, start_year, end_year, output_file = NULL) {
  
  # Filter data for the comparison period
  period_data <- catch %>%
    filter(year >= start_year & year <= end_year)%>% group_by(.draw)%>%
    summarise_at(.vars = c("count", "population", ".epred"), 
                 .funs = list(mean=mean))%>%
    mutate(baseline_ir=.epred_mean/(population_mean/100000))
  colnames(period_data)<-c(".draw", "baseline_count", "baseline_population", "baseline_.epred", "baseline_ir")
  
  # Check if we have data for the requested period
  if (nrow(period_data) == 0) {
    stop(paste("No data available for period", start_year, "to", end_year))
  }
 # Join the data for the baseline period with the data for all other years. This is because you do all calcualtions on the draws THEN average
 comb<-left_join(catch, period_data, by=c(".draw"))%>%
   group_by(year)%>%
 # calculate the IR for each draw and year
   mutate(
     raw_ir=count/(population/100000),
     est_ir=.epred/(population/100000))%>%
 # calculate relative risk and percent change by draw
   mutate(relative_risk=  est_ir/baseline_ir,
          percent_change= ((est_ir / baseline_ir) - 1) * 100)%>%
 # extract estimates from the draws
   summarise(
      count=mean(count),
      population=mean(population),
      raw_ir=mean(raw_ir),
      baseline_count=mean(baseline_count),
      baseline_population=mean(baseline_population),
      baseline_.epred_lower_hdi = (hdi(baseline_.epred, credMass = 0.95)[1]),
      baseline_.epred_upper_hdi = (hdi(baseline_.epred, credMass = 0.95)[2]),
      baseline_.epred_est=mean(baseline_.epred),
      baseline_ir_lower_hdi = (hdi(baseline_ir, credMass = 0.95)[1]),
      baseline_ir_upper_hdi = (hdi(baseline_ir, credMass = 0.95)[2]),
      baseline_ir=mean(baseline_ir),
      est_ir_lower_hdi = (hdi(est_ir, credMass = 0.95)[1]),
      est_ir_upper_hdi = (hdi(est_ir, credMass = 0.95)[2]),
      est_ir=mean(est_ir),
      relative_risk_lower_hdi = (hdi(relative_risk, credMass = 0.95)[1]),
      relative_risk_upper_hdi = (hdi(relative_risk, credMass = 0.95)[2]),
      relative_risk_est=mean(relative_risk),
      percent_change_lower_hdi = (hdi(percent_change, credMass = 0.95)[1]),
      percent_change_upper_hdi = (hdi(percent_change, credMass = 0.95)[2]),
      percent_change_est=mean(percent_change))%>% # should we do the mean or median?
    mutate(comparison_period = paste0(start_year, "-", end_year))
  # Calculate relative risks for each year in the dataset relative to the baseline period
  # latest_year <- max(catchir_data$year) $ if you only want the more recent year, you can modify the code to use "latest_year"
  # Get the most recent year's data
  # latest_data <- catchir_data %>% filter(year == latest_year)
  
  # Round numeric columns for readability
  result <- comb %>%
    mutate(across(where(is.numeric), ~round(., 2)))
  
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
