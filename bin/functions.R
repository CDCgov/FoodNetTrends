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

# LINPREAD_DRAW_FN function
LINPREAD_DRAW_FN <- function(data, model) {
  # Prepare newdata: convert to tibble, ungroup, add a row identifier,
  # and force the Population column to be numeric.
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
}

# Implementation of CATCHMENT function
CATCHMENT <- function(draws) {
  # Group by relevant variables and calculate summary statistics
  catchment_data <- draws %>%
    group_by(Year, State, .draw) %>%
    summarise(
      pred_incidence = mean(pred_incidence),
      .groups = "drop"
    ) %>%
    # Calculate HDI intervals for each Year/State combination
    group_by(Year, State) %>%
    summarise(
      mean_incidence = mean(pred_incidence),
      median_incidence = median(pred_incidence),
      lower_hdi = hdi(pred_incidence, credMass = 0.95)[1],
      upper_hdi = hdi(pred_incidence, credMass = 0.95)[2],
      .groups = "drop"
    )
  
  return(catchment_data)
}

# Implementation of LINPRED_TO_CATCHIR function
LINPRED_TO_CATCHIR <- function(catchment_data) {
  # Format the data for output
  ir_data <- catchment_data %>%
    mutate(
      Year = as.integer(Year),
      # Round numeric values to 2 decimal places
      mean_incidence = round(mean_incidence, 2),
      median_incidence = round(median_incidence, 2),
      lower_hdi = round(lower_hdi, 2),
      upper_hdi = round(upper_hdi, 2)
    ) %>%
    # Arrange by Year and State for better readability
    arrange(Year, State)
  
  return(ir_data)
}

# New function: Plot site-specific trends
PLOT_SITE_TRENDS <- function(catchir_data, pathogen, outDir) {
  # Create a plot for each state showing trends over time
  p <- ggplot(catchir_data, aes(x = Year, y = median_incidence)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower_hdi, ymax = upper_hdi), alpha = 0.3) +
    facet_wrap(~ State, scales = "free_y") +
    labs(
      title = paste("Site-Specific Trends for", pathogen),
      subtitle = "Median incidence with 95% HDI intervals",
      y = "Incidence per 100,000 population",
      x = "Year"
    ) +
    theme_minimal() +
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
    group_by(Year) %>%
    summarise(
      median_incidence = mean(median_incidence),
      lower_hdi = mean(lower_hdi),
      upper_hdi = mean(upper_hdi),
      .groups = "drop"
    )
  
  # Create the plot
  p <- ggplot(overall_data, aes(x = Year, y = median_incidence)) +
    geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = lower_hdi, ymax = upper_hdi), alpha = 0.3) +
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

