#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("haven"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("tidybayes"))
suppressPackageStartupMessages(library("HDInterval"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))

# Set warning options
options(warn = -1)

sink("output_log.txt")  # Redirect output to a log file

# Dynamically determine number of cores to use based on available resources
available_cores <- parallel::detectCores()
model_cores <- min(available_cores, as.integer(Sys.getenv("NCPUS", available_cores)))

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
parser$add_argument("--use_splines", type = "logical", default = TRUE,
                    help = "Use spline-based model instead of linear model")
opts <- parser$parse_args()

# Make options globally available
assign("opts", opts, envir = .GlobalEnv)

# Load helper functions
source("/scicomp/home-pure/smn9/FoodNetTrends/bin/functions.R")

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
  outDir        <- opts$outDir
  pathogen_arg  <- opts$pathogen
  preprocessed  <- opts$preprocessed
  cleanFile     <- opts$cleanFile
  use_splines   <- opts$use_splines
} else {
  # Debug mode defaults
  mmwrFile      <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/mmwr9623_Jan2024.sas7bdat"
  censusFile_B  <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623.sas7bdat"
  censusFile_P  <- "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623_para.sas7bdat"
  projID        <- "20240705"
  travel        <- CLEAN_LIST(c("NO", "UNKNOWN", "YES"))
  cidt          <- CLEAN_LIST(c("CIDT+", "CX+", "PARASITIC"))
  outDir        <- "./"
  pathogen_arg  <- opts$pathogen
  preprocessed  <- opts$preprocessed
  cleanFile     <- opts$cleanFile
  use_splines   <- TRUE
}

# Create output directory if needed
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

# Define model based on use_splines flag
if (use_splines) {
  message("Using spline-based model")
  CURRENT_MODEL <- function(data) {
    model <- brm(
      count ~ s(Year, by = State) + State + offset(log(Population)),
      data = data,
      family = negbinomial(),
      chains = 4,
      iter = 2000,
      cores = model_cores,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    )
    return(model)
  }
} else {
  message("Using linear model")
  CURRENT_MODEL <- function(data) {
    model <- brm(
      count ~ Year + State + offset(log(Population)),
      data = data,
      family = negbinomial(),
      chains = 4,
      iter = 2000,
      cores = model_cores,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    )
    return(model)
  }
}

if (nrow(mmwrdata) == 0) {
  stop("No data remaining after filtering. Check the input or the --pathogen parameter.")
}

message("Fitting model for pathogen: ", unique(mmwrdata$Pathogen))
model_fit <- try(CURRENT_MODEL(mmwrdata), silent = TRUE)
if (inherits(model_fit, "try-error")) {
  stop(paste("Error fitting model for", unique(mmwrdata$Pathogen), ":", model_fit[1])) # Extract error message
} else {
  saveRDS(model_fit, file = file.path(outDir, paste0(unique(mmwrdata$Pathogen), "_brm.Rds")))
  message("Model fitting complete. Model saved.")
}

################################################################################
# Post-Model Processing
################################################################################
message("-- POST-MODEL PROCESSING --")

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

      # Generate visualizations if ggplot2 and gridExtra are available
      if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("gridExtra", quietly = TRUE)) {
        tryCatch({
          # Create site-specific trends plot
          site_plot <- ggplot(catchir.linpred, aes(x = Year, y = median_incidence)) +
            geom_line(linewidth = 1) +
            geom_ribbon(aes(ymin = lower_hdi, ymax = upper_hdi), alpha = 0.3) +
            facet_wrap(~ State, scales = "free_y") +
            labs(
              title = paste("Site-Specific Trends for", unique(mmwrdata$Pathogen)),
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

          # Save the site-specific plot
          site_plot_file <- file.path(outDir, paste0(unique(mmwrdata$Pathogen), "_site_trends.png"))
          ggsave(site_plot_file, site_plot, width = 12, height = 8, dpi = 300)
          message("Site-specific trends plot saved to ", site_plot_file)

          # Calculate overall incidence by year
          overall_data <- catchir.linpred %>%
            group_by(Year) %>%
            summarise(
              median_incidence = mean(median_incidence),
              lower_hdi = mean(lower_hdi),
              upper_hdi = mean(upper_hdi),
              .groups = "drop"
            )

          # Create overall trend plot
          overall_plot <- ggplot(overall_data, aes(x = Year, y = median_incidence)) +
            geom_line(linewidth = 1.5) +
            geom_ribbon(aes(ymin = lower_hdi, ymax = upper_hdi), alpha = 0.3) +
            labs(
              title = paste("Overall Trend for", unique(mmwrdata$Pathogen)),
              subtitle = "Median incidence with 95% HDI intervals",
              y = "Incidence per 100,000 population",
              x = "Year"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5)
            )

          # Save the overall plot
          overall_plot_file <- file.path(outDir, paste0(unique(mmwrdata$Pathogen), "_overall_trend.png"))
          ggsave(overall_plot_file, overall_plot, width = 10, height = 6, dpi = 300)
          message("Overall trend plot saved to ", overall_plot_file)

          # Combine the plots
          combined_plot <- gridExtra::grid.arrange(overall_plot, site_plot,
                                                  ncol = 1, heights = c(1, 2))

          # Save the combined plot
          combined_plot_file <- file.path(outDir, paste0(unique(mmwrdata$Pathogen), "_combined.png"))
          ggsave(combined_plot_file, combined_plot, width = 12, height = 14, dpi = 300)
          message("Combined visualization saved to ", combined_plot_file)
        }, error = function(e) {
          message("Error generating visualizations: ", e$message)
        })
      } else {
        message("Skipping visualizations: ggplot2 or gridExtra not available.")
      }
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
print(sessionInfo())

# Close the log file
sink()

