# functions.R

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(HDInterval)
library(brms)

# Confirmation that functions.R has been sourced
print("functions.R has been successfully sourced.")

# PATH_ANALYSIS function
PATH_ANALYSIS <- function(mmwrdata, census, outBase) {
  print("Starting PATH_ANALYSIS function.")
  
  # Validate presence of required columns in mmwrdata
  required_mmwr_cols <- c("Pathogen", "STEC_Class", "State", "Year", "Cste")
  missing_mmwr_cols <- setdiff(required_mmwr_cols, names(mmwrdata))
  if (length(missing_mmwr_cols) > 0) {
    stop(paste("Missing columns in mmwrdata:", paste(missing_mmwr_cols, collapse = ", ")))
  }
  
  # Validate presence of required columns in census
  required_census_cols <- c("Year", "State", "Population")
  missing_census_cols <- setdiff(required_census_cols, names(census))
  if (length(missing_census_cols) > 0) {
    stop(paste("Missing columns in census:", paste(missing_census_cols, collapse = ", ")))
  }
  
  # Print column names for debugging
  print("Column names in mmwrdata:")
  print(names(mmwrdata))
  print("Column names in census:")
  print(names(census))
  
  # Define the list of pathogens excluding Listeria
  pathogens <- c("CAMPYLOBACTER", "CYCLOSPORA",
                "SALMONELLA", "SHIGELLA", "STEC", "VIBRIO", "YERSINIA")
  
  # Create separate datasets for each group of pathogens
  ## All pathogens but Cyclospora and Listeria
  df1 <- mmwrdata %>%
    filter(Pathogen %in% pathogens & Pathogen != "CYCLOSPORA")
  
  ## Make a dataset for STEC O157
  df2 <- mmwrdata %>%
    filter(Pathogen == "STEC" & STEC_Class == "STEC O157") %>%
    mutate(Pathogen = "STEC O157")
  
  ## Make a dataset for STEC NONO157
  df3 <- mmwrdata %>%
    filter(STEC_Class %in% c("STEC NONO157", "STEC O AG UNDET")) %>%
    mutate(Pathogen = "STEC NONO157")
  
  # Make a Listeria dataset
  df4 <- mmwrdata %>%
    filter(Pathogen == "LISTERIA" & Cste == "YES")
  
  # Combine datasets using bind_rows without .id to preserve 'Pathogen' column
  selectDf <- bind_rows(df1, df2, df3, df4)
  
  # Clean up intermediate data frames
  rm(df1, df2, df3, df4)
  
  # Debugging: Check the structure of selectDf after binding
  print("Structure of selectDf after bind_rows:")
  print(str(selectDf))
  
  # Revise year
  selectDf <- selectDf %>%
    mutate(
      yearn = as.numeric(as.character(Year)),
      year = as.factor(Year)
    )
  
  # Split the dataset by pathogen and summarize counts using split()
  selectDf <- selectDf %>%
    split(.$Pathogen) %>%
    lapply(function(x) {
      x <- x %>%
        group_by(year, State) %>%
        summarize(count = n(), .groups = "drop") %>%
        # Ensure that State-year combinations with 0 cases are represented
        complete(year, State = unique(mmwrdata$State), fill = list(count = 0)) %>%
        distinct() # Keep only unique rows
      return(x)
    }) %>%
    # Convert from a list of datasets to a single dataset
    bind_rows(.id = "Pathogen")
  
  # Debugging: Check the structure of selectDf after summarizing
  print("Structure of selectDf after summarizing:")
  print(str(selectDf))
  
  # Drop year-State combinations from the dataset for years
  # before the given state entered the FoodNet catchment
  selectDf <- selectDf %>%
    mutate(year = as.numeric(as.character(year))) %>%
    filter(
      (State == "CA") | 
      (State == "CO" & year >= 2001) |
      (State == "CT") | 
      (State == "GA") |
      (State == "MD" & year >= 1998) |
      (State == "MN") | 
      (State == "NM" & year >= 2004) |
      (State == "NY" & year >= 1998) | 
      (State == "OR") |
      (State == "TN" & year >= 2000)
    )
  
  # Prepare census data for joining
  census_prepared <- census %>%
    select(Year, State, Population) %>%
    mutate(
      year = as.numeric(as.character(Year)),
      State = as.character(State)
    ) %>%
    group_by(year, State) %>%
    summarize(Population = sum(Population), .groups = "drop")
  
  # Ensure data types match
  selectDf <- selectDf %>%
    mutate(
      State = as.character(State),
      year = as.numeric(as.character(year))
    )
  
  census_prepared <- census_prepared %>%
    mutate(
      State = as.character(State),
      year = as.numeric(as.character(year))
    )
  
  # Debugging: Check the structure before merging
  print("Structure of selectDf before merging:")
  print(str(selectDf))
  print("Structure of census_prepared before merging:")
  print(str(census_prepared))
  
  # Merge the aggregated FoodNet data with census data using left_join
  mergedDf <- selectDf %>%
    left_join(census_prepared, by = c("year", "State")) %>%
    distinct()
  
  # Debugging: Check the structure of mergedDf
  print("Structure of mergedDf after merging:")
  print(str(mergedDf))
  
  # Ensure 'yearn' is present and numeric
  mergedDf <- mergedDf %>%
    mutate(yearn = year)  # 'year' is already numeric
  
  # Debugging: Confirm 'yearn' is correctly assigned
  print("Checking 'yearn' column:")
  print(head(mergedDf$yearn))
  
  # Save merged data using paste0 to avoid directory issues
  saveFile <- paste0(outBase, "mmwrdata.csv")
  write.csv(mergedDf, saveFile, row.names = FALSE)
  
  print(paste("Merged data saved to:", saveFile))
  
  return(mergedDf)
}

# PROPOSED_BM function (Actual Model Fitting)
PROPOSED_BM <- function(data) {
  print("Entering PROPOSED_BM function.")
  
  # Check if Population column exists to avoid errors in the model
  if (!"Population" %in% names(data)) {
    stop("Population column is missing in the data.")
  }
  
  # Check if 'yearn' is numeric
  if (!is.numeric(data$yearn)) {
    stop("'yearn' must be numeric.")
  }
  
  # Define the model formula
  formula <- bf(count ~ s(yearn, bs = "ps") + offset(log(Population)), family = poisson())
  
  print("Model formula defined.")
  
  # Fit the model
  model <- brm(
    formula,
    data = data,
    chains = 4,
    cores = 4,
    iter = 2000,    # Default number of iterations
    warmup = 1000,  # Default number of warmup iterations
    seed = 123       # For reproducibility
  )
  
  print("Model fitting completed.")
  
  return(model)
}

# LINPRED_DRAW_FN function
LINPRED_DRAW_FN <- function(data, model) {
  print("Entering LINPRED_DRAW_FN function.")
  
  linpred <- data %>%
    add_predicted_draws(model, re_formula = NA, allow_new_levels = TRUE)
  
  print("Predicted draws added.")
  
  return(linpred)
}
