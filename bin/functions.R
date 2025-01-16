##################################################################
# LIBRARY
##################################################################
require(tidyverse)
require(tidybayes)
LOAD_PACKAGES <- function(pkgs) {
  package_vctr <- pkgs %>% stringr::str_split('\n') %>% unlist()
  invisible(lapply(package_vctr, library, character.only = TRUE))
}

##################################################################
# Helper function for safe file writing
##################################################################
SAFE_WRITE <- function(data, file_path) {
  tryCatch({
    if (endsWith(file_path, ".csv")) {
      if (file.exists(file_path)) {
        # Append without column names
        write.table(data, file = file_path, append = F, quote = TRUE, sep = ",",
                    col.names = FALSE, row.names = FALSE)
      } else {
        # Write with column names
        write.table(data, file = file_path, append = FALSE, quote = TRUE, sep = ",",
                    col.names = TRUE, row.names = FALSE)
      }
    } else if (endsWith(file_path, ".Rds")) {
      saveRDS(data, file = file_path)
    }
    message("Successfully wrote file: ", file_path)
  }, error = function(e) {
    message("Failed to write file: ", file_path)
    message("Error: ", e$message)
  })
}

# What does this function do?
CLEANUP_TEMP <- function(temp_dir) {
  tryCatch({
    unlink(temp_dir, recursive = TRUE, force = TRUE)
    message("Temporary directory cleaned: ", temp_dir)
  }, error = function(e) {
    message("Failed to clean temporary directory: ", temp_dir)
    message("Error: ", e$message)
  })
}

##################################################################
# REFORMATTING STRINGS
##################################################################
CLEAN_LIST <- function(input_string) {
  cleanedString <- gsub('[\\[\\]\"]', '', input_string)
  cleanedList <- strsplit(cleanedString, ",")[[1]]
  return(cleanedList)
}

##################################################################
# Creating a dataset for bacterial pathogen incidence by year and state
##################################################################
PATH_ANALYSIS<-function(mmwrdata,census){
  # Create a list of pathogen names excluding Listeria
  pathogens<-c("CAMPYLOBACTER", "CYCLOSPORA", 
               "SALMONELLA", "SHIGELLA", "STEC", "VIBRIO", "YERSINIA")
  
  # Create dataset that includes CAMPYLOBACTER, SALMONELLA, SHIGELLA, STEC, VIBRIO, and YERSINIA. We will create a separate dataset for Listeria because there are additional exclusions needed for ensuring Listeria cases meet the CSTE definition
  # We will need to create separate datasets for Cyclospora because a different census dataset is needed (using function below)
  ## Datasets for all pathogens but Cyclospora and Listeria
  df1=mmwrdata %>% filter(pathogen %in% pathogens
                          & pathogen!="CYCLOSPORA")
  ## Create datasets for STEC O157 and non-STEC O157
  df2=mmwrdata%>% filter(pathogen == "STEC" & stec_class=="STEC O157")%>%
    mutate(pathogen="STEC O157")
  df3=mmwrdata%>% filter(stec_class=="STEC NONO157" | stec_class== "STEC O AG UNDET")%>%
    mutate(pathogen="STEC NONO157")
  
  # make a Listeria dataset 
  df4=mmwrdata%>% filter(pathogen == "LISTERIA" & cste=="YES") # Basically for a Listeria infection to be included it needs to meet the CSTE definition
  
  # join together the four pathogen-specific datasets, then remove these four separate datasets
  selectDf=gtools::smartbind(df1,df2,df3,df4)
  remove(df1,df2,df3,df4)
  
  # revise year
  selectDf = selectDf %>%
    mutate(yearn=as.numeric(as.character(year)),
           year=as.factor(year))
  
  # split the dataset by pathogen
  # Create a function and map (aka loop) over each 
  # pathogen-specific dataset to summarize disease counts 
  # by year and state
  selectDf=selectDf %>%
    # split the dataset by pathogen
    split(.$pathogen)%>% 
    map(function(x) {
      x<-x%>%
        group_by(year, state) %>%
        dplyr::summarize(count = n())%>%
        ungroup()%>%
        # Make sure that state-year combinations with 0 cases 
        # for the given pathogen are represented in the dataset
        complete(year, state = unique(mmwrdata$state),
                 fill = list(count = 0))%>%
        as.data.frame()%>%distinct() # keep onlyunique rows
    }) %>%
    # convert from a list of datasets to a single dataset 
    bind_rows(.id = "pathogen")
  
  # Drop year-state combinations from the dataset for years 
  # before the given state entered the FoodNet catchment
  selectDf = selectDf %>%
    mutate(year=as.numeric(as.character(year))) %>%
    subset((state=="CA") | (state=="CO" & year>=2001) | 
             (state=="CT") | (state=="GA") | 
             (state=="MD" & year>=1998) | 
             (state=="MN") | (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | 
             (state=="TN" & year>=2000))
  
  # Join the aggregated FoodNet data to the census data
  mergedDf = selectDf%>% merge(census %>%
                                 filter(pathogentype=="Bacterial") %>%
                                 dplyr::select(year, state, 
                                               population) %>%
                                 mutate(year=as.numeric(
                                   as.character(year))), 
                               by=c("year", "state"), 
                               all.x=TRUE) %>%
    distinct() %>%
    as.data.frame()
  
  saveFile=paste0(outBase,"mmwrdata.csv")
  write.csv(mergedDf,saveFile)
  
  return(mergedDf)
}

##################################################################
# Creating a dataset for parasitic (i.e., Cyclospora) pathogen incidence by year and state
##################################################################
CYCLOSPORA_ANALYSIS<-function(mmwrdata,census){
  para<-mmwrdata %>% 
    filter(pathogen=="CYCLOSPORA")%>%
    mutate(yearn=as.numeric(as.character(year))) %>%
    split(.$pathogen)%>%
    map(function(x) {
      x<-x%>%
        group_by(year, state) %>%
        dplyr::summarize(count = n())%>%
        ungroup()%>%
        complete(year, state = unique(mmwrdata$state), 
                 fill = list(count = 0))%>% 
        as.data.frame()%>%distinct()}) %>%
    bind_rows(.id = "pathogen") %>%
    mutate(year=as.numeric(as.character(year))) %>%
    subset((state=="CA") | (state=="CO" & year>=2001) | 
             (state=="CT") | (state=="GA") | 
             (state=="MD" & year>=1998) | (state=="MN") | 
             (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | 
             (state=="TN" & year>=2000)) %>%
    left_join(census%>%filter(pathogentype=="Parasitic") %>%
                mutate(year=as.numeric(as.character(year))), 
              by=c("year", "state")) %>%
    as.data.frame()
  return(para)
}

##################################################################
# Creating a dataset for with incidence for each of the most frequently reported Salmonella serotypes
##################################################################
SALMONELLA_ANALYSIS<-function(pathDf,mmwrdata,census){
  pop = pathDf %>%
    filter(pathogen == "CAMPYLOBACTER")  %>%
    ungroup
  
  # adjust year to factor
  pop$year<-as.factor(pop$year)
  mmwrdata$year<-as.factor(mmwrdata$year)
  
  # filter most common
  sal.rank = mmwrdata %>%
    filter(pathogen=="SALMONELLA")%>%
    filter(serotypesummary!="NOT SEROTYPED" & serotypesummary!="") %>%
    filter(year==max(as.numeric(as.character(mmwrdata$year)), 
                     na.rm=TRUE))%>%
    group_by(serotypesummary)%>%
    dplyr::summarize(count = n())%>%
    as.data.frame()
  sal.rank<-sal.rank[order(sal.rank$count, decreasing = T),] %>%
    filter(serotypesummary!="NOT SEROTYPED" & serotypesummary!="")
  sal.rank<-c(sal.rank$serotypesummary[1:10]) # change the last number if you want estimates for the 20 most frequently reported serotypes instead of the 10 most frequently reported serotypes
  
  
  sal = mmwrdata%>% 
    filter(pathogen=="SALMONELLA" & serotypesummary %in% sal.rank) %>%
    mutate(pathogen=serotypesummary, 
           yearn=as.numeric(as.character(year)))%>%
    mutate(year=droplevels(year))%>%
    split(.$pathogen)%>%
    map(function(x) {
      x<-x%>%
        group_by(year, state) %>%
        dplyr::summarize(count = n())%>%
        ungroup()%>%
        complete(year, state = unique(mmwrdata$state), 
                 fill = list(count = 0))%>% 
        as.data.frame()%>%distinct()
    }) %>%
    bind_rows(.id = "pathogen") %>%
    mutate(year=as.numeric(as.character(year)))%>%
    subset((state=="CA") | (state=="CO" & year>=2001) | 
             (state=="CT") | (state=="GA") | 
             (state=="MD" & year>=1998) | (state=="MN") | 
             (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | 
             (state=="TN" & year>=2000))%>%
    left_join(pop%>%dplyr::select(year, state, population) %>%
                mutate(year=as.numeric(as.character(year))), 
              by=c("year", "state")) %>%
    as.data.frame()
  
  return(sal)
}


##################################################################
# BRM Modeling
##################################################################
PROPOSED_BM <- function(bact) {
  # Define output file for the model
  output_file <- paste0(outBase, "brm_model.Rds")
  
  # Log model configuration
  cat("Starting BRM model with the following configuration:\n")
  cat(paste0("  - Chains: ", 4, "\n"))
  cat(paste0("  - Iterations: ", 5001, "\n")) # 101 is too few iterations for convergence. 5001 was needed for less common etiologic agents.
  cat(paste0("  - Cores: ", 9, "\n"))
  cat(paste0("  - Seed: ", 47, "\n"))
  cat(paste0("  - Adapt Delta: ", 0.99, "\n"))
  cat(paste0("  - Max Tree Depth: ", 16, "\n"))
  cat(paste0("  - Output File: ", output_file, "\n"))
  
  tryCatch({
    # Start fitting the model
    cat("Fitting the BRM model...\n")
    proposed <- brm(
      count ~ s(yearn, by = state) + state + offset(log(population)),
      data = bact,
      family = "negbinomial",
      save_pars = save_pars(all = TRUE),
      chains = 4,
      iter = 2501, # iterations may need to be tuned for each model, but needs to be sufficiently high to ensure convergence I tested it out and 2501 worked well for most things. Could we add a way for this to be changeable when specifying parameters for running the Nextflow model 
      cores = 9,
      seed = 47,
      control = list(adapt_delta = 0.999, max_treedepth = 19)
    )
    
    # Save the model and log success
    saveRDS(proposed, output_file)
    cat("Successfully saved BRM model to: ", output_file, "\n")
  }, error = function(e) {
    # Log error if model fitting or saving fails
    cat("Failed to run or save BRM model: ", output_file, "\n")
    cat("Error: ", e$message, "\n")
  })
  
  # Return the model object
  cat("Returning the fitted BRM model.\n")
  return(proposed)
}

##################################################################
# LINPRED
##################################################################
LINPREAD_DRAW_FN <- function(data, model) {
  # Ensure the model is valid
  if (!inherits(model, "brmsfit")) {
    stop("Error: `model` must be a valid `brmsfit` object.")
  }
  
  # Generate posterior predictive draws
  d.prop.fdraws <- add_linpred_draws(
    object = model,  # Bayesian model
    newdata = data,  # Data for predictions
    ndraws = NULL,   # Use all posterior draws
    transform = TRUE, # Transform predictions to response (if you want it on the link scale set to FALSE)
    value = ".linpred" # Column for predicted values
  )
  
  # Ensure .draw column exists
  if (!".draw" %in% colnames(d.prop.fdraws)) {
    stop("Error: `.draw` column is missing in posterior draws.")
  }
  
  # Add predicted incidence (pinc) and actual incidence (inc)
  d.prop.fdraws$pinc <- (d.prop.fdraws$.linpred) / (d.prop.fdraws$population / 100000)
  d.prop.fdraws$inc <- (d.prop.fdraws$count) / (d.prop.fdraws$population / 100000)
  
  # Generate unique identifier for state and year
  d.prop.fdraws$code <- paste(d.prop.fdraws$state, d.prop.fdraws$yearn)
  
  return(d.prop.fdraws)
}
##################################################################
# Catchment
##################################################################

SAVE_PLOT <- function(p, label, pathogen_name) {
  # Include pathogen name in the file name
  save_file <- paste0(outBase, pathogen_name, "_", label, ".png")
  tryCatch({
    ggsave(save_file, p)
    message("Successfully saved plot: ", save_file)
  }, error = function(e) {
    message("Failed to save plot: ", save_file)
    message("Error: ", e$message)
  })
}

CATCHMENT <- function(prop_linpred) {
  # Check if .draw exists in the input
  if (!".draw" %in% colnames(prop_linpred)) {
    stop("Error: `.draw` column is missing in input to CATCHMENT.")
  }
  
  # Compute .value
  value <- prop_linpred %>%
    reshape2::dcast(yearn + .draw ~ state, value.var = ".linpred") %>%
    mutate(.value = rowSums(across(-c(yearn, .draw)), na.rm = TRUE))
  
  # Compute count
  count <- prop_linpred %>%
    reshape2::dcast(yearn + .draw ~ state, value.var = "count") %>%
    mutate(count = rowSums(across(-c(yearn, .draw)), na.rm = TRUE))
  
  # Compute population
  pop <- prop_linpred %>%
    reshape2::dcast(yearn + .draw ~ state, value.var = "population") %>%
    mutate(population = rowSums(across(-c(yearn, .draw)), na.rm = TRUE))
  
  # Merge results
  catch <- value %>%
    select(yearn, .draw, .value) %>%  # Include .draw explicitly
    left_join(pop %>% select(yearn, .draw, population), by = c("yearn", ".draw")) %>%
    left_join(count %>% select(yearn, .draw, count), by = c("yearn", ".draw"))
  
  # Ensure .draw remains intact
  if (!".draw" %in% colnames(catch)) {
    stop("Error: `.draw` column was dropped during CATCHMENT processing.")
  }
  
  # Save plots with pathogen name
  p <- ggplot(value, aes(x = .value, color = "red")) +
    geom_density() +
    geom_vline(data = count, aes(xintercept = count), color = "black") +
    facet_wrap(~yearn)
  SAVE_PLOT(p, "catch_value", pathogen_name)
  
  p <- ggplot(catch, aes(x = .value, color = "red")) +
    geom_density() +
    geom_vline(data = count, aes(xintercept = count), color = "black") +
    facet_wrap(~yearn)
  SAVE_PLOT(p, "catch", pathogen_name)
  
  p <- ggplot(catch) + geom_boxplot(aes(group = yearn, x = yearn, y = .value), color = "black") +
    geom_point(aes(x = yearn, y = count), color = "red") +
    theme_minimal()
  SAVE_PLOT(p, "catch_boxplot", pathogen_name)
  
  p <- ggplot(catch) +
    geom_violin(aes(group = yearn, x = yearn, y = .value), color = "black") +
    geom_point(aes(x = yearn, y = count), color = "red") +
    theme_minimal()
  SAVE_PLOT(p, "catch_vplot", pathogen_name)
  
  catch$inc <- (catch$count / (catch$population / 100000))
  catch$pinc <- (catch$.value / (catch$population / 100000))
  
  p <- ggplot(catch) + geom_boxplot(aes(group = yearn, x = yearn, y = pinc), color = "black") +
    geom_point(aes(x = yearn, y = inc), color = "red") +
    theme_minimal()
  SAVE_PLOT(p, "catch2_boxplot", pathogen_name)
  
  p <- ggplot(catch) + geom_violin(aes(group = yearn, x = yearn, y = pinc), color = "black") +
    geom_point(aes(x = yearn, y = inc), color = "red") +
    theme_minimal()
  SAVE_PLOT(p, "catch2_vplot", pathogen_name)
  
  # Return the final dataset with .draw intact
  return(catch)
}

##################################################################
# Catchment IR
##################################################################
LINPRED_TO_CATCHIR <- function(catchments) {
  # Check for required columns
  if (!all(c("count", "population", ".value", "yearn", ".draw") %in% colnames(catchments))) {
    stop("Error: Required columns (count, population, .value, yearn, .draw) are missing in `catchments`.")
  }
  
  # Summarize data with LL, median, and UL, keeping .draw
  raw_data <- catchments %>%
    group_by(yearn, .draw) %>%  # Keep .draw for posterior tracking
    summarise(
      count = mean(count, na.rm = TRUE),
      population = mean(population, na.rm = TRUE),
      .value = mean(.value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pinc = (.value) / (population / 100000),  # Predicted incidence
      inc = (count) / (population / 100000)    # Observed incidence
    )
  
  # Calculate summary statistics (LL, median, UL) # this
  summary_stats <- raw_data %>%
    group_by(yearn) %>%
    summarise(
      count_LL = quantile(count, probs = 0.025, na.rm = TRUE),
      count_median = median(count, na.rm = TRUE),
      count_UL = quantile(count, probs = 0.975, na.rm = TRUE),
      population_LL = quantile(population, probs = 0.025, na.rm = TRUE),
      population_median = median(population, na.rm = TRUE),
      population_UL = quantile(population, probs = 0.975, na.rm = TRUE),
      .value_LL = quantile(.value, probs = 0.025, na.rm = TRUE),
      .value_median = median(.value, na.rm = TRUE),
      .value_UL = quantile(.value, probs = 0.975, na.rm = TRUE),
      pinc_LL = quantile(pinc, probs = 0.025, na.rm = TRUE),
      pinc_median = median(pinc, na.rm = TRUE),
      pinc_UL = quantile(pinc, probs = 0.975, na.rm = TRUE),
      inc_LL = quantile(inc, probs = 0.025, na.rm = TRUE),
      inc_median = median(inc, na.rm = TRUE),
      inc_UL = quantile(inc, probs = 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Join summary statistics back to raw_data
  combined_data <- raw_data %>%
    left_join(summary_stats, by = "yearn")
  
  # Return combined data
  return(combined_data)
}
##################################################################
# IR COMP
##################################################################
IR_COMP <- function(catch, year1, year2, fileBase, target, travel, culture) {
  ref <- catch %>%
    filter(yearn <= year2 & yearn >= year1) %>%
    group_by(.draw) %>%
    summarise(across(c("count", "population", ".value"), list(mean = mean)))
  colnames(ref) <- c(".draw", "ref_count", "ref_pop", "ref_value")
  
  comb <- left_join(catch %>% select(yearn, .draw, .value, population, count), ref, by = c(".draw")) %>%
    mutate(
      pinc = (.value) / (population / 100000),
      inc = (count) / (population / 100000),
      ref.ir = (ref_count) / (ref_pop / 100000),
      ref.est.ir = (ref_value) / (ref_pop / 100000),
      rr = round(pinc / ref.est.ir, 2),
      pct.change = round(((pinc - ref.est.ir) / ref.est.ir) * 100, 2)
    )
  
  comb2 <- comb %>%
    group_by(yearn) %>%
    summarise(across(c("pinc", "inc", "ref.ir", "ref.est.ir", "ref_pop", "rr", "pct.change", ".value", "ref_value", "population", "ref_pop", "count", "ref_count"),
                     list(LL = ~quantile(., probs = 0.025, na.rm = TRUE),
                          median = median,
                          UL = ~quantile(., probs = 0.975, na.rm = TRUE)))) %>%
    select(yearn, .value_median, count_median, ref_value_median, ref_count_median, population_median, ref_pop_median, pinc_median, pinc_LL, pinc_UL, inc_median, ref.est.ir_median, ref.est.ir_LL, ref.est.ir_UL, ref.ir_median, rr_median, rr_LL, rr_UL, pct.change_median, pct.change_LL, pct.change_UL) %>%
    dplyr::rename(inc = inc_median)
  
  if (nrow(comb2) > 0) {
    comb2 <- comb2 %>%
      mutate(
        pathogen = target,
        travel = travel,
        culture = culture,
        startyear = year2
      )
  }
  
  saveFile <- paste0(outBase, fileBase)
  
  # Check for duplicates
  if (file.exists(saveFile)) {
    tryCatch({
      existing_data <- read.csv(saveFile, stringsAsFactors = FALSE)
      existing_data <- existing_data %>% select(names(comb2)) # Align columns
      comb2 <- anti_join(comb2, existing_data, by = colnames(comb2))
    }, error = function(e) {
      message("Error reading existing file, skipping duplicate check: ", e$message)
    })
  }
  
  # Write the data
  write.table(comb2, saveFile, append = TRUE, quote = TRUE, sep = ",", qmethod = "double", col.names = !file.exists(saveFile), row.names = FALSE)
  message("Results saved to: ", saveFile)
}