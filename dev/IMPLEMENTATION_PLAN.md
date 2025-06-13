# Implementation Plan: New Features

## Overview
This document outlines the implementation plan for two new configurable features requested by the client.

**Note**: Configuration files use CSV format instead of JSON because the current Singularity container does not include JSON parsing packages (jsonlite, rjson). CSV format leverages existing readr/tidyverse capabilities without requiring container modifications.

## Feature 1: Configurable Serotype Recoding

### Current State
- **Location**: `bin/preprocess.R` lines 56-62
- **Method**: Hardcoded list in `seroList` variable
- **Logic**: 
  ```r
  seroList <- c("NOT SPECIATED", "UNKNOWN", "PARTIAL SERO", "NOT SERO", "")
  mmwrdata$sero2 <- ifelse(mmwrdata$sero1 %in% seroList, "Missing", mmwrdata$sero1)
  mmwrdata$sero2 <- ifelse(grepl("UNDET", mmwrdata$sero2), "Missing", mmwrdata$sero2)
  ```

### Implementation Plan

#### 1. Add Command Line Parameter
- **Parameter**: `--serotype-config` (optional)
- **Type**: File path to CSV configuration
- **Default**: Use existing hardcoded behavior if not provided

#### 2. Configuration File Format (CSV)
```csv
serotype_value,replacement,pathogen,match_type,notes
NOT SPECIATED,Missing,all,exact,Default missing value
UNKNOWN,Missing,all,exact,Default missing value
PARTIAL SERO,Missing,all,exact,Default missing value
NOT SERO,Missing,all,exact,Default missing value
"",Missing,all,exact,Empty string missing value
UNDET,Missing,all,contains,Pattern-based rule
ROUGH,Missing,SALMONELLA,exact,Salmonella-specific
NONTYPEABLE,Missing,SALMONELLA,exact,Salmonella-specific
TYPHIMURIUM VAR 5-,TYPHIMURIUM,SALMONELLA,exact,Custom mapping
```

**Column definitions:**
- `serotype_value`: Original serotype value to match
- `replacement`: What to replace it with ("Missing" for non-informative values)
- `pathogen`: Apply rule to specific pathogen or "all" for universal rules
- `match_type`: "exact" for exact match, "contains" for substring matching
- `notes`: Optional documentation (ignored by code)

#### 3. Code Changes
1. **Add parameter parsing** in preprocess.R argument parser
2. **Create configuration reader function**:
   ```r
   read_serotype_config <- function(config_path = NULL) {
     if (is.null(config_path)) {
       # Return default configuration as data frame
       return(data.frame(
         serotype_value = c("NOT SPECIATED", "UNKNOWN", "PARTIAL SERO", "NOT SERO", "", "UNDET"),
         replacement = c("Missing", "Missing", "Missing", "Missing", "Missing", "Missing"),
         pathogen = rep("all", 6),
         match_type = c("exact", "exact", "exact", "exact", "exact", "contains"),
         stringsAsFactors = FALSE
       ))
     }
     # Read and validate CSV
     config <- read.csv(config_path, stringsAsFactors = FALSE)
     validate_serotype_config(config)
     return(config)
   }
   
   validate_serotype_config <- function(config) {
     required_cols <- c("serotype_value", "replacement", "pathogen", "match_type")
     missing_cols <- setdiff(required_cols, names(config))
     if (length(missing_cols) > 0) {
       stop("Missing required columns in serotype config: ", paste(missing_cols, collapse = ", "))
     }
     # Validate match_type values
     valid_match_types <- c("exact", "contains")
     invalid_types <- setdiff(config$match_type, valid_match_types)
     if (length(invalid_types) > 0) {
       stop("Invalid match_type values: ", paste(invalid_types, collapse = ", "))
     }
   }
   ```
3. **Replace hardcoded logic** with dynamic configuration application:
   ```r
   apply_serotype_config <- function(data, config, current_pathogen = NULL) {
     # Filter config for current pathogen
     if (!is.null(current_pathogen)) {
       applicable_config <- config[
         config$pathogen == "all" | config$pathogen == current_pathogen, 
       ]
     } else {
       applicable_config <- config
     }
     
     # Apply exact matches first
     exact_rules <- applicable_config[applicable_config$match_type == "exact", ]
     for (i in 1:nrow(exact_rules)) {
       data$sero2[data$sero1 == exact_rules$serotype_value[i]] <- exact_rules$replacement[i]
     }
     
     # Apply pattern matches
     pattern_rules <- applicable_config[applicable_config$match_type == "contains", ]
     for (i in 1:nrow(pattern_rules)) {
       data$sero2[grepl(pattern_rules$serotype_value[i], data$sero2)] <- pattern_rules$replacement[i]
     }
     
     return(data)
   }
   ```
4. **Update Nextflow pipeline** to pass parameter through

#### 4. Testing Strategy
- Create sample config files for edge cases
- Verify backward compatibility when no config provided
- Test pathogen-specific rules
- Validate CSV parsing and error handling
- Test both exact and pattern matching rules

---

## Feature 2: Configurable Catchment Definitions

### Current State
- **Location**: `bin/functions.R` in PATH_ANALYSIS, CYCLOSPORA_ANALYSIS, SALMONELLA_ANALYSIS
- **Method**: Hardcoded subset() calls
- **Logic**:
  ```r
  subset((state=="CA") | (state=="CO" & year>=2001) | (state=="CT") | 
         (state=="GA") | (state=="MD" & year>=1998) | (state=="MN") | 
         (state=="NM" & year>=2004) | (state=="NY" & year>=1998) | 
         (state=="OR") | (state=="TN" & year>=2000))
  ```

### Implementation Plan

#### 1. Add Function Parameter
- **Parameter**: `catchment_config` (optional)
- **Type**: File path to CSV configuration
- **Default**: Use existing hardcoded FoodNet definitions

#### 2. Configuration File Format (CSV)
```csv
state,start_year,end_year,pathogen_type,notes
CA,1996,2024,both,Original FoodNet site
CO,2001,2024,both,Joined 2001
CT,1996,2024,both,Original FoodNet site
GA,1996,2024,both,Original FoodNet site
MD,1998,2024,both,Joined 1998
MN,1996,2024,both,Original FoodNet site
NM,2004,2024,both,Joined 2004
NY,1998,2024,both,Joined 1998
OR,1996,2024,both,Original FoodNet site
TN,2000,2024,both,Joined 2000
```

#### 3. Code Changes

##### A. Add Configuration Reader Function
```r
read_catchment_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    # Return default FoodNet configuration
    return(data.frame(
      state = c("CA", "CO", "CT", "GA", "MD", "MN", "NM", "NY", "OR", "TN"),
      start_year = c(1996, 2001, 1996, 1996, 1998, 1996, 2004, 1998, 1996, 2000),
      end_year = rep(2024, 10),
      pathogen_type = rep("both", 10)
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
  # Additional validation...
}
```

##### B. Create Dynamic Filtering Function
```r
apply_catchment_filter <- function(data, catchment_config, pathogen_type = "both") {
  # Filter config for relevant pathogen type
  if (pathogen_type != "both") {
    relevant_config <- catchment_config[
      catchment_config$pathogen_type %in% c("both", pathogen_type), 
    ]
  } else {
    relevant_config <- catchment_config
  }
  
  # Build dynamic filter
  valid_combinations <- data.frame()
  for (i in 1:nrow(relevant_config)) {
    state_data <- data[data$state == relevant_config$state[i], ]
    year_filtered <- state_data[
      state_data$year >= relevant_config$start_year[i] & 
      state_data$year <= relevant_config$end_year[i], 
    ]
    valid_combinations <- rbind(valid_combinations, year_filtered)
  }
  
  return(valid_combinations)
}
```

##### C. Update Analysis Functions
```r
# Before (hardcoded):
selectDf <- selectDf %>% 
  subset((state=="CA") | (state=="CO" & year>=2001) | ...)

# After (configurable):
selectDf <- apply_catchment_filter(selectDf, catchment_config, "bacterial")
```

##### D. Update Function Signatures
```r
# Add catchment_config parameter to:
PATH_ANALYSIS <- function(mmwrdata, census, catchment_config = NULL)
CYCLOSPORA_ANALYSIS <- function(mmwrdata, census, catchment_config = NULL)
SALMONELLA_ANALYSIS <- function(mmwrdata, census, catchment_config = NULL)
```

##### E. Update trendy.R
- Add `--catchment-config` parameter
- Pass configuration to analysis functions
- Update Nextflow pipeline integration

#### 4. Testing Strategy
- Test with different state/year combinations
- Verify bacterial vs parasitic pathogen handling
- Test edge cases (missing states, invalid years)
- Ensure backward compatibility
- Validate CSV parsing and error handling

---

## Implementation Priority and Timeline

### Phase 1: Serotype Configuration (Easier)
- **Estimated time**: 2-3 hours
- **Risk level**: Low
- **Dependencies**: None

### Phase 2: Catchment Configuration (More Complex)
- **Estimated time**: 4-6 hours  
- **Risk level**: Medium
- **Dependencies**: Understanding of business logic
- **Recommendation**: Schedule Teams session with client first

### Phase 3: Testing and Documentation
- **Estimated time**: 2-3 hours
- **Deliverables**: 
  - Sample configuration files
  - Updated README documentation
  - Parameter documentation in help text

---

## Risk Mitigation

### Serotype Feature
- **Risk**: Breaking existing serotype logic
- **Mitigation**: Maintain exact backward compatibility when no config provided

### Catchment Feature  
- **Risk**: Complex business logic edge cases
- **Mitigation**: Teams session with client to clarify requirements
- **Risk**: Different pathogen type handling
- **Mitigation**: Thorough testing with both bacterial and parasitic datasets

### Both Features
- **Risk**: Pipeline integration complexity
- **Mitigation**: Test with existing Nextflow infrastructure
- **Risk**: Parameter validation and error handling
- **Mitigation**: Comprehensive input validation functions

---

## Success Criteria

1. **Backward Compatibility**: Pipeline works identically when no config files provided
2. **Configuration Flexibility**: Users can define custom serotype and catchment rules
3. **Error Handling**: Clear error messages for invalid configurations
4. **Documentation**: Updated help text and example config files
5. **Testing**: All existing functionality verified post-implementation