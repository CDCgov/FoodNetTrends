# Data Configuration Options for FoodNet Trends Pipeline

**Important**: These are CSV data configuration files, NOT Nextflow configuration files.

The FoodNet Trends pipeline now supports optional CSV configuration files to customize serotype recoding rules and catchment area definitions. This document explains how to use these features.

## Table of Contents
- [Serotype Configuration](#serotype-configuration)
- [Catchment Configuration](#catchment-configuration)
- [Usage Examples](#usage-examples)

## Serotype Configuration

The serotype configuration feature allows you to customize which serotype values are recoded as "Missing" or mapped to other values.

### Default Behavior

By default, the following serotype values are recoded to "Missing":
- NOT SPECIATED
- UNKNOWN
- PARTIAL SERO
- NOT SERO
- (empty string)
- Any value containing "UNDET"

### Custom Configuration

Create a CSV file with the following columns:

| Column | Description | Required |
|--------|-------------|----------|
| serotype_value | The original serotype value to match | Yes |
| replacement | What to replace it with | Yes |
| pathogen | Apply to specific pathogen or "all" | Yes |
| match_type | "exact" or "contains" | Yes |
| notes | Optional documentation | No |

Example configuration file (`analysis_configs/examples/serotype_config.csv`):

```csv
serotype_value,replacement,pathogen,match_type,notes
NOT SPECIATED,Missing,all,exact,Default missing value
UNKNOWN,Missing,all,exact,Default missing value  
PARTIAL SERO,Missing,all,exact,Default missing value
NOT SERO,Missing,all,exact,Default missing value
"",Missing,all,exact,Empty string missing value
UNDET,Missing,all,contains,Pattern-based rule for undetermined
ROUGH,Missing,SALMONELLA,exact,Salmonella-specific non-informative value
NONTYPEABLE,Missing,SALMONELLA,exact,Salmonella-specific non-informative value
```

### Using the Configuration

#### With Nextflow
```bash
nextflow run main.nf \
  --mmwrFile data/mmwr.sas7bdat \
  --censusFileB data/census_b.sas7bdat \
  --censusFileP data/census_p.sas7bdat \
  --serotype_config config/my_serotype_rules.csv
```

#### With R Script Directly
```bash
Rscript bin/preprocess.R \
  --mmwrFile data/mmwr.sas7bdat \
  --outputFile clean_mmwr.csv \
  --serotype-config config/my_serotype_rules.csv
```

## Catchment Configuration

The catchment configuration feature allows you to define custom geographic and temporal boundaries for analysis.

### Default Behavior

By default, the pipeline uses the standard FoodNet catchment areas:
- CA: 1996-present
- CO: 2001-present
- CT: 1996-present
- GA: 1996-present
- MD: 1998-present
- MN: 1996-present
- NM: 2004-present
- NY: 1998-present
- OR: 1996-present
- TN: 2000-present

### Custom Configuration

Create a CSV file with the following columns:

| Column | Description | Required |
|--------|-------------|----------|
| state | Two-letter state code | Yes |
| start_year | First year to include | Yes |
| end_year | Last year to include (use 9999 for no end) | Yes |
| pathogen_type | "bacterial", "parasitic", or "both" | No |
| notes | Optional documentation | No |

Example configuration file (`analysis_configs/examples/catchment_config.csv`):

```csv
state,start_year,end_year,pathogen_type,notes
CA,1996,9999,both,Original FoodNet site
CO,2001,9999,both,Joined 2001
CT,1996,9999,both,Original FoodNet site
GA,1996,9999,both,Original FoodNet site
MD,1998,9999,both,Joined 1998
MN,1996,9999,both,Original FoodNet site
NM,2004,9999,both,Joined 2004
NY,1998,9999,both,Joined 1998
OR,1996,9999,both,Original FoodNet site
TN,2000,9999,both,Joined 2000
```

### Using the Configuration

#### With Nextflow
```bash
nextflow run main.nf \
  --mmwrFile data/mmwr.sas7bdat \
  --censusFileB data/census_b.sas7bdat \
  --censusFileP data/census_p.sas7bdat \
  --catchment_config config/my_catchment_areas.csv
```

#### With R Script Directly
```bash
Rscript bin/trendy.R \
  --mmwrFile data/mmwr.sas7bdat \
  --censusFileB data/census_b.sas7bdat \
  --censusFileP data/census_p.sas7bdat \
  --projID myproject \
  --catchment-config config/my_catchment_areas.csv
```

## Usage Examples

### Example 1: Custom Analysis for Research Study

If you're conducting a study focusing on specific states and years:

1. Create a custom catchment configuration:
```csv
state,start_year,end_year,pathogen_type,notes
CA,2010,2020,both,Study period California
NY,2010,2020,both,Study period New York
TX,2015,2020,both,Added Texas for comparison
```

2. Run the pipeline:
```bash
nextflow run main.nf \
  --mmwrFile data/mmwr.sas7bdat \
  --censusFileB data/census_b.sas7bdat \
  --censusFileP data/census_p.sas7bdat \
  --catchment_config config/study_catchment.csv
```

### Example 2: Pathogen-Specific Serotype Rules

For a Salmonella-focused analysis with custom serotype mapping:

1. Create serotype configuration:
```csv
serotype_value,replacement,pathogen,match_type,notes
TYPHIMURIUM VAR 5-,TYPHIMURIUM,SALMONELLA,exact,Consolidate variants
I 4,[5],12:i:-,TYPHIMURIUM,SALMONELLA,exact,Monophasic variant
ROUGH,Missing,SALMONELLA,exact,Non-typeable
```

2. Run the pipeline:
```bash
nextflow run main.nf \
  --mmwrFile data/mmwr.sas7bdat \
  --censusFileB data/census_b.sas7bdat \
  --censusFileP data/census_p.sas7bdat \
  --serotype_config config/salmonella_serotypes.csv \
  --pathogen SALMONELLA
```

### Example 3: Both Configurations Together

```bash
nextflow run main.nf \
  --mmwrFile data/mmwr.sas7bdat \
  --censusFileB data/census_b.sas7bdat \
  --censusFileP data/census_p.sas7bdat \
  --serotype_config config/custom_serotypes.csv \
  --catchment_config config/custom_catchment.csv \
  --projID custom_analysis_2024
```

## Best Practices

1. **Test with Small Datasets**: When using custom configurations, test with a subset of data first.

2. **Document Your Rules**: Use the notes column to explain why certain rules exist.

3. **Version Control**: Keep your configuration files in version control alongside your analysis code.

4. **Validation**: The pipeline validates configuration files and will provide clear error messages if issues are found.

5. **Backward Compatibility**: If no configuration files are provided, the pipeline uses the original hardcoded defaults.