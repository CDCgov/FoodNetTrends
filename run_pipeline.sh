#!/bin/bash
#==============================================================================
# FoodNetTrends Pipeline v1.0.0-rc.1 - Interactive Execution Script
#==============================================================================
#
# OVERVIEW FOR MAINTAINERS:
# This is the main user interface for the FoodNetTrends analysis pipeline.
# It provides an interactive menu system that guides epidemiologists through
# parameter selection and executes the full Bayesian spline analysis workflow.
#
# KEY RESPONSIBILITIES:
# 1. Data Discovery: Automatically finds MMWR and census files
# 2. Parameter Collection: Interactive menus for pathogen/state/serotype selection
# 3. Configuration Management: Save/load analysis configurations
# 4. Resource Optimization: HPC profile selection for different data scales
# 5. Pipeline Execution: Coordinates Nextflow workflow with proper parameters
# 6. Quality Validation: Checks outputs and validates spline trend generation
# 7. Progress Tracking: Real-time monitoring with shell/R coordination
#
# DESIGN PHILOSOPHY:
# - User-friendly: No command-line expertise required
# - Data-driven: Discovers available options from actual data files
# - Production-ready: Comprehensive error handling and validation
# - HPC-optimized: Resource profiles for different computational requirements
#
# If --get-command-only is provided, just print the command that would be run
if [[ "$1" == "--get-command-only" ]]; then
    GET_COMMAND_ONLY=true
    shift
else
    GET_COMMAND_ONLY=false
fi
#
# Features:
#   - Interactive parameter collection and validation
#   - HPC-optimized resource profiles for different analysis scales
#   - Automatic file discovery and validation
#   - Memory and CPU optimization for Bayesian modeling
#   - Support for preprocessing, full analysis, and dashboard-only modes
#
# Usage:
#   ./run_pipeline.sh
#
# Dependencies:
#   - Nextflow ≥ 24.10.4
#   - Singularity ≥ 4.1.4
#   - Access to HPC environment with SGE scheduler
#   - Required data files (MMWR and census data)
#
# Last updated: 2025-05-22
#==============================================================================

# Initialize log file for error tracking
error_log="foodnet_errors.log"
echo "$(date): Starting FoodNetTrends Analysis Pipeline v1.0.0-rc.1" > "$error_log"

# =============================================================================
# DYNAMIC PROGRESS TRACKING SYSTEM
# =============================================================================
#
# MAINTAINER NOTE: This section implements a sophisticated progress tracking
# system that coordinates between the shell script and R analysis components.
#
# PROGRESS COMPONENTS:
# 1. Shell-based progress bars: Visual feedback during file operations
# 2. R coordination: Signals between shell and R scripts for MCMC tracking
# 3. Stage completion: Milestone-based progress for complex workflows
# 4. Error integration: Progress updates even when components fail
#
# DESIGN FEATURES:
# - In-place updates: Uses \r to overwrite lines for clean display
# - Cross-platform: Works on both HPC and local environments
# - Fail-safe: Continues working even if progress.R is unavailable
# - User-friendly: Clear visual feedback for long-running Bayesian analyses
# ============================================================================="

# Dynamic Progress Bar Function
# MAINTAINER NOTE: Creates visual progress bars that update in-place using \r
# Provides real-time feedback during long-running operations like file discovery
# and Bayesian model fitting. Essential for user experience on HPC systems.
show_pipeline_progress() {
    local current=$1
    local total=$2
    local stage_name="$3"
    local message="$4"
    local width=50
    
    # Calculate percentage and bar length
    local percent=$((current * 100 / total))
    local filled=$((current * width / total))
    
    # Build progress bar with Unicode blocks
    local bar=""
    for ((i=0; i<filled; i++)); do bar+="█"; done
    for ((i=filled; i<width; i++)); do bar+="░"; done
    
    # Create progress line
    local progress_line
    if [[ -n "$message" ]]; then
        progress_line=$(printf "\r[%s] %3d%% | %s - %s" "$bar" "$percent" "$stage_name" "$message")
    else
        progress_line=$(printf "\r[%s] %3d%% | %s" "$bar" "$percent" "$stage_name")
    fi
    
    # Display on terminal with proper output handling
    printf "%s" "$progress_line" >/dev/tty 2>/dev/null || printf "%s" "$progress_line"
}

# Clear progress bar 
clear_pipeline_progress() {
    printf "\r%*s\r" 100 "" >/dev/tty 2>/dev/null || printf "\r%*s\r" 100 ""
}

# Complete a progress stage
complete_pipeline_stage() {
    local stage_name="$1"
    clear_pipeline_progress
    echo "✓ $stage_name completed" >/dev/tty 2>/dev/null || echo "✓ $stage_name completed"
}

# Monitor file operations with progress
monitor_file_operation() {
    local operation_name="$1"
    local check_function="$2"
    local max_wait="${3:-300}"  # Default 5 minutes
    
    echo "Starting $operation_name..." >/dev/tty 2>/dev/null || echo "Starting $operation_name..."
    
    local elapsed=0
    local dots=0
    
    while ! eval "$check_function" && [[ $elapsed -lt $max_wait ]]; do
        dots=$(( (elapsed / 2) % 4 ))
        local dot_string=""
        for ((i=0; i<dots; i++)); do dot_string+="."; done
        
        show_pipeline_progress $elapsed $max_wait "$operation_name" "$dot_string"
        sleep 2
        elapsed=$((elapsed + 2))
    done
    
    if eval "$check_function"; then
        complete_pipeline_stage "$operation_name"
        return 0
    else
        clear_pipeline_progress
        echo "✗ $operation_name timed out after ${max_wait}s" >/dev/tty 2>/dev/null || echo "✗ $operation_name timed out"
        return 1
    fi
}

# Execute pipeline stages with coordinated progress
execute_pipeline_stages() {
    local stages=("$@")
    local total=${#stages[@]}
    local current=0
    
    echo ""
    echo "Executing FoodNetTrends Pipeline (${total} stages):" >/dev/tty 2>/dev/null || echo "Executing FoodNetTrends Pipeline (${total} stages):"
    echo ""
    
    for stage_info in "${stages[@]}"; do
        local stage_name="${stage_info%%:*}"
        local stage_function="${stage_info##*:}"
        current=$((current + 1))
        
        show_pipeline_progress $current $total "Stage $current" "$stage_name"
        
        # Execute stage function
        if declare -f "$stage_function" >/dev/null; then
            "$stage_function"
            local exit_code=$?
            if [[ $exit_code -ne 0 ]]; then
                clear_pipeline_progress
                echo "✗ Stage $current ($stage_name) failed with exit code $exit_code" >/dev/tty 2>/dev/null
                return $exit_code
            fi
        else
            echo "Warning: Stage function '$stage_function' not found" >/dev/tty 2>/dev/null
        fi
        
        sleep 0.5  # Brief pause to show completion
    done
    
    complete_pipeline_stage "All pipeline stages"
    echo ""
}

# =============================================================================
# OUTPUT VALIDATION FUNCTIONS
# =============================================================================

# Validate that expected outputs were generated
validate_analysis_outputs() {
    local pathogen="$1"
    local output_dir="$2"
    local missing_files=()
    local warning_files=()
    
    echo "Validating analysis outputs for $pathogen..."
    echo "Looking in directory: ${output_dir}/"
    echo "Directory contents:"
    ls -la "${output_dir}/" 2>/dev/null || echo "  Directory does not exist"
    
    # Check for critical output files
    local expected_files=(
        "${pathogen}_brm.Rds"
        "${pathogen}_IRCatch.csv"
        "${pathogen}_summary.txt"
    )
    
    # Check for visualization files (should exist unless error)
    local viz_files=(
        "${pathogen}_spline_trend.png"
        "${pathogen}_state_spline_trends.png"
        "${pathogen}_foodnettrends_comparison.png"
    )
    
    # Check critical files - look in output directory
    for file in "${expected_files[@]}"; do
        local file_path="${output_dir}/${file}"
        if [[ ! -f "$file_path" ]]; then
            missing_files+=("$file")
        elif [[ ! -s "$file_path" ]]; then
            warning_files+=("$file (empty)")
        fi
    done
    
    # Check visualization files - look in output directory
    for file in "${viz_files[@]}"; do
        local file_path="${output_dir}/${file}"
        if [[ ! -f "$file_path" ]]; then
            # Check for error versions
            local error_file="${file%.*}_error.${file##*.}"
            if [[ -f "${output_dir}/${error_file}" ]]; then
                warning_files+=("$file (error plot generated)")
            else
                missing_files+=("$file")
            fi
        fi
    done
    
    # Report results
    if [[ ${#missing_files[@]} -eq 0 && ${#warning_files[@]} -eq 0 ]]; then
        echo "✓ All expected outputs generated successfully for $pathogen"
        return 0
    else
        echo "⚠ Output validation issues for $pathogen:"
        
        for file in "${missing_files[@]}"; do
            echo "  ✗ Missing: $file"
        done
        
        for file in "${warning_files[@]}"; do
            echo "  ⚠ Warning: $file"
        done
        
        return 1
    fi
}

# Check if spline trends were generated (main fix verification)
verify_spline_trends() {
    local pathogen="$1"
    local output_dir="$2"
    
    local spline_files=(
        "${pathogen}_spline_trend.png"
        "${pathogen}_foodnettrends_comparison.png"
    )
    
    local spline_generated=0
    for file in "${spline_files[@]}"; do
        if [[ -f "${output_dir}/${file}" ]]; then
            spline_generated=1
            break
        fi
    done
    
    if [[ $spline_generated -eq 1 ]]; then
        echo "✓ Spline trend visualizations generated - 'spikey' graph issue should be resolved"
        return 0
    else
        echo "✗ No spline trend plots found - spikey graph issue may persist"
        return 1
    fi
}

# =============================================================================
# CONFIGURATION MANAGEMENT
# =============================================================================

# Save current configuration to file
save_configuration() {
    local config_file="foodnet_config_$(date +%Y%m%d_%H%M%S).sh"
    local config_dir="./configs"
    
    # Create configs directory if it doesn't exist
    mkdir -p "$config_dir" 2>/dev/null
    local full_config_path="$config_dir/$config_file"
    
    cat > "$full_config_path" <<EOF
#!/bin/bash
# FoodNetTrends Configuration - Generated $(date)
# This file can be sourced to reproduce the same analysis configuration

# File paths
export mmwrFile="$mmwrFile"
export censusFileB="$censusFileB" 
export censusFileP="$censusFileP"
export preprocessed_metadata="$preprocessed_metadata"

# Analysis parameters
export pathogens="$pathogens"
export travel="$travel"
export cidt="$cidt"
export stec_serogroups="$stec_serogroups"
export salmonella_serotypes="$salmonella_serotypes"
export states="$states"

# Resource settings
export cores="$cores"
export memory="$memory"
export queue="$queue"
export runtime="$runtime"

# Output settings
export outDir="$outDir"
export projID="$projID"
export enable_dashboard="$enable_dashboard"
export background="$background"

# Workflow mode
export workflow_mode="$workflow_mode"

echo "Configuration loaded from: $full_config_path"
echo "Generated on: $(date)"
EOF

    chmod +x "$full_config_path"
    echo "✓ Configuration saved to: $full_config_path"
    echo "  To reuse: source $full_config_path && ./run_pipeline.sh"
}

# Display configuration summary
show_configuration_summary() {
    echo ""
    echo "======== CONFIGURATION SUMMARY ========="
    echo "Input Files:"
    echo "  MMWR Data: $mmwrFile"
    echo "  Census Bacterial: ${censusFileB:-'(not specified)'}"
    echo "  Census Parasitic: ${censusFileP:-'(not specified)'}"
    if [[ -n "$preprocessed_metadata" ]]; then
        echo "  Metadata: $preprocessed_metadata"
    fi
    echo ""
    echo "Analysis Settings:"
    echo "  Pathogens: $pathogens"
    echo "  Travel Filter: $travel"
    echo "  CIDT Filter: $cidt"
    if [[ -n "$stec_serogroups" ]]; then
        echo "  STEC Serogroups: $stec_serogroups"
    fi
    if [[ -n "$salmonella_serotypes" ]]; then
        echo "  Salmonella Serotypes: $salmonella_serotypes"
    fi
    echo "  States: ${states:-'ALL'}"
    echo ""
    echo "Resources:"
    echo "  Cores: $cores"
    echo "  Memory: ${memory}GB"
    echo "  Queue: $queue"
    echo "  Runtime: $runtime"
    echo ""
    echo "Output:"
    echo "  Directory: $outDir"
    echo "  Project ID: $projID"
    echo "  Dashboard: $enable_dashboard"
    echo "  Background: $background"
    echo "========================================="
    echo ""
}

# =============================================================================
# FILE VALIDATION FUNCTIONS
# =============================================================================

# Early validation of all required files
validate_required_files() {
    local missing_files=()
    local validation_errors=()
    
    echo "Validating input files..."
    show_pipeline_progress 1 4 "File Validation" "Checking MMWR file"
    
    # Check MMWR file
    if [[ -n "$mmwrFile" ]]; then
        if [[ ! -f "$mmwrFile" ]]; then
            missing_files+=("MMWR data file: $mmwrFile")
        elif [[ ! -r "$mmwrFile" ]]; then
            validation_errors+=("MMWR file not readable: $mmwrFile")
        fi
    else
        missing_files+=("MMWR data file (not specified)")
    fi
    
    show_pipeline_progress 2 4 "File Validation" "Checking census files"
    
    # Check census files
    if [[ -n "$censusFileB" ]]; then
        if [[ ! -f "$censusFileB" ]]; then
            missing_files+=("Bacterial census file: $censusFileB")
        elif [[ ! -r "$censusFileB" ]]; then
            validation_errors+=("Bacterial census file not readable: $censusFileB")
        fi
    fi
    
    if [[ -n "$censusFileP" ]]; then
        if [[ ! -f "$censusFileP" ]]; then
            missing_files+=("Parasitic census file: $censusFileP")
        elif [[ ! -r "$censusFileP" ]]; then
            validation_errors+=("Parasitic census file not readable: $censusFileP")
        fi
    fi
    
    show_pipeline_progress 3 4 "File Validation" "Checking metadata files"
    
    # Check metadata file if specified
    if [[ -n "$preprocessed_metadata" && ! -f "$preprocessed_metadata" ]]; then
        validation_errors+=("Metadata file not found: $preprocessed_metadata")
    fi
    
    show_pipeline_progress 4 4 "File Validation" "Validation complete"
    complete_pipeline_stage "File validation"
    
    # Report validation results
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        echo ""
        echo "✗ Missing required files:"
        for file in "${missing_files[@]}"; do
            echo "  - $file"
        done
        echo ""
        return 1
    fi
    
    if [[ ${#validation_errors[@]} -gt 0 ]]; then
        echo ""
        echo "✗ File validation errors:"
        for error in "${validation_errors[@]}"; do
            echo "  - $error"
        done
        echo ""
        return 1
    fi
    
    echo "✓ All required files validated successfully"
    return 0
}

# =============================================================================
# DATA DISCOVERY FUNCTIONS  
# =============================================================================

# Data-Driven Pathogen Discovery Function
# MAINTAINER NOTE: This function implements the core data discovery system
# that replaced hardcoded pathogen lists. It examines actual data files to
# determine which pathogens are available for analysis.
#
# DISCOVERY STRATEGY:
# 1. First tries metadata JSON file (fastest, if available)
# 2. Falls back to direct CSV analysis using R
# 3. Uses case-insensitive matching for robust pathogen detection
# 4. Returns space-separated list of available pathogens
#
# CRITICAL: This function is essential for user experience - it ensures
# users only see pathogens that actually exist in their data files.
discover_pathogens() {
    local mmwr_file="$1"
    local metadata_file="$2"
    
    show_pipeline_progress 1 3 "Data Discovery" "Analyzing pathogens"
    
    # Try metadata file first
    if [[ -f "$metadata_file" ]]; then
        local metadata_pathogens=$(parse_metadata_json "$metadata_file" "pathogens")
        if [[ -n "$metadata_pathogens" ]]; then
            echo "$metadata_pathogens"
            return 0
        fi
    fi
    
    # Fallback: scan MMWR file directly using R for robust CSV parsing
    if [[ -f "$mmwr_file" ]]; then
        show_pipeline_progress 2 3 "Data Discovery" "Scanning MMWR data"
        
        # Use R to properly parse CSV and find pathogen column
        local discovered_pathogens=$(Rscript -e "
            tryCatch({
                data <- read.csv('$mmwr_file', stringsAsFactors = FALSE, nrows = 1000)
                
                # Look for pathogen column (case insensitive)
                pathogen_col <- NULL
                for (col in names(data)) {
                    if (grepl('^pathogen$|^organism$|^etiology$', col, ignore.case = TRUE)) {
                        pathogen_col <- col
                        break
                    }
                }
                
                if (is.null(pathogen_col)) {
                    # Try first few columns that might contain pathogen data
                    for (i in 1:min(5, ncol(data))) {
                        col_values <- unique(toupper(as.character(data[[i]])))
                        col_values <- col_values[col_values != '' & !is.na(col_values)]
                        
                        # Check if this looks like pathogen data
                        known_pathogens <- c('SALMONELLA', 'CAMPYLOBACTER', 'SHIGA', 'STEC', 'ECOLI', 'CYCLOSPORA', 'LISTERIA', 'VIBRIO', 'YERSINIA')
                        matches <- sum(sapply(known_pathogens, function(p) any(grepl(p, col_values))))
                        
                        if (matches > 0) {
                            pathogen_col <- names(data)[i]
                            break
                        }
                    }
                }
                
                if (!is.null(pathogen_col)) {
                    pathogens <- unique(toupper(as.character(data[[pathogen_col]])))
                    pathogens <- pathogens[pathogens != '' & !is.na(pathogens)]
                    
                    # Clean up pathogen names - keep only alphanumeric and common symbols
                    pathogens <- gsub('[^A-Z0-9 _-]', '', pathogens)
                    pathogens <- pathogens[nchar(pathogens) > 0 & nchar(pathogens) < 50]
                    
                    cat(paste(pathogens, collapse=','))
                } else {
                    cat('')
                }
            }, error = function(e) {
                cat('')
            })
        " 2>/dev/null)
        
        if [[ -n "$discovered_pathogens" ]]; then
            echo "$discovered_pathogens"
            return 0
        fi
    fi
    
    # Emergency fallback - provide common FoodNet pathogens for manual selection
    show_pipeline_progress 3 3 "Data Discovery" "Using fallback pathogen list"
    echo "SALMONELLA,CAMPYLOBACTER,SHIGA,STEC,CYCLOSPORA,LISTERIA,VIBRIO,YERSINIA"
    return 1
}

# Discover available states from data files  
discover_states() {
    local mmwr_file="$1"
    local metadata_file="$2"
    
    show_pipeline_progress 1 3 "Data Discovery" "Analyzing states"
    
    # Try metadata file first
    if [[ -f "$metadata_file" ]]; then
        local metadata_states=$(parse_metadata_json "$metadata_file" "states")
        if [[ -n "$metadata_states" ]]; then
            echo "$metadata_states"
            return 0
        fi
    fi
    
    # Fallback: scan MMWR file for state column
    if [[ -f "$mmwr_file" ]]; then
        show_pipeline_progress 2 3 "Data Discovery" "Scanning state data"
        
        # Find state column (common names: state, STATE, st, ST)
        local state_col=$(head -1 "$mmwr_file" | tr ',' '\n' | grep -i -n -E '^state$|^st$' | cut -d: -f1 | head -1)
        
        if [[ -n "$state_col" ]]; then
            # Extract unique states from identified column
            local discovered_states=$(awk -F',' -v col="$state_col" 'NR>1 && $col!="" {states[toupper($col)]++} END {for(s in states) printf "%s,", s}' "$mmwr_file" | sed 's/,$//')
            
            if [[ -n "$discovered_states" ]]; then
                echo "$discovered_states"
                return 0
            fi
        fi
    fi
    
    # Emergency fallback
    echo ""
    return 1
}

# Discover available serotypes/serogroups for specific pathogens
discover_pathogen_subtypes() {
    local pathogen="$1"
    local mmwr_file="$2" 
    local metadata_file="$3"
    local subtype="$4"  # "serotype" or "serogroup"
    
    # Try metadata first
    if [[ -f "$metadata_file" ]]; then
        local metadata_key="${pathogen,,}_${subtype}s"  # e.g., "stec_serogroups"
        local metadata_subtypes=$(parse_metadata_json "$metadata_file" "$metadata_key")
        if [[ -n "$metadata_subtypes" ]]; then
            echo "$metadata_subtypes"
            return 0
        fi
    fi
    
    # Fallback: scan MMWR file using robust R-based CSV parsing
    if [[ -f "$mmwr_file" ]]; then
        local discovered_subtypes=$(Rscript -e "
            tryCatch({
                data <- read.csv('$mmwr_file', stringsAsFactors = FALSE, nrows = 2000)
                
                # Find pathogen column (same logic as main pathogen discovery)
                pathogen_col <- NULL
                for (col in names(data)) {
                    if (grepl('^pathogen\$|^organism\$|^etiology\$', col, ignore.case = TRUE)) {
                        pathogen_col <- col
                        break
                    }
                }
                
                if (is.null(pathogen_col)) {
                    # Try first few columns for pathogen data
                    for (i in 1:min(5, ncol(data))) {
                        col_values <- unique(toupper(as.character(data[[i]])))
                        col_values <- col_values[col_values != '' & !is.na(col_values)]
                        known_pathogens <- c('SALMONELLA', 'CAMPYLOBACTER', 'SHIGA', 'STEC', 'ECOLI', 'CYCLOSPORA', 'LISTERIA', 'VIBRIO', 'YERSINIA')
                        matches <- sum(sapply(known_pathogens, function(p) any(grepl(p, col_values))))
                        if (matches > 0) {
                            pathogen_col <- names(data)[i]
                            break
                        }
                    }
                }
                
                # Find subtype column (serotype/serogroup)
                subtype_col <- NULL
                if ('$subtype' == 'serotype') {
                    for (col in names(data)) {
                        if (grepl('serotype|sero(?!group)', col, ignore.case = TRUE, perl = TRUE)) {
                            subtype_col <- col
                            break
                        }
                    }
                } else if ('$subtype' == 'serogroup') {
                    for (col in names(data)) {
                        if (grepl('serogroup|sero_group', col, ignore.case = TRUE)) {
                            subtype_col <- col
                            break
                        }
                    }
                }
                
                if (!is.null(pathogen_col) && !is.null(subtype_col)) {
                    # Filter data for this pathogen
                    pathogen_data <- data[toupper(data[[pathogen_col]]) == toupper('$pathogen'), ]
                    
                    if (nrow(pathogen_data) > 0) {
                        # Extract unique subtypes
                        subtypes <- unique(as.character(pathogen_data[[subtype_col]]))
                        subtypes <- subtypes[subtypes != '' & !is.na(subtypes)]
                        
                        # Clean up subtypes - remove obvious garbage
                        subtypes <- gsub('[\"\\\\]', '', subtypes)  # Remove quotes and backslashes
                        subtypes <- subtypes[nchar(subtypes) > 0 & nchar(subtypes) < 100]
                        subtypes <- subtypes[!grepl('TRAVEL|FAMILY|PFGE|DAYCARE|ONSET|MOM|FATHER|DIALYSIS|UNKNOWN|^[0-9]+\$', subtypes)]
                        
                        if (length(subtypes) > 0) {
                            cat(paste(subtypes, collapse=','))
                        } else {
                            cat('')
                        }
                    } else {
                        cat('')
                    }
                } else {
                    cat('')
                }
            }, error = function(e) {
                cat('')
            })
        " 2>/dev/null)
        
        if [[ -n "$discovered_subtypes" ]]; then
            echo "$discovered_subtypes"
            return 0
        fi
    fi
    
    # Return empty if nothing found
    echo ""
    return 1
}

# JSON parsing function using Python
parse_metadata_json() {
    local json_file="$1"
    local field="$2"
    
    if [[ ! -f "$json_file" ]]; then
        echo ""
        return 1
    fi
    
    # Use Python for reliable JSON parsing
    if command -v python3 >/dev/null 2>&1; then
        # echo "DEBUG: Python available, parsing field '$field' from file '$json_file'" >&2
        
        local python_result=$(python3 -c "
import json
import sys

try:
    with open('$json_file', 'r') as f:
        data = json.load(f)
    
    field_name = '$field'
    # print(f'DEBUG: Parsing field: {field_name}', file=sys.stderr)
    
    # Handle specific serotype/serogroup fields
    if field_name == 'salmonella_serotypes':
        # Try _names field first (preferred)
        if 'salmonella_serotype_names' in data:
            # print('DEBUG: Using salmonella_serotype_names field', file=sys.stderr)
            serotype_names = data['salmonella_serotype_names']
            # Limit to first 5 for cleaner display
            if len(serotype_names) > 5:
                result = serotype_names[:5] + [f'...and {len(serotype_names) - 5} others']
            else:
                result = serotype_names
            print('|'.join(result))
        else:
            # print('DEBUG: No salmonella serotype data found', file=sys.stderr)
            print('')
    
    elif field_name == 'stec_serogroups':
        # Try _names field first (preferred)
        if 'stec_serogroup_names' in data:
            # print('DEBUG: Using stec_serogroup_names field', file=sys.stderr)
            serogroup_names = data['stec_serogroup_names']
            print('|'.join(serogroup_names))
        # Fallback to object keys
        elif 'stec_serogroups' in data and isinstance(data['stec_serogroups'], dict):
            # print('DEBUG: Using stec_serogroups object keys', file=sys.stderr)
            serogroup_names = list(data['stec_serogroups'].keys())
            print('|'.join(serogroup_names))
        else:
            # print('DEBUG: No STEC serogroup data found', file=sys.stderr)
            print('')
    
    # Generic field handling for other fields (pathogens, states, etc.)
    elif field_name in data:
        field_data = data[field_name]
        if isinstance(field_data, list):
            # Array field - join with commas
            print(','.join(str(x) for x in field_data))
        elif isinstance(field_data, dict):
            # Object field - return keys
            print(','.join(field_data.keys()))
        else:
            # Single value
            print(str(field_data))
    else:
        # print(f'DEBUG: Field {field_name} not found in metadata', file=sys.stderr)
        print('')

except Exception as e:
    # print(f'DEBUG: Python Error: {e}', file=sys.stderr)
    print('')
")
        
        if [[ -n "$python_result" ]]; then
            # echo "DEBUG: Python parsing succeeded, returning: '$python_result'" >&2
            echo "$python_result"
            return 0
        else
            # echo "DEBUG: Python parsing failed or returned empty" >&2
        fi
    else
        # echo "DEBUG: Python not available" >&2
    fi
    
    # Fallback: Simple bash parsing for basic JSON arrays
    if [[ "$field" == "pathogens" || "$field" == "states" ]]; then
        local bash_result=$(grep -o "\"$field\"[[:space:]]*:[[:space:]]*\[[^]]*\]" "$json_file" | \
                           sed 's/.*\[\(.*\)\].*/\1/' | \
                           sed 's/"//g' | \
                           tr -d ' ' | \
                           sed 's/,/, /g')
        
        if [[ -n "$bash_result" ]]; then
            echo "$bash_result"
            return 0
        fi
    fi
    
    # Enhanced bash fallback for serotype/serogroup objects
    if [[ "$field" == "salmonella_serotypes" || "$field" == "stec_serogroups" ]]; then
        # Try to extract from object keys first
        local bash_result=$(grep -o "\"$field\"[[:space:]]*:[[:space:]]*{[^}]*}" "$json_file" | \
                           sed 's/.*{\(.*\)}.*/\1/' | \
                           grep -o '"[^"]*"[[:space:]]*:' | \
                           sed 's/"//g' | sed 's/[[:space:]]*://' | \
                           tr '\n' ',' | sed 's/,$//')
        
        if [[ -n "$bash_result" ]]; then
            echo "$bash_result"
            return 0
        fi
        
        # Try alternative field names
        local alt_field="${field%s}_names"  # Convert "serotypes" to "serotype_names"
        local alt_result=$(grep -o "\"$alt_field\"[[:space:]]*:[[:space:]]*\[[^]]*\]" "$json_file" | \
                          sed 's/.*\[\(.*\)\].*/\1/' | \
                          sed 's/"//g' | \
                          tr -d ' ' | \
                          sed 's/,/, /g')
        
        if [[ -n "$alt_result" ]]; then
            echo "$alt_result"
            return 0
        fi
    fi
    
    # If all else fails, return empty
    echo ""
    return 1
}

# Process user-entered serotypes with flexible separators and case handling
process_user_serotypes() {
    local user_input="$1"
    local available_serotypes="$2"  # Pipe-separated list from metadata
    
    # Convert to uppercase
    user_input=$(echo "$user_input" | tr '[:lower:]' '[:upper:]')
    
    # Replace multiple separators (comma, pipe, space) with a single pipe
    # Also handle combinations like ", " or " | "
    user_input=$(echo "$user_input" | sed 's/[,|]/ /g' | tr -s ' ' | tr ' ' '|')
    
    # Convert available serotypes to array for validation
    IFS='|' read -ra AVAILABLE_ARRAY <<< "$available_serotypes"
    
    # Process user input
    IFS='|' read -ra USER_ARRAY <<< "$user_input"
    
    # Validate and build final list
    local valid_serotypes=()
    local invalid_serotypes=()
    
    for serotype in "${USER_ARRAY[@]}"; do
        serotype=$(echo "$serotype" | xargs)  # Trim whitespace
        if [[ -n "$serotype" ]]; then
            # Check if serotype exists in available list
            local found=0
            for available in "${AVAILABLE_ARRAY[@]}"; do
                if [[ "${available^^}" == "${serotype^^}" ]]; then
                    # Use the properly cased version from metadata
                    valid_serotypes+=("$available")
                    found=1
                    break
                fi
            done
            if [[ $found -eq 0 ]]; then
                invalid_serotypes+=("$serotype")
            fi
        fi
    done
    
    # Display results
    if [[ ${#valid_serotypes[@]} -gt 0 ]]; then
        echo ""
        echo "Selected serotypes:"
        echo "  $(IFS=' | '; echo "${valid_serotypes[*]}")"
    fi
    
    if [[ ${#invalid_serotypes[@]} -gt 0 ]]; then
        echo ""
        echo "⚠️  Warning: The following serotypes were not found in the dataset:"
        echo "  $(IFS=' | '; echo "${invalid_serotypes[*]}")"
    fi
    
    # Return comma-separated list for pipeline (maintaining compatibility)
    if [[ ${#valid_serotypes[@]} -gt 0 ]]; then
        IFS=',' 
        echo "${valid_serotypes[*]}"
    else
        echo ""
    fi
}

# Set up paths, files
# Make data directory path environment-aware
if [ -d "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/" ]; then
  DEFAULT_DATA_DIR="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/"
elif [ -d "/project/foodnet/data" ]; then
  DEFAULT_DATA_DIR="/project/foodnet/data"
else
  DEFAULT_DATA_DIR="./data"
  # Create local data directory if it doesn't exist
  mkdir -p "$DEFAULT_DATA_DIR"
fi
outDir="output"  # Default to "output" directory in current location

# Setup modules based on environment
if command -v module &> /dev/null; then
    module purge
    module load nextflow/24.10.4
    module load singularity/4.1.4
    module load java/17.0.6
else
    echo "Warning: Module system not detected. Dependencies must be in PATH."
    echo "$(date): Module system not detected" >> "$error_log"
fi

# Set up TMPDIR if not already defined
TMPDIR=${TMPDIR:-/scicomp/scratch/$(whoami)}
mkdir -p "$TMPDIR/nextflow" 2>/dev/null

# Create timestamp for ID
timestamp=$(date +%Y%m%d_%H%M%S)

# Display welcome banner
echo "========================================="
echo "   FoodNetTrends Analysis Pipeline      "
echo "         v1.0.0-rc.1                    "
echo "========================================="
echo ""

# Search for available metadata files first
echo "Searching for preprocessed data metadata files..."
mapfile -t found_json < <(find . -type f -name "*_metadata.json" 2>/dev/null | grep -E "(preprocessed|output)" | sort -r)

if [[ ${#found_json[@]} -gt 0 ]]; then
    echo ""
    echo "Found ${#found_json[@]} preprocessed dataset(s):"
    echo ""
    for i in "${!found_json[@]}"; do
        # Extract key info from metadata
        json_file="${found_json[$i]}"
        mod_time=$(stat -c "%y" "$json_file" 2>/dev/null | cut -d' ' -f1,2 | cut -d'.' -f1)
        pathogen_info=$(parse_metadata_json "$json_file" "pathogens" | tr ',' ' ' | wc -w)
        state_info=$(parse_metadata_json "$json_file" "states" | tr ',' ' ' | wc -w)
        year_info=$(parse_metadata_json "$json_file" "years")
        
        printf "%2d) %s\n" $((i+1)) "$json_file"
        printf "    Modified: %s | %s pathogens | %s states | Years: %s\n" \
               "$mod_time" "$pathogen_info" "$state_info" "$year_info"
        echo ""
    done
    echo " 0) Run new preprocessing (start from raw data)"
    echo ""
    read -p "Select dataset to use [1]: " metadata_choice
    metadata_choice=${metadata_choice:-1}
    
    if [[ "$metadata_choice" == "0" ]]; then
        input_method=2  # Manual preprocessing
    elif [[ "$metadata_choice" -ge 1 && "$metadata_choice" -le ${#found_json[@]} ]]; then
        input_method=1  # Use metadata
        selected_metadata="${found_json[$((metadata_choice-1))]}"
    else
        echo "Invalid selection. Starting manual preprocessing."
        input_method=2
    fi
else
    echo "No preprocessed data found."
    echo ""
    echo "Select data input method:"
    echo "1) Manual file selection (raw data → preprocessing → analysis)"
    echo "2) Load saved configuration"
    read -p "Enter selection [1]: " input_method
    
    # Remap choices since we don't have metadata option
    if [[ "$input_method" == "2" ]]; then
        input_method=4  # Configuration
    else
        input_method=2  # Manual
    fi
fi

# Map input method to workflow mode
if [[ "$input_method" == "1" || "$input_method" == "3" ]]; then
    workflow_mode=2  # Use existing preprocessed data
elif [[ "$input_method" == "2" ]]; then
    workflow_mode=1  # Full preprocessing
elif [[ "$input_method" == "4" ]]; then
    # Load saved configuration
    echo ""
    echo "Available configurations:"
    if [[ -d "./configs" ]]; then
        ls -1 ./configs/*.sh 2>/dev/null | nl -w2 -s') '
        echo ""
        read -p "Select configuration file number: " config_choice
        config_file=$(ls -1 ./configs/*.sh 2>/dev/null | sed -n "${config_choice}p")
        if [[ -f "$config_file" ]]; then
            echo "Loading configuration from: $config_file"
            source "$config_file"
            # Configuration loaded, skip to parameter validation
            workflow_mode=2
        else
            echo "Invalid selection. Starting manual setup."
            input_method=2
            workflow_mode=1
        fi
    else
        echo "No saved configurations found. Starting manual setup."
        input_method=2
        workflow_mode=1
    fi
fi

# Initialize variables
preprocessed_data=""
preprocessed_metadata=""
mmwrFile=""
censusFileB=""
censusFileP=""

# Handle metadata-first approach
if [[ "$input_method" == "1" ]]; then
    # Use the selected metadata file
    preprocessed_metadata="$selected_metadata"
    
    echo ""
    echo "Loading data from metadata: $preprocessed_metadata"
    echo ""
    
    # Extract all file paths from metadata
    metadata_dir=$(dirname "$preprocessed_metadata")
    
    # Load MMWR data file
    mmwr_filename=$(parse_metadata_json "$preprocessed_metadata" "output_file")
    if [[ -n "$mmwr_filename" ]]; then
        mmwrFile="${metadata_dir}/${mmwr_filename}"
        if [[ ! -f "$mmwrFile" ]]; then
            # Try without directory prefix
            mmwrFile="$mmwr_filename"
        fi
        if [[ -f "$mmwrFile" ]]; then
            echo "✓ MMWR data file: $mmwrFile"
            preprocessed_data="$mmwrFile"
        else
            echo "✗ MMWR data file not found: $mmwr_filename"
            exit 1
        fi
    fi
    
    # Load original source file paths from metadata
    original_mmwr=$(parse_metadata_json "$preprocessed_metadata" "source_file")
    if [[ -n "$original_mmwr" ]]; then
        echo "✓ Original MMWR source: $original_mmwr"
    fi
    
    # Load original census files (these will be used by the pipeline)
    original_census_b=$(parse_metadata_json "$preprocessed_metadata" "census_file_bacterial_original")
    if [[ -n "$original_census_b" ]]; then
        censusFileB="$original_census_b"
        if [[ -f "$censusFileB" ]]; then
            echo "✓ Original bacterial census: $censusFileB"
        else
            echo "⚠️  Warning: Original bacterial census not found at: $censusFileB"
            # Fall back to preprocessed census if original not found
            census_b_file=$(parse_metadata_json "$preprocessed_metadata" "census_file_bacterial_preprocessed")
            if [[ -n "$census_b_file" ]]; then
                censusFileB="${metadata_dir}/${census_b_file}"
                if [[ -f "$censusFileB" ]]; then
                    echo "   Using preprocessed census instead: $censusFileB"
                else
                    echo "✗ Neither original nor preprocessed bacterial census found"
                    exit 1
                fi
            fi
        fi
    else
        # For older metadata files without original paths, use preprocessed
        census_b_file=$(parse_metadata_json "$preprocessed_metadata" "census_file_bacterial_preprocessed")
        if [[ -n "$census_b_file" ]]; then
            censusFileB="${metadata_dir}/${census_b_file}"
            if [[ -f "$censusFileB" ]]; then
                echo "✓ Bacterial census (preprocessed): $censusFileB"
            else
                echo "✗ Bacterial census not found: $census_b_file"
                exit 1
            fi
        fi
    fi
    
    original_census_p=$(parse_metadata_json "$preprocessed_metadata" "census_file_parasitic_original")
    if [[ -n "$original_census_p" ]]; then
        censusFileP="$original_census_p"
        if [[ -f "$censusFileP" ]]; then
            echo "✓ Original parasitic census: $censusFileP"
        else
            echo "⚠️  Warning: Original parasitic census not found at: $censusFileP"
            # Fall back to preprocessed census if original not found
            census_p_file=$(parse_metadata_json "$preprocessed_metadata" "census_file_parasitic_preprocessed")
            if [[ -n "$census_p_file" ]]; then
                censusFileP="${metadata_dir}/${census_p_file}"
                if [[ -f "$censusFileP" ]]; then
                    echo "   Using preprocessed census instead: $censusFileP"
                else
                    echo "✗ Neither original nor preprocessed parasitic census found"
                    exit 1
                fi
            fi
        fi
    else
        # For older metadata files without original paths, use preprocessed
        census_p_file=$(parse_metadata_json "$preprocessed_metadata" "census_file_parasitic_preprocessed")
        if [[ -n "$census_p_file" ]]; then
            censusFileP="${metadata_dir}/${census_p_file}"
            if [[ -f "$censusFileP" ]]; then
                echo "✓ Parasitic census (preprocessed): $censusFileP"
            else
                echo "✗ Parasitic census not found: $census_p_file"
                exit 1
            fi
        fi
    fi
    
    # Skip directly to parameter collection
    workflow_mode=2
fi

# Mode 1: Preprocessing
if [[ "$workflow_mode" == "1" ]]; then
    echo ""
    echo "======== Input Files ========"
    
    # Default data file
    defaultMmwrFile="${DEFAULT_DATA_DIR}/mmwr9624_May2025.sas7bdat"
    
    read -p "MMWR data file [${defaultMmwrFile}]: " mmwrFile
    mmwrFile=${mmwrFile:-$defaultMmwrFile}
    
    # Validate file exists
    if [ ! -f "${mmwrFile}" ]; then
        echo "Error: MMWR file does not exist: ${mmwrFile}"
        echo "$(date): Missing MMWR file: ${mmwrFile}" >> "$error_log"
        echo "Exiting."
        exit 1
    fi
    
    # Census files for preprocessing - REQUIRED
    # Default census data files
    defaultCensusFileB="${DEFAULT_DATA_DIR}/cen9624.sas7bdat"
    defaultCensusFileP="${DEFAULT_DATA_DIR}/cen9624_para.sas7bdat"
    
    echo ""
    echo "IMPORTANT: Census files are REQUIRED for this pipeline."
    echo "These files contain population data needed to calculate incidence rates."
    
    read -p "Census file (bacterial) [${defaultCensusFileB}]: " censusFileB
    censusFileB=${censusFileB:-$defaultCensusFileB}
    
    read -p "Census file (parasitic) [${defaultCensusFileP}]: " censusFileP
    censusFileP=${censusFileP:-$defaultCensusFileP}
    
    # Check census files - required for analysis
    if [ ! -f "${censusFileB}" ]; then
        echo "ERROR: Census bacterial file does not exist: ${censusFileB}"
        echo "$(date): Missing required census bacterial file: ${censusFileB}" >> "$error_log"
        echo "Census files are required for accurate rate calculations."
        exit 1
    fi
    
    if [ ! -f "${censusFileP}" ]; then
        echo "ERROR: Census parasitic file does not exist: ${censusFileP}"
        echo "$(date): Missing required census parasitic file: ${censusFileP}" >> "$error_log"
        echo "Census files are required for accurate rate calculations."
        exit 1
    fi
    
    # Set output location
    echo ""
    echo "======== Output Settings ========"
    
    # Organize preprocessed data in a clean directory structure
    timestamp=$(date +%Y%m%d_%H%M%S)
    defaultPreprocessedDir="preprocessed/${timestamp}"
    read -p "Output directory for preprocessed data [${defaultPreprocessedDir}]: " preprocessedDir
    preprocessedDir=${preprocessedDir:-$defaultPreprocessedDir}
    
    # Default output base name (derived from input filename)
    fileBasename=$(basename "${mmwrFile}" | sed 's/\.[^.]*$//')
    defaultOutputBase="foodnet_data_${fileBasename}"
    read -p "Base name for output files [${defaultOutputBase}]: " outputBase
    outputBase=${outputBase:-$defaultOutputBase}
    
    # Ask about metadata generation
    echo ""
    read -p "Generate metadata JSON? (y/n) [y]: " generate_metadata
    generate_metadata=${generate_metadata:-y}
    if [[ "$generate_metadata" =~ ^[Yy]$ ]]; then
        metadata_param="--generateMetadata true"
    else
        metadata_param="--generateMetadata false"
    fi
    
    echo "Running preprocessing..."
    
    # Ask about preprocessing resource allocation
    echo ""
    echo "======== Preprocessing Resource Allocation ========"
    echo "Select resource level for preprocessing:"
    echo "1) Standard resources (4 cores, 8GB memory)"
    echo "2) High-performance resources (16 cores, 32GB memory)"
    echo "3) Maximum resources (32 cores, 64GB memory)"
    echo "4) Auto-select based on data size"
    read -p "Enter selection [4]: " preproc_resource
    preproc_resource=${preproc_resource:-4}
    
    # If auto-select, check file size to determine resource allocation
    if [[ "$preproc_resource" == "4" ]]; then
        mmwr_size=$(stat -c%s "${mmwrFile}" 2>/dev/null || stat -f%z "${mmwrFile}" 2>/dev/null || echo "0")
        # Convert to MB
        mmwr_size_mb=$((mmwr_size / 1024 / 1024))
        echo "Detected MMWR file size: ${mmwr_size_mb}MB"
        
        if [[ $mmwr_size_mb -lt 100 ]]; then
            # Small file (<100MB)
            preproc_resource=1
            echo "Auto-selected: Standard resources for small dataset"
        elif [[ $mmwr_size_mb -lt 500 ]]; then
            # Medium file (100MB-500MB)
            preproc_resource=2
            echo "Auto-selected: High-performance resources for medium dataset"
        else
            # Large file (>500MB)
            preproc_resource=3
            echo "Auto-selected: Maximum resources for large dataset"
        fi
    fi
    
    # Ask if preprocessing should resume a previous run
    echo ""
    echo "Resume previous preprocessing if it exists?"
    echo "1) No, start preprocessing from scratch"
    echo "2) Yes, resume from last successful step"
    read -p "Enter selection [1]: " preproc_resume
    preproc_resume=${preproc_resume:-1}
    
    preproc_resume_flag=""
    if [[ "$preproc_resume" == "2" ]]; then
        preproc_resume_flag="-resume"
        echo "Preprocessing will resume from last successful step if possible"
    fi
    
    # Set preprocessing resources based on selection
    case $preproc_resource in
        1)  # Standard resources
            preproc_cores=4
            preproc_memory="8.GB"
            echo "Using standard resources for preprocessing"
            ;;
        2)  # High-performance resources
            preproc_cores=16
            preproc_memory="32.GB"
            echo "Using high-performance resources for preprocessing"
            ;;
        3)  # Maximum resources
            preproc_cores=32
            preproc_memory="64.GB"
            echo "Using maximum resources for preprocessing"
            ;;
        *)  # Default to high-performance
            preproc_cores=16
            preproc_memory="32.GB"
            echo "Using high-performance resources for preprocessing"
            ;;
    esac
    
    # Run the preprocessing workflow with HPC optimization
    preprocess_cmd="nextflow run main.nf -profile singularity,production -entry PREPROCESS_ONLY ${preproc_resume_flag} \
      --mmwrFile \"${mmwrFile}\" \
      --censusFileB \"${censusFileB}\" \
      --censusFileP \"${censusFileP}\" \
      --outdir \"${preprocessedDir}\" \
      --outputBase \"${outputBase}\" \
      --cores ${preproc_cores} \
      -process.memory ${preproc_memory} \
      ${metadata_param}"
    
    echo "$(date): Running preprocessing command: ${preprocess_cmd}" >> "$error_log"
    
    if ! eval $preprocess_cmd; then
        echo "Preprocessing failed."
        echo "Check .nextflow.log for details."
        echo "$(date): Preprocessing failed, check .nextflow.log" >> "$error_log"
        exit 1
    fi
    
    # Set paths to preprocessed data and metadata in the new directory structure
    # Files are published to preprocessed/ subdirectory by the Nextflow module
    preprocessed_data="${preprocessedDir}/preprocessed/${outputBase}.csv"
    preprocessed_metadata="${preprocessedDir}/preprocessed/${outputBase}_metadata.json"
    
    # Check for metadata in alternate location (for backward compatibility)
    if [ ! -f "${preprocessed_metadata}" ] && [ -f "${preprocessedDir}/preprocessed/metadata/${outputBase}_metadata.json" ]; then
        echo "Found metadata file in alternate location, using it instead."
        preprocessed_metadata="${preprocessedDir}/preprocessed/metadata/${outputBase}_metadata.json"
    fi
    
    # Check if files were created
    if [ ! -f "${preprocessed_data}" ]; then
        echo "Preprocessing failed to create expected CSV file: ${preprocessed_data}"
        echo "$(date): Missing expected output CSV: ${preprocessed_data}" >> "$error_log"
        echo "Check logs for details."
        exit 1
    fi
    
    echo "Preprocessing complete!"
    echo "- Cleaned data: ${preprocessed_data}"
    
    # Continue to analysis?
    echo ""
    read -p "Continue to analysis? (y/n) [y]: " continue_analysis
    continue_analysis=${continue_analysis:-y}
    if [[ ! "$continue_analysis" =~ ^[Yy]$ ]]; then
        echo "Preprocessing complete. Exiting."
        exit 0
    fi
    
    # Set the MMWR file for analysis to the preprocessed data
    mmwrFile=$preprocessed_data
    
    # Update census files to use preprocessed (aggregated) versions created during preprocessing
    # This ensures analysis uses state-level data instead of county-level raw data
    preprocessed_dir=$(dirname "${preprocessed_data}")
    preprocessed_basename=$(basename "${preprocessed_data}" .csv)
    
    potential_census_b="${preprocessed_dir}/${preprocessed_basename}_census_bacterial.csv"
    potential_census_p="${preprocessed_dir}/${preprocessed_basename}_census_parasitic.csv"
    
    if [[ -f "$potential_census_b" ]]; then
        echo "Switching to preprocessed bacterial census file: $potential_census_b"
        censusFileB="$potential_census_b"
    fi
    
    if [[ -f "$potential_census_p" ]]; then
        echo "Switching to preprocessed parasitic census file: $potential_census_p"
        censusFileP="$potential_census_p"
    fi
    
    # CRITICAL CHECK: Ensure cleaned census files are being used for analysis
    # This prevents the pipeline from using raw county-level data which causes join explosion
    if [[ ! -f "$potential_census_b" || ! -f "$potential_census_p" ]]; then
        echo ""
        echo "ERROR: Preprocessed census files not found!"
        echo "Expected files:"
        echo "  - Bacterial: $potential_census_b"
        echo "  - Parasitic: $potential_census_p"
        echo ""
        echo "The preprocessing step should have created state-level aggregated census files."
        echo "Using raw county-level census data will cause analysis failures."
        echo ""
        echo "Possible solutions:"
        echo "1. Check if preprocessing completed successfully"
        echo "2. Verify census files were created in the preprocessing output"
        echo "3. Re-run preprocessing if census files are missing"
        echo ""
        echo "Cannot proceed with analysis using raw census data."
        exit 1
    fi
    
    echo "✓ Using preprocessed (state-level) census files for analysis"
    
# Mode 2: Use existing preprocessed data
elif [[ "$workflow_mode" == "2" ]]; then
    echo ""
    echo "======== Input Files ========"

    # Search for preprocessed CSV files in the organized directory structure
    echo "Searching for preprocessed CSV files..."
    mapfile -t found_csv < <(find . -type f \( \
        -path "./preprocessed/*/foodnet_data_*.csv" -o \
        -path "preprocessed_*/preprocessed/*.csv" -o \
        -path "preprocessed_*/*.csv" -o \
        -name "*mmwr*.csv" \
        \) 2>/dev/null)
    if [[ ${#found_csv[@]} -gt 0 ]]; then
        echo "Found the following preprocessed CSV files:"
        for i in "${!found_csv[@]}"; do
            printf "%2d) %s\n" $((i+1)) "${found_csv[$i]}"
        done
        echo "$(( ${#found_csv[@]} + 1 ))) Enter a file path manually"
        read -p "Select a CSV file [1]: " csv_choice
        csv_choice=${csv_choice:-1}
        if [[ "$csv_choice" -ge 1 && "$csv_choice" -le ${#found_csv[@]} ]]; then
            preprocessed_data="${found_csv[$((csv_choice-1))]}"
        else
            read -p "Preprocessed data file: " preprocessed_data
        fi
    else
        echo ""
        echo "No preprocessed CSV files found in expected locations."
        echo "Searched in:"
        echo "  - ./preprocessed/*/foodnet_data_*.csv"
        echo "  - ./preprocessed_*/preprocessed/*.csv (legacy)"
        echo "  - ./preprocessed_*/*.csv (legacy)"
        echo ""
        echo "Options:"
        echo "1) Go back and run full preprocessing workflow"
        echo "2) Enter preprocessed data file path manually"
        read -p "Select option [1]: " no_data_choice
        no_data_choice=${no_data_choice:-1}
        
        if [[ "$no_data_choice" == "1" ]]; then
            echo "Returning to workflow selection..."
            echo ""
            # Reset workflow mode to trigger main menu
            workflow_mode=""
            exec "$0" "$@"
        else
            read -p "Preprocessed data file: " preprocessed_data
        fi
    fi

    # Validate file exists
    if [ ! -f "${preprocessed_data}" ]; then
        echo "Error: Preprocessed data file does not exist: ${preprocessed_data}"
        echo "$(date): Missing preprocessed data file: ${preprocessed_data}" >> "$error_log"
        echo "Exiting."
        exit 1
    fi

    # Set MMWR file to preprocessed data
    mmwrFile=$preprocessed_data
    
    # Try to discover census files in the same directory as preprocessed data
    preprocessed_dir=$(dirname "${preprocessed_data}")
    preprocessed_basename=$(basename "${preprocessed_data}" .csv)
    
    # Look for preprocessed census files with expected naming pattern
    potential_census_b="${preprocessed_dir}/${preprocessed_basename}_census_bacterial.csv"
    potential_census_p="${preprocessed_dir}/${preprocessed_basename}_census_parasitic.csv"
    
    if [[ -f "$potential_census_b" ]]; then
        echo "Found preprocessed bacterial census file: $potential_census_b"
        censusFileB="$potential_census_b"
    fi
    
    if [[ -f "$potential_census_p" ]]; then
        echo "Found preprocessed parasitic census file: $potential_census_p"
        censusFileP="$potential_census_p"
    fi
    
    # CRITICAL CHECK: Ensure preprocessed census files are available for analysis
    # This prevents the pipeline from using raw county-level data which causes join explosion
    if [[ ! -f "$potential_census_b" || ! -f "$potential_census_p" ]]; then
        echo ""
        echo "ERROR: Preprocessed census files not found!"
        echo "Expected files:"
        echo "  - Bacterial: $potential_census_b"
        echo "  - Parasitic: $potential_census_p"
        echo ""
        echo "When using existing preprocessed data, the census files should have been"
        echo "preprocessed to state-level aggregation in the same directory."
        echo ""
        echo "Possible solutions:"
        echo "1. Use Option 1 (Full preprocessing workflow) instead"
        echo "2. Manually preprocess census files to state level"
        echo "3. Check if files exist with different naming patterns"
        echo ""
        echo "Cannot proceed with analysis using potentially raw census data."
        exit 1
    fi
    
    echo "✓ Using preprocessed (state-level) census files for analysis"

    # Search for preprocessed JSON metadata files
    echo ""
    echo "Searching for preprocessed metadata JSON files..."
    mapfile -t found_json < <(find . -type f -path "*/preprocessed/*.json" -o -path "preprocessed_*/preprocessed/*.json" 2>/dev/null)
    if [[ ${#found_json[@]} -gt 0 ]]; then
        echo "Found the following preprocessed metadata JSON files:"
        for i in "${!found_json[@]}"; do
            printf "%2d) %s\n" $((i+1)) "${found_json[$i]}"
        done
        echo "$(( ${#found_json[@]} + 1 ))) Enter a file path manually or leave blank for none"
        read -p "Select a metadata file [${#found_json[@]}+1]: " json_choice
        json_choice=${json_choice:-$(( ${#found_json[@]} + 1 ))}
        if [[ "$json_choice" -ge 1 && "$json_choice" -le ${#found_json[@]} ]]; then
            preprocessed_metadata="${found_json[$((json_choice-1))]}"
        else
            read -p "Enter path to metadata file (leave blank if none): " preprocessed_metadata
        fi
    else
        echo "No preprocessed metadata JSON files found. Enter path manually or leave blank for none."
        read -p "Enter path to metadata file (leave blank if none): " preprocessed_metadata
    fi

    # Validate metadata file exists if provided
    if [ -n "${preprocessed_metadata}" ] && [ ! -f "${preprocessed_metadata}" ]; then
        # Try to find metadata in alternate location (in metadata subdirectory)
        metadata_dir=$(dirname "${preprocessed_metadata}")
        metadata_basename=$(basename "${preprocessed_metadata}")
        alt_metadata_path="${metadata_dir}/metadata/${metadata_basename}"

        if [ -f "${alt_metadata_path}" ]; then
            echo "Found metadata file in alternate location: ${alt_metadata_path}"
            preprocessed_metadata="${alt_metadata_path}"
        else
            echo "Warning: Metadata file does not exist: ${preprocessed_metadata}"
            echo "$(date): Missing metadata file: ${preprocessed_metadata}" >> "$error_log"
            read -p "Continue anyway? (y/n) [n]: " continue_choice
            continue_choice=${continue_choice:-n}
            if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
                echo "Exiting."
                exit 1
            fi
            # Clear metadata file if continuing without it
            preprocessed_metadata=""
        fi
    fi
fi

# Initialize default pathogen and state lists
ALL_PATHOGENS="SALMONELLA,CAMPYLOBACTER,SHIGA,STEC,CYCLOSPORA,LISTERIA,VIBRIO,YERSINIA"
ALL_STATES=""

# Load metadata if available
has_metadata=false
available_serotypes=""
available_serogroups=""

if [[ -n "${preprocessed_metadata}" && -f "${preprocessed_metadata}" ]]; then
    has_metadata=true
    
    # If not already loaded from metadata-first approach, load census files
    if [[ "$input_method" != "1" ]]; then
        echo ""
        echo "======== Loading from Metadata ========"
        preprocessed_dir=$(dirname "${preprocessed_data}")
        
        # Load census files from metadata
        census_bacterial_filename=$(parse_metadata_json "${preprocessed_metadata}" "census_file_bacterial_preprocessed")
        if [[ -n "$census_bacterial_filename" ]]; then
            census_bacterial_path="${preprocessed_dir}/${census_bacterial_filename}"
            if [[ -f "$census_bacterial_path" ]]; then
                censusFileB="$census_bacterial_path"
                echo "✓ Bacterial census loaded from metadata"
            fi
        fi
        
        census_parasitic_filename=$(parse_metadata_json "${preprocessed_metadata}" "census_file_parasitic_preprocessed")
        if [[ -n "$census_parasitic_filename" ]]; then
            census_parasitic_path="${preprocessed_dir}/${census_parasitic_filename}"
            if [[ -f "$census_parasitic_path" ]]; then
                censusFileP="$census_parasitic_path"
                echo "✓ Parasitic census loaded from metadata"
            fi
        fi
    fi
    
    # Extract available data from metadata
    metadata_pathogens=$(parse_metadata_json "${preprocessed_metadata}" "pathogens")
    if [ -n "$metadata_pathogens" ]; then
        ALL_PATHOGENS="$metadata_pathogens"
    else
        ALL_PATHOGENS="SALMONELLA,CAMPYLOBACTER,SHIGA,STEC,CYCLOSPORA,LISTERIA,VIBRIO,YERSINIA"
    fi
    
    metadata_states=$(parse_metadata_json "${preprocessed_metadata}" "states")
    if [ -n "$metadata_states" ]; then
        ALL_STATES="$metadata_states"
    fi
    
    # Load serotypes/serogroups from metadata
    # echo "DEBUG: Parsing salmonella_serotypes from: ${preprocessed_metadata}" >&2
    available_serotypes=$(parse_metadata_json "${preprocessed_metadata}" "salmonella_serotypes")
    # echo "DEBUG: Raw salmonella result: '$available_serotypes' (length: ${#available_serotypes})" >&2
    
    # echo "DEBUG: Parsing stec_serogroups" >&2
    available_serogroups=$(parse_metadata_json "${preprocessed_metadata}" "stec_serogroups")
    # echo "DEBUG: Raw stec result: '$available_serogroups' (length: ${#available_serogroups})" >&2
    
    if [[ -n "$available_serotypes" ]]; then
        echo "✓ Salmonella serotypes available in metadata: $available_serotypes"
    else
        echo "⚠️  No Salmonella serotypes found in metadata"
    fi
    if [[ -n "$available_serogroups" ]]; then
        echo "✓ STEC serogroups available in metadata: $available_serogroups"
    else
        echo "⚠️  No STEC serogroups found in metadata"
    fi
fi

# Ensure key pathogens are always available (for backward compatibility)
if [[ ! "$ALL_PATHOGENS" == *"CAMPYLOBACTER"* ]]; then
    ALL_PATHOGENS="$ALL_PATHOGENS,CAMPYLOBACTER"
fi
if [[ ! "$ALL_PATHOGENS" == *"CYCLOSPORA"* ]]; then
    ALL_PATHOGENS="$ALL_PATHOGENS,CYCLOSPORA"
fi

# Clean up any leading commas from pathogen list
ALL_PATHOGENS=$(echo "$ALL_PATHOGENS" | sed 's/^,//')

# Validate all files after census file discovery is complete
echo ""
echo "======== File Validation ========"
if ! validate_required_files; then
    echo "Cannot proceed with missing or invalid files. Please check the file paths and try again."
    exit 1
fi

# Get travel status filter
echo ""
echo "======== Travel Status Filter ========"
echo "Select travel status to include:"
echo "1) All travel statuses (NO, UNKNOWN, YES)"
echo "2) Only non-travel related (NO only)"
echo "3) Non-travel and unknown (NO, UNKNOWN)"
echo "4) Travel-related only (YES only)"
echo "5) Custom selection"
read -p "Enter selection [1]: " travel_choice
travel_choice=${travel_choice:-1}

case $travel_choice in
    1) travel="NO,UNKNOWN,YES" ;;
    2) travel="NO" ;;
    3) travel="NO,UNKNOWN" ;;
    4) travel="YES" ;;
    5)
        echo "Enter comma-separated travel statuses (NO,UNKNOWN,YES):"
        read -p "Travel statuses: " travel
        # Default if empty
        travel=${travel:-NO,UNKNOWN,YES}
        ;;
    *)
        echo "Invalid selection. Using default (All travel statuses)."
        travel="NO,UNKNOWN,YES"
        ;;
esac

# Get CIDT/culture method filter
echo ""
echo "======== CIDT/Culture Method Filter ========"
echo "Select CIDT/culture methods to include:"
echo "1) All methods (CIDT+, CX+, PARASITIC)"
echo "2) Culture positive only (CX+)"
echo "3) CIDT positive only (CIDT+)"
echo "4) Custom selection"
read -p "Enter selection [1]: " cidt_choice
cidt_choice=${cidt_choice:-1}

case $cidt_choice in
    1) cidt="CIDT+,CX+,PARASITIC" ;;
    2) cidt="CX+" ;;
    3) cidt="CIDT+" ;;
    4)
        echo "Enter comma-separated CIDT/culture methods (CIDT+,CX+,PARASITIC):"
        read -p "CIDT/culture methods: " cidt
        # Default if empty
        cidt=${cidt:-CIDT+,CX+,PARASITIC}
        ;;
    *)
        echo "Invalid selection. Using default (All methods)."
        cidt="CIDT+,CX+,PARASITIC"
        ;;
esac

# Get pathogen selection
echo ""
echo "======== Pathogen Selection ========"

# Discover available pathogens from data
echo "Discovering available pathogens from your data..."
available_pathogens=$(discover_pathogens "$mmwrFile" "$preprocessed_metadata")
complete_pipeline_stage "Pathogen discovery"

# Validate discovered pathogens - check for garbage data
valid_pathogens=""
if [[ -n "$available_pathogens" ]]; then
    IFS=',' read -ra PATHOGEN_ARRAY <<< "$available_pathogens"
    for p in "${PATHOGEN_ARRAY[@]}"; do
        if [[ -n "$p" ]]; then
            # Filter out obvious garbage - check length and exclude comment-like text
            clean_p=$(echo "$p" | tr -d '"' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            
            # Simple validation: reasonable length and not obviously a comment
            if [[ ${#clean_p} -ge 3 ]] && [[ ${#clean_p} -le 30 ]] && \
               [[ ! "$clean_p" =~ TRAVEL ]] && [[ ! "$clean_p" =~ FAMILY ]] && \
               [[ ! "$clean_p" =~ PFGE ]] && [[ ! "$clean_p" =~ DAYCARE ]] && \
               [[ ! "$clean_p" =~ ONSET ]] && [[ ! "$clean_p" =~ MOM ]] && \
               [[ ! "$clean_p" =~ FATHER ]] && [[ ! "$clean_p" =~ DIALYSIS ]] && \
               [[ ! "$clean_p" =~ "=" ]] && [[ ! "$clean_p" =~ "\"" ]]; then
                if [[ -n "$valid_pathogens" ]]; then
                    valid_pathogens="$valid_pathogens,$clean_p"
                else
                    valid_pathogens="$clean_p"
                fi
            fi
        fi
    done
fi

if [[ -n "$valid_pathogens" ]]; then
    echo "Available pathogens in this dataset:"
    IFS=',' read -ra PATHOGEN_ARRAY <<< "$valid_pathogens"
    for p in "${PATHOGEN_ARRAY[@]}"; do
        if [[ -n "$p" ]]; then
            echo "- $p"
        fi
    done
    available_pathogens="$valid_pathogens"
    echo ""
    
    echo "1) Run ALL available pathogens"
    echo "2) Select specific pathogens"
    read -p "Enter selection [2]: " pathogen_mode
    pathogen_mode=${pathogen_mode:-2}
    
    if [[ "$pathogen_mode" == "1" ]]; then
        pathogens="$available_pathogens"
        echo "Selected: ALL available pathogens"
    else
        echo ""
        echo "Enter pathogens to analyze (comma-separated from list above):"
        read -p "Pathogens: " pathogens
        
        # Validate selection against available pathogens
        if [[ -z "$pathogens" ]]; then
            # Default to first two pathogens if available
            IFS=',' read -ra DEFAULT_ARRAY <<< "$available_pathogens"
            if [[ ${#DEFAULT_ARRAY[@]} -ge 2 ]]; then
                pathogens="${DEFAULT_ARRAY[0]},${DEFAULT_ARRAY[1]}"
                echo "Using default: $pathogens"
            else
                pathogens="$available_pathogens"
                echo "Using all available: $pathogens"
            fi
        fi
    fi
else
    echo "⚠️  Could not automatically discover valid pathogens from data."
    echo ""
    echo "This may be due to:"
    echo "- Unexpected data format or column structure"
    echo "- Missing pathogen column in the data file"
    echo "- Data quality issues"
    echo ""
    echo "Common FoodNet pathogens include:"
    echo "- SALMONELLA"
    echo "- CAMPYLOBACTER" 
    echo "- SHIGA (or STEC)"
    echo "- CYCLOSPORA"
    echo "- LISTERIA"
    echo "- VIBRIO"
    echo "- YERSINIA"
    echo ""
    echo "Please enter pathogens manually from your data:"
    read -p "Pathogens (comma-separated): " pathogens
    
    if [[ -z "$pathogens" ]]; then
        echo "Error: No pathogens specified. Cannot proceed with analysis."
        exit 1
    fi
    
    # Clean up manually entered pathogens
    pathogens=$(echo "$pathogens" | tr '[:lower:]' '[:upper:]' | sed 's/[[:space:]]//g')
fi

# Advanced pathogen configuration - serotypes and serogroups
echo ""
echo "======== Advanced Pathogen Configuration ========"

# STEC Serogroup Selection (only if STEC is selected)
if [[ "$pathogens" == *"STEC"* ]]; then
    echo ""
    echo "--- STEC Serogroup Selection ---"
    
    # Discover available STEC serogroups from data
    echo "Discovering STEC serogroups in your data..."
    available_stec_serogroups=$(discover_pathogen_subtypes "STEC" "$mmwrFile" "$preprocessed_metadata" "serogroup")
    
    if [[ -n "$available_stec_serogroups" ]]; then
        echo "Available STEC serogroups in dataset: $available_stec_serogroups"
    fi
    
    echo "STEC analysis uses serogroups (epidemiological standard):"
    echo "1) ALL"
    echo "2) O157"
    echo "3) NON-O157"
    read -p "Select STEC serogroups [1]: " stec_serogroup_choice
    stec_serogroup_choice=${stec_serogroup_choice:-1}
    
    case $stec_serogroup_choice in
        1) stec_serogroups="ALL" ;;
        2) stec_serogroups="O157" ;;
        3) stec_serogroups="NON-O157" ;;
    esac
    echo "Selected STEC serogroups: $stec_serogroups"
else
    stec_serogroups=""
fi

# Salmonella Serotype Selection (only if SALMONELLA is selected)
if [[ "$pathogens" == *"SALMONELLA"* ]]; then
    echo ""
    echo "--- Salmonella Serotype Selection ---"
    
    # Use metadata serotypes if available, otherwise discover from data
    if [[ -n "$available_serotypes" ]]; then
        # Parse serotypes and filter out the "...and X others" entry if present
        IFS='|' read -ra RAW_SEROTYPES_ARRAY <<< "$available_serotypes"
        SALMONELLA_SEROTYPES_ARRAY=()
        total_count=""
        
        for serotype in "${RAW_SEROTYPES_ARRAY[@]}"; do
            if [[ "$serotype" =~ ^\.\.\. ]]; then
                # Extract total count from "...and X others"
                total_count=$(echo "$serotype" | grep -o '[0-9]\+' | head -1)
            else
                SALMONELLA_SEROTYPES_ARRAY+=("$serotype")
            fi
        done
        
        total_serotypes=${#SALMONELLA_SEROTYPES_ARRAY[@]}
        if [[ -n "$total_count" ]]; then
            total_serotypes=$((total_serotypes + total_count))
        fi
        
        echo ""
        echo "Salmonella serotypes in your data ($total_serotypes total):"
        echo ""
        
        # Show hierarchical selection
        echo "1) Top 5 most common serotypes:"
        if [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 5 ]]; then
            echo "   ${SALMONELLA_SEROTYPES_ARRAY[0]} | ${SALMONELLA_SEROTYPES_ARRAY[1]} | ${SALMONELLA_SEROTYPES_ARRAY[2]} | ${SALMONELLA_SEROTYPES_ARRAY[3]} | ${SALMONELLA_SEROTYPES_ARRAY[4]}"
        else
            echo "   $(IFS=' | '; echo "${SALMONELLA_SEROTYPES_ARRAY[*]}")"
        fi
        echo ""
        
        if [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 10 ]]; then
            echo "2) Top 10 most common serotypes:"
            echo "   + ${SALMONELLA_SEROTYPES_ARRAY[5]}, ${SALMONELLA_SEROTYPES_ARRAY[6]}, ${SALMONELLA_SEROTYPES_ARRAY[7]}, ${SALMONELLA_SEROTYPES_ARRAY[8]}, ${SALMONELLA_SEROTYPES_ARRAY[9]}"
            echo ""
        fi
        
        if [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 15 ]]; then
            echo "3) Top 15 most common serotypes:"
            echo "   + ${SALMONELLA_SEROTYPES_ARRAY[10]}, ${SALMONELLA_SEROTYPES_ARRAY[11]}, ${SALMONELLA_SEROTYPES_ARRAY[12]}, ${SALMONELLA_SEROTYPES_ARRAY[13]}, ${SALMONELLA_SEROTYPES_ARRAY[14]}"
            echo ""
        fi
        
        if [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 20 ]]; then
            echo "4) Top 20 most common serotypes:"
            echo "   + ${SALMONELLA_SEROTYPES_ARRAY[15]}, ${SALMONELLA_SEROTYPES_ARRAY[16]}, ${SALMONELLA_SEROTYPES_ARRAY[17]}, ${SALMONELLA_SEROTYPES_ARRAY[18]}, ${SALMONELLA_SEROTYPES_ARRAY[19]}"
            echo ""
        fi
        
        # Determine next option numbers
        next_opt=2
        [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 10 ]] && ((next_opt++))
        [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 15 ]] && ((next_opt++))
        [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 20 ]] && ((next_opt++))
        
        echo "$next_opt) All serotypes (no filtering) - $total_serotypes total serotypes"
        echo "$((next_opt + 1))) Enter specific serotypes manually"
        echo ""
        
        read -p "Select Salmonella serotype analysis [$next_opt]: " sal_sero_choice
        sal_sero_choice=${sal_sero_choice:-$next_opt}
        
        case $sal_sero_choice in
            1)
                # Top 5
                if [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 5 ]]; then
                    salmonella_serotypes="${SALMONELLA_SEROTYPES_ARRAY[0]},${SALMONELLA_SEROTYPES_ARRAY[1]},${SALMONELLA_SEROTYPES_ARRAY[2]},${SALMONELLA_SEROTYPES_ARRAY[3]},${SALMONELLA_SEROTYPES_ARRAY[4]}"
                else
                    salmonella_serotypes=$(IFS=','; echo "${SALMONELLA_SEROTYPES_ARRAY[*]}")
                fi
                ;;
            2)
                # Top 10 (if available)
                if [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 10 ]]; then
                    salmonella_serotypes="${SALMONELLA_SEROTYPES_ARRAY[0]},${SALMONELLA_SEROTYPES_ARRAY[1]},${SALMONELLA_SEROTYPES_ARRAY[2]},${SALMONELLA_SEROTYPES_ARRAY[3]},${SALMONELLA_SEROTYPES_ARRAY[4]},${SALMONELLA_SEROTYPES_ARRAY[5]},${SALMONELLA_SEROTYPES_ARRAY[6]},${SALMONELLA_SEROTYPES_ARRAY[7]},${SALMONELLA_SEROTYPES_ARRAY[8]},${SALMONELLA_SEROTYPES_ARRAY[9]}"
                else
                    salmonella_serotypes="ALL"
                fi
                ;;
            3)
                # Top 15 (if available)
                if [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 15 ]]; then
                    salmonella_serotypes=""
                    for i in {0..14}; do
                        [[ $i -gt 0 ]] && salmonella_serotypes+=","
                        salmonella_serotypes+="${SALMONELLA_SEROTYPES_ARRAY[$i]}"
                    done
                else
                    salmonella_serotypes="ALL"
                fi
                ;;
            4)
                # Top 20 (if available)
                if [[ ${#SALMONELLA_SEROTYPES_ARRAY[@]} -ge 20 ]]; then
                    salmonella_serotypes=""
                    for i in {0..19}; do
                        [[ $i -gt 0 ]] && salmonella_serotypes+=","
                        salmonella_serotypes+="${SALMONELLA_SEROTYPES_ARRAY[$i]}"
                    done
                else
                    salmonella_serotypes="ALL"
                fi
                ;;
            $next_opt)
                # All serotypes
                salmonella_serotypes="ALL"
                ;;
            $((next_opt + 1)))
                # Manual entry
                read -p "Enter Salmonella serotypes (comma, pipe, or space separated): " user_serotypes
                if [[ -z "$user_serotypes" ]]; then
                    salmonella_serotypes="ALL"
                else
                    # Process user input with validation
                    salmonella_serotypes=$(process_user_serotypes "$user_serotypes" "$available_serotypes")
                    if [[ -z "$salmonella_serotypes" ]]; then
                        echo "No valid serotypes found. Using ALL serotypes."
                        salmonella_serotypes="ALL"
                    fi
                fi
                ;;
            *)
                # Default to all
                salmonella_serotypes="ALL"
                ;;
        esac
    else
        # Fallback to data discovery if no metadata
        salmonella_serotype_list=$(discover_pathogen_subtypes "SALMONELLA" "$mmwrFile" "$preprocessed_metadata" "serotype")
        if [[ -n "$salmonella_serotype_list" ]]; then
            IFS='|' read -ra SALMONELLA_SEROTYPES_ARRAY <<< "$salmonella_serotype_list"
            echo "Detected Salmonella serotypes in data:"
            for i in "${!SALMONELLA_SEROTYPES_ARRAY[@]}"; do
                printf "%2d) %s\n" $((i+1)) "${SALMONELLA_SEROTYPES_ARRAY[$i]}"
            done
            echo "$(( ${#SALMONELLA_SEROTYPES_ARRAY[@]} + 1 ))) Use all serotypes (no filtering)"
            read -p "Select serotypes or $(( ${#SALMONELLA_SEROTYPES_ARRAY[@]} + 1 )) for all [$(( ${#SALMONELLA_SEROTYPES_ARRAY[@]} + 1 ))]: " sal_sero_choice
            sal_sero_choice=${sal_sero_choice:-$(( ${#SALMONELLA_SEROTYPES_ARRAY[@]} + 1 ))}
            
            if [[ "$sal_sero_choice" -eq $(( ${#SALMONELLA_SEROTYPES_ARRAY[@]} + 1 )) ]]; then
                salmonella_serotypes="ALL"
            else
                # Convert indices to serotype names
                salmonella_serotypes=""
                IFS=',' read -ra IDX <<< "$sal_sero_choice"
                for idx in "${IDX[@]}"; do
                    idx=$((idx-1))
                    if [[ $idx -ge 0 && $idx -lt ${#SALMONELLA_SEROTYPES_ARRAY[@]} ]]; then
                        salmonella_serotypes+="${SALMONELLA_SEROTYPES_ARRAY[$idx]},"
                    fi
                done
                salmonella_serotypes=$(echo "$salmonella_serotypes" | sed 's/,*$//')
            fi
        else
            echo "No serotype information available."
            salmonella_serotypes="ALL"
        fi
    fi
    echo "Selected Salmonella serotypes: $salmonella_serotypes"
else
    salmonella_serotypes=""
fi

# Get state selection
echo ""
echo "======== State Selection ========"
echo "Available states in this dataset:"

# Display each state
IFS=',' read -ra STATE_ARRAY <<< "$ALL_STATES"
for s in "${STATE_ARRAY[@]}"; do
    echo "- $s"
done
echo ""

echo "1) Use ALL available states"
echo "2) Select specific states"
read -p "Enter selection [1]: " state_mode
state_mode=${state_mode:-1}

if [[ "$state_mode" == "1" ]]; then
    # Use all states
    echo "Selected: ALL states"
    states="$ALL_STATES"
else
    # Ask for specific states
    echo ""
    echo "Enter states to analyze (comma-separated with NO spaces)"
    read -p "Leave blank for all states: " states
    states=${states:-"$ALL_STATES"}
fi

# Calculate pathogen count to recommend appropriate resources
IFS=',' read -ra PATHOGEN_COUNT_ARRAY <<< "$pathogens"
pathogen_count=${#PATHOGEN_COUNT_ARRAY[@]}

# NEW: HPC performance profile selection
echo ""
echo "======== HPC Performance Profile ========"
echo "Select a performance profile based on your needs:"
echo "1) Quick Test - Minimal resources for testing (~1 hour, 4GB RAM, 2 cores)"
echo "2) Standard Run - Balanced resources (~4 hours, 16GB RAM, 8 cores)"
echo "3) Production Run - Full resources for accurate results (~8 hours, 32GB RAM, 16 cores)"
echo "4) Maximum Performance - Highest resource allocation (~12 hours, 64GB RAM, 32 cores)"
echo "5) Stable (RECOMMENDED) - Proven Stan settings that work reliably (~3 hours, 16GB RAM, 8 cores)"
echo "6) Custom Resource Configuration"

# Set default based on pathogen count - default to stable (most reliable)
if [[ $pathogen_count -le 1 ]]; then
    # For single pathogen, suggest stable
    default_profile=5
elif [[ $pathogen_count -le 3 ]]; then
    # For 2-3 pathogens, suggest stable
    default_profile=5
else
    # For 4+ pathogens, suggest stable (proven to work)
    default_profile=5
fi

read -p "Enter selection [$default_profile]: " performance_profile
performance_profile=${performance_profile:-$default_profile}

# Set parameters based on performance profile
case $performance_profile in
    1) # Quick Test
        chains=2
        iterations=100
        adapt_delta=0.8
        max_treedepth=8
        cores=2
        memory="4.GB"
        flag="test"
        echo ""
        echo "Quick Test Profile Selected"
        echo "Chains: $chains"
        echo "Iterations: $iterations"
        echo "Adapt delta: $adapt_delta"
        echo "Max treedepth: $max_treedepth"
        echo "Cores: $cores"
        echo "Memory: 4GB"
        ;;
    2) # Standard Run
        chains=4
        iterations=1000
        adapt_delta=0.95
        max_treedepth=10
        cores=8
        memory="16.GB"
        flag="standard"
        echo ""
        echo "Standard Run Profile Selected"
        echo "Chains: $chains"
        echo "Iterations: $iterations"
        echo "Adapt delta: $adapt_delta"
        echo "Max treedepth: $max_treedepth"
        echo "Cores: $cores"
        echo "Memory: 16GB"
        ;;
    3) # Production Run
        chains=8
        iterations=2000
        adapt_delta=0.99
        max_treedepth=12
        cores=16
        memory="32.GB"
        flag="production"
        echo ""
        echo "Production Run Profile Selected"
        echo "Chains: $chains"
        echo "Iterations: $iterations"
        echo "Adapt delta: $adapt_delta"
        echo "Max treedepth: $max_treedepth"
        echo "Cores: $cores"
        echo "Memory: 32GB"
        ;;
    4) # Maximum Performance
        chains=16
        iterations=5000
        adapt_delta=0.99
        max_treedepth=15
        cores=32
        memory="64.GB"
        flag="maximum"
        echo ""
        echo "Maximum Performance Profile Selected"
        echo "Chains: $chains"
        echo "Iterations: $iterations"
        echo "Adapt delta: $adapt_delta"
        echo "Max treedepth: $max_treedepth"
        echo "Cores: $cores"
        echo "Memory: 64GB"
        ;;
    5) # Stable Profile - Ultra-conservative settings that actually work
        chains=2
        iterations=200
        adapt_delta=0.8
        max_treedepth=8
        cores=2
        memory="8.GB"
        flag="stable"
        echo ""
        echo "Stable Profile Selected (RECOMMENDED)"
        echo "Chains: $chains"
        echo "Iterations: $iterations"
        echo "Adapt delta: $adapt_delta"
        echo "Max treedepth: $max_treedepth"
        echo "Cores: $cores"
        echo "Memory: 8GB"
        echo "This profile uses ultra-conservative settings proven to work without Stan crashes."
        ;;
    6) # Custom Configuration
        echo ""
        echo "======== Custom Resource Configuration ========"
        echo "IMPORTANT: For optimal performance in HPC environments, ensure:"
        echo "- Cores should be divisible by chains for optimal parallelization"
        echo "- Memory allocation should account for Stan's overhead (~0.5-1GB per chain)"
        echo "- Consider setting higher adapt_delta values (0.95-0.99) for complex models"
        echo ""
        
        # Get MCMC parameters
        echo "Recommended settings based on pathogen count ($pathogen_count pathogens):"
        echo ""
        echo "IMPORTANT: Each pathogen gets the FULL number of chains you specify!"
        echo "Total chains = pathogens × chains per pathogen"
        echo ""
        
        # Base recommendations on pathogen count - INVERSE relationship
        if [[ $pathogen_count -le 1 ]]; then
            echo "Single pathogen: Can use more chains for better convergence"
            echo "Recommended: 8-16 chains, 5000 iterations, 16 cores"
            default_chains=8
            default_iterations=5000
            default_cores=16
        elif [[ $pathogen_count -le 3 ]]; then
            echo "Multiple pathogens: Moderate chains to balance quality and resources"
            echo "Recommended: 4-8 chains, 3000 iterations, 24 cores"
            echo "Total chains will be: $pathogen_count pathogens × chains = $(($pathogen_count * 4))-$(($pathogen_count * 8)) chains"
            default_chains=4
            default_iterations=3000
            default_cores=24
        else
            echo "Many pathogens: Reduce chains per pathogen to avoid resource exhaustion"
            echo "Recommended: 2-4 chains, 2000 iterations, 32 cores"
            echo "Total chains will be: $pathogen_count pathogens × chains = $(($pathogen_count * 2))-$(($pathogen_count * 4)) chains"
            default_chains=2
            default_iterations=2000
            default_cores=32
        fi
        
        read -p "Number of chains [$default_chains]: " chains
        chains=${chains:-$default_chains}
        
        read -p "Number of iterations [$default_iterations]: " iterations
        iterations=${iterations:-$default_iterations}
        
        read -p "Adapt delta (0.0-1.0) [0.99]: " adapt_delta
        adapt_delta=${adapt_delta:-0.99}
        
        read -p "Max treedepth [15]: " max_treedepth
        max_treedepth=${max_treedepth:-15}
        
        read -p "Number of cores [$default_cores]: " cores
        cores=${cores:-$default_cores}
        
        # Calculate optimal memory based on cores, chains, and iterations
        optimal_memory=$((cores * 2))
        chains_memory=$((chains * 1))
        # Use bc if available, otherwise use shell arithmetic
        if command -v bc >/dev/null 2>&1; then
            iterations_factor=$(echo "scale=2; ${iterations}/1000" | bc 2>/dev/null || echo "1")
            iterations_memory=$(echo "scale=0; ${iterations_factor} * 4" | bc 2>/dev/null || echo "4")
        else
            # Fallback to shell arithmetic (less precise but works everywhere)
            iterations_factor=$((iterations / 1000))
            iterations_factor=${iterations_factor:-1}
            iterations_memory=$((iterations_factor * 4))
        fi
        
        # Ensure iterations_memory is a number even if bc fails
        if ! [[ "$iterations_memory" =~ ^[0-9]+$ ]]; then
            iterations_memory=4
        fi
        
        # Use the larger of the memory calculations with a minimum of cores*2
        suggested_memory=$((optimal_memory > chains_memory ? optimal_memory : chains_memory))
        suggested_memory=$((suggested_memory > iterations_memory ? suggested_memory : iterations_memory))
        
        echo "Memory recommendation: Based on cores, chains, and iterations, recommend at least ${suggested_memory}GB"
        read -p "Memory in GB [${suggested_memory}]: " memory_gb
        memory_gb=${memory_gb:-$suggested_memory}
        memory="${memory_gb}.GB"
        
        flag="custom"
        ;;
    *)
        echo "Invalid selection. Using Standard Run profile."
        chains=4
        iterations=1000
        adapt_delta=0.95
        max_treedepth=10
        cores=8
        memory="16.GB"
        flag="standard"
        ;;
esac

# Get resume option
echo ""
echo "Resume previous failed run?"
echo "1) No, start fresh"
echo "2) Yes, resume from last successful step"
read -p "Enter selection [1]: " resume_choice
resume_choice=${resume_choice:-1}

resume_flag=""
if [[ "$resume_choice" == "2" ]]; then
    resume_flag="-resume"
    echo "Run will resume from last successful step"
fi

# Get background execution preference
echo ""
read -p "Run in background? (y/n) [n]: " bg_choice
bg_choice=${bg_choice:-n}
if [[ "$bg_choice" =~ ^[Yy]$ ]]; then
    background=true
else
    background=false
fi

# Get output directory
echo ""
read -p "Output directory [${outDir}]: " user_outdir
outDir=${user_outdir:-$outDir}
# Trim any leading/trailing spaces
outDir=$(echo "$outDir" | xargs)

# Dashboard preferences - simplified and automatic
echo ""
echo "======== Dashboard Generation ========"
echo "An interactive HTML dashboard will be automatically generated with run details included."

# Automatically generate dashboard title with timestamp for identification
dashboard_title="FoodNetTrends Analysis - Run ${projID}"
dashboard_params="--enable_dashboard true --dashboard_title \"$dashboard_title\""

# Tell the user what's happening
echo "Dashboard will be created with title: \"$dashboard_title\""
echo "Project ID: ${projID} (use this if you need to regenerate the dashboard later)"

# Provide instructions for fallback dashboard generation if needed
echo ""
echo "NOTE: If dashboard generation fails, you can regenerate it after completion with:"
echo "  ./dashboard.sh ${projID}"

# Add warning about resource usage for heavy analyses
if [[ $pathogen_count -gt 4 || $chains -gt 8 || $iterations -gt 3000 ]]; then
    echo ""
    echo "====== PERFORMANCE WARNING ======"
    echo "You have selected a resource-intensive configuration:"
    echo "- Pathogens: $pathogen_count"
    echo "- Chains: $chains"
    echo "- Iterations: $iterations"
    echo "- Memory: $memory"
    echo ""
    echo "This analysis may take considerable time and resources."
    echo "Consider submitting as a background job and checking logs periodically."
    echo "==============================="
    
    # Default to background mode for very heavy analyses
    if [[ $background == false && ($chains -gt 12 || $iterations -gt 5000) ]]; then
        echo "Recommending background execution for this heavy analysis."
        read -p "Run in background? (y/n) [y]: " bg_recommendation
        bg_recommendation=${bg_recommendation:-y}
        if [[ "$bg_recommendation" =~ ^[Yy]$ ]]; then
            background=true
            echo "Running in background mode enabled."
        fi
    fi
fi

# Generate unique project identifier for output organization
# Users can override this by setting projID before running the script
if [[ -z "$projID" ]]; then
    projID="$timestamp"
    echo "Using auto-generated project ID: $projID"
else
    echo "Using user-specified project ID: $projID"
fi

# Build the command
cmd="nextflow run main.nf -profile singularity,production"
if [[ -n "$resume_flag" ]]; then
    cmd="$cmd $resume_flag"
fi

# Add HPC configuration optimizations
cmd="$cmd -process.memory $memory"
cmd="$cmd -process.cpus $cores"
cmd="$cmd -executor.queueSize 100"
cmd="$cmd -executor.submitRateLimit '10/1min'"

# Only include census file parameters if they have non-empty values
if [[ -n "${censusFileB}" && "${censusFileB}" != "true" ]]; then
    cmd="$cmd --censusFileB \"${censusFileB}\""
fi

if [[ -n "${censusFileP}" && "${censusFileP}" != "true" ]]; then
    cmd="$cmd --censusFileP \"${censusFileP}\""
fi

cmd="$cmd --travel \"${travel}\""
cmd="$cmd --cidt \"${cidt}\""
cmd="$cmd --iterations ${iterations}"
cmd="$cmd --chains ${chains}"
cmd="$cmd --adapt_delta ${adapt_delta}"
cmd="$cmd --max_treedepth ${max_treedepth}"
cmd="$cmd --cores ${cores}"
cmd="$cmd --seed 123"
cmd="$cmd --outdir \"${outDir}\""
cmd="$cmd --projID \"${projID}\""
cmd="$cmd --pathogen \"${pathogens}\""
cmd="$cmd --states \"${states}\""

# Add MMWR file parameter based on workflow mode
if [[ "$workflow_mode" == "3" ]]; then
    # Using preprocessed data
    cmd="$cmd --mmwrFile \"$mmwrFile\" --preprocessed true --cleanFile \"$mmwrFile\""
    
    # Add metadata file if available
    if [[ -n "$preprocessed_metadata" ]]; then
        cmd="$cmd --metadata \"$preprocessed_metadata\""
    fi
else
    # Using raw data
    cmd="$cmd --mmwrFile \"$mmwrFile\""
    
    # If this is workflow_mode 1, we already preprocessed
    if [[ "$workflow_mode" == "1" ]]; then
        cmd="$cmd --preprocessed true --cleanFile \"$mmwrFile\""
        if [[ -n "$preprocessed_metadata" ]]; then
            cmd="$cmd --metadata \"$preprocessed_metadata\""
        fi
    fi
fi

# Add dashboard parameters
cmd="$cmd $dashboard_params"

# Add serotype and serogroup parameters to command BEFORE background wrapper
if [[ -n "$stec_serogroups" ]]; then
    cmd="$cmd --stec_serogroups \"$stec_serogroups\""
fi
if [[ -n "$salmonella_serotypes" ]]; then
    cmd="$cmd --salmonella_serotypes \"$salmonella_serotypes\""
fi

# Add background option if needed (AFTER all parameters are added)
if [[ "$background" == true ]]; then
    # Create logs directory if it doesn't exist
    mkdir -p logs
    log_file="logs/foodnet_run_${timestamp}.log"
    cmd="nohup $cmd > \"${log_file}\" 2>&1 &"
    echo "Process will run in background with log: ${log_file}"
fi

# Analysis levels: STEC uses serogroups, Salmonella uses serotypes

# Display summary
# Show configuration summary and save
show_configuration_summary

echo "Save this configuration for future use?"
read -p "Save configuration? (y/n) [y]: " save_config
save_config=${save_config:-y}

if [[ "$save_config" =~ ^[Yy]$ ]]; then
    save_configuration
fi

echo ""
read -p "Proceed with this configuration? (y/n) [y]: " proceed
proceed=${proceed:-y}

if [[ ! "$proceed" =~ ^[Yy]$ ]]; then
    echo "Analysis cancelled by user."
    exit 0
fi

echo ""
echo "========= Analysis Summary ==========="
echo "Mode: $flag"
if [[ -n "$resume_flag" ]]; then
    echo "Resume previous run: Yes"
else
    echo "Resume previous run: No"
fi

if [[ "$workflow_mode" == "1" || "$workflow_mode" == "3" ]]; then
    echo "Using preprocessed data: Yes"
else
    echo "Using preprocessed data: No"
fi

echo "Pathogens: $pathogens"
echo "States: $states"
echo "Travel status: $travel"
echo "CIDT/culture methods: $cidt"

echo "Input files:"
echo "  Data file: $mmwrFile"
echo "  Census file (bacterial): $censusFileB" 
echo "  Census file (parasitic): $censusFileP"
if [[ -n "$preprocessed_metadata" ]]; then
    echo "  Metadata: $preprocessed_metadata"
fi

echo "Chains: $chains"
echo "Iterations: $iterations"
echo "Adapt delta: $adapt_delta"
echo "Max treedepth: $max_treedepth"
echo "Cores: $cores"
echo "Memory: $memory"

echo "Output directory: $outDir"
echo "Run in background: $([ "$background" == true ] && echo "Yes" || echo "No")"
echo ""
echo "Command to run:"
# Clean command for display (remove extra spaces and newlines)
display_cmd=$(echo "$cmd" | tr '\n' ' ' | sed 's/  */ /g')
echo "$display_cmd"
echo ""

# If we're just getting the command, print it and exit
if [[ "$GET_COMMAND_ONLY" == "true" ]]; then
    echo "$cmd"
    exit 0
fi

# Get confirmation from user
read -p "Execute command? (y/n) [y]: " execute
execute=${execute:-y}

if [[ "$execute" =~ ^[Yy]$ ]]; then
    # Final HPC resource optimization check
    echo ""
    echo "========== HPC Resource Validation =========="
    
    # Check if cores are divisible by chains for optimal performance
    if [ $((cores % chains)) -ne 0 ]; then
        echo "WARNING: Number of cores ($cores) is not divisible by chains ($chains)."
        echo "For optimal parallelization, consider using a multiple of chains."
        suggested_cores=$((chains * (cores / chains + 1)))
        if [ $suggested_cores -le $cores ]; then
            suggested_cores=$cores
        fi
        echo "Suggested correction: $suggested_cores cores"
        
        read -p "Adjust cores to $suggested_cores? (y/n) [y]: " adjust_cores
        adjust_cores=${adjust_cores:-y}
        if [[ "$adjust_cores" =~ ^[Yy]$ ]]; then
            cores=$suggested_cores
            cmd=$(echo "$cmd" | sed "s/--cores [0-9]*/--cores $cores/")
            echo "Adjusted cores to $cores"
        fi
    else
        echo "OPTIMAL: Cores ($cores) are perfectly divisible by chains ($chains)."
    fi
    
    # Check for optimal memory-to-core ratio on HPC
    # Use bc if available, otherwise use shell arithmetic
    if command -v bc >/dev/null 2>&1; then
        mem_per_core=$(echo "scale=2; ${memory_gb}/${cores}" | bc 2>/dev/null || echo "0")
        low_memory=$(echo "$mem_per_core < 1.5" | bc -l 2>/dev/null || echo "0")
        if [[ "$low_memory" == "1" ]]; then
            echo "⚠️  WARNING: Low memory per core ratio ($mem_per_core GB/core)"
            echo "   This may cause Stan initialization failures or slow performance"
            echo "   HPC environments typically perform best with 2-4GB per core"
            suggested_memory=$((cores * 2))
            echo "   Suggested memory: ${suggested_memory}GB"
            
            read -p "Adjust memory to ${suggested_memory}GB? (y/n) [y]: " adjust_memory
            adjust_memory=${adjust_memory:-y}
            if [[ "$adjust_memory" =~ ^[Yy]$ ]]; then
                memory_gb=$suggested_memory
                memory="${memory_gb}.GB"
                cmd=$(echo "$cmd" | sed "s/-process.memory [0-9]*\.[G|M]B/-process.memory $memory/")
                echo "Adjusted memory to $memory"
            fi
        else
            echo "OPTIMAL: Memory-to-core ratio is good (${mem_per_core}GB per core)."
        fi
    else
        # Fallback for systems without bc
        mem_per_core_x10=$((memory_gb * 10 / cores))
        if (( mem_per_core_x10 < 15 )); then
            echo "⚠️  WARNING: Low memory per core ratio"
            echo "   This may cause Stan initialization failures or slow performance"
            echo "   HPC environments typically perform best with 2-4GB per core"
            suggested_memory=$((cores * 2))
            echo "   Suggested memory: ${suggested_memory}GB"
            
            read -p "Adjust memory to ${suggested_memory}GB? (y/n) [y]: " adjust_memory
            adjust_memory=${adjust_memory:-y}
            if [[ "$adjust_memory" =~ ^[Yy]$ ]]; then
                memory_gb=$suggested_memory
                memory="${memory_gb}.GB"
                cmd=$(echo "$cmd" | sed "s/-process.memory [0-9]*\.[G|M]B/-process.memory $memory/")
                echo "Adjusted memory to $memory"
            fi
        else
            echo "OPTIMAL: Memory-to-core ratio is good."
        fi
    fi
    
    echo "Starting analysis..."
    echo "$(date): Executing command: $cmd" >> "$error_log"
    
    # Clean up command string to prevent eval issues
    cmd=$(echo "$cmd" | tr '\n' ' ' | sed 's/  */ /g')
    
    # Execute command using bash -c instead of eval to handle complex parameter strings
    if ! bash -c "$cmd"; then
        echo "Error running analysis command."
        echo "Check .nextflow.log for details."
        echo "$(date): Command execution failed" >> "$error_log"
        exit 1
    fi
    
    # Skip validation for background jobs
    if [[ "$background" == true ]]; then
        echo ""
        echo "======== Analysis Running in Background ========"
        echo "Process is running in background. Check log file for progress:"
        echo "  $log_file"
        echo ""
        echo "To monitor progress:"
        echo "  tail -f $log_file"
        echo ""
        echo "To check if still running:"
        echo "  ps aux | grep nextflow"
        echo ""
    else
        # Validate outputs after successful execution (foreground only)
        echo ""
        echo "======== Output Validation ========"
    echo "Output directory structure:"
    ls -la "$outDir" 2>/dev/null || echo "Output directory does not exist: $outDir"
    
    # Construct full output path using project ID subdirectory
    outDir="$outDir/$projID"
    echo "Validating outputs in: $outDir"
    ls -la "$outDir" 2>/dev/null || echo "Output directory does not exist: $outDir"
    
    # Parse pathogens to validate each one
    IFS=',' read -ra PATHOGEN_LIST <<< "$pathogens"
    validation_passed=true
    
    for pathogen_name in "${PATHOGEN_LIST[@]}"; do
        pathogen_name=$(echo "$pathogen_name" | tr -d ' ')  # Remove whitespace
        
        # Validate outputs for this pathogen
        if ! validate_analysis_outputs "$pathogen_name" "$outDir"; then
            validation_passed=false
        fi
        
        # Verify spline trend visualizations were generated
        verify_spline_trends "$pathogen_name" "$outDir"
    done
    
    if [[ "$validation_passed" == true ]]; then
        echo ""
        echo "🎉 Analysis completed successfully with all expected outputs!"
        echo "📊 Smooth trend visualizations should resolve the 'spikey graph' issue."
    else
        echo ""
        echo "⚠ Analysis completed but some outputs may be missing or incomplete."
        echo "Check the validation messages above for details."
    fi
    
    # Apply standard file permissions for HPC compatibility
    # All pipeline outputs must have 755 permissions per infrastructure requirements
    echo ""
    echo "Setting file permissions to 755 for HPC compatibility..."
    find "$outDir" -type f -exec chmod 755 {} + 2>/dev/null
    echo "✓ File permissions updated"
    fi  # End of background check
else
    echo "Execution canceled."
    echo "$(date): User canceled execution" >> "$error_log"
fi