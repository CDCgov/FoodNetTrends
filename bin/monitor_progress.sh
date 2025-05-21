#!/usr/bin/env bash
# =========================================================================
# FoodNet Trends Pathogen Analysis Progress Monitor
# =========================================================================
#
# This script provides a user-friendly interface for monitoring the progress
# of pathogen analysis jobs, regardless of environment (HPC, local, etc.)
#
# The script displays real-time progress information by reading progress files
# created by the analysis pipeline, providing a clear view of current status.
#
# Usage:
#   ./monitor_progress.sh [options]
#
# Options:
#   -d, --dir DIR       Directory to monitor for progress files (default: current dir)
#   -p, --pathogen NAME Only show progress for specific pathogen
#   -f, --follow        Continuously update (like 'tail -f')
#   -i, --interval SEC  Refresh interval in seconds (default: 5)
#   -h, --help          Show this help message
#
# =========================================================================

# Default settings
MONITOR_DIR="."
SPECIFIC_PATHOGEN=""
FOLLOW_MODE=false
REFRESH_INTERVAL=5
USE_COLORS=1  # Default to using colors

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -d|--dir)
      MONITOR_DIR="$2"
      shift 2
      ;;
    -p|--pathogen)
      SPECIFIC_PATHOGEN="$2"
      shift 2
      ;;
    -f|--follow)
      FOLLOW_MODE=true
      shift
      ;;
    -i|--interval)
      REFRESH_INTERVAL="$2"
      shift 2
      ;;
    --no-color)
      USE_COLORS=0
      # Redefine all color variables to be empty
      RED=''
      GREEN=''
      YELLOW=''
      BLUE=''
      MAGENTA=''
      CYAN=''
      GRAY=''
      BOLD=''
      RESET=''
      shift
      ;;
    -h|--help)
      echo -e "${BOLD}FoodNet Trends Progress Monitor${RESET}"
      echo
      echo "Usage: $0 [options]"
      echo
      echo "Options:"
      echo "  -d, --dir DIR       Directory to monitor for progress files (default: current dir)"
      echo "  -p, --pathogen NAME Only show progress for specific pathogen"
      echo "  -f, --follow        Continuously update (like 'tail -f')"
      echo "  -i, --interval SEC  Refresh interval in seconds (default: 5)"
      echo "  --no-color          Disable colored output (useful for non-interactive terminals)"
      echo "  -h, --help          Show this help message"
      echo
      echo "Example:"
      echo "  $0 --dir /path/to/results --follow"
      echo "  $0 --pathogen SALMONELLA --follow --interval 10"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use -h or --help for usage information"
      exit 1
      ;;
  esac
done

# Check if monitor directory exists
if [ ! -d "$MONITOR_DIR" ]; then
  echo -e "${RED}Error: Directory '$MONITOR_DIR' not found${RESET}"
  exit 1
fi

# ANSI color codes for prettier output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
GRAY='\033[0;37m'
BOLD='\033[1m'
RESET='\033[0m'

# Check if we're in a terminal that supports colors
if [ -t 1 ]; then
  USE_COLORS=1
else
  USE_COLORS=0
  # Redefine all color variables to be empty if not in a terminal
  RED=''
  GREEN=''
  YELLOW=''
  BLUE=''
  MAGENTA=''
  CYAN=''
  GRAY=''
  BOLD=''
  RESET=''
fi

# Function to get progress color based on percentage
get_progress_color() {
  local percent=$1
  if [ $percent -lt 25 ]; then
    echo -e "${RED}"
  elif [ $percent -lt 50 ]; then
    echo -e "${YELLOW}"
  elif [ $percent -lt 75 ]; then
    echo -e "${CYAN}"
  else
    echo -e "${GREEN}"
  fi
}

# Function to draw a progress bar
draw_progress_bar() {
  local percent=$1
  local width=40
  local filled=$((percent * width / 100))
  local empty=$((width - filled))
  
  # Progress bar color based on completion
  local color=$(get_progress_color $percent)
  
  printf "${color}["
  printf "%0.s#" $(seq 1 $filled)
  printf "%0.s." $(seq 1 $empty)
  printf "] %3d%%${RESET}" $percent
}

# Function to display progress information
display_progress() {
  local clear_screen=true
  if [ "$1" == "no_clear" ]; then
    clear_screen=false
  fi
  
  if $clear_screen; then
    clear
    echo -e "${BOLD}FoodNet Trends Pathogen Analysis Progress Monitor${RESET}"
    echo -e "${GRAY}Monitoring directory: $MONITOR_DIR${RESET}"
    echo -e "${GRAY}Press Ctrl+C to exit${RESET}"
    echo "════════════════════════════════════════════════════════════════════════"
  fi
  
  # Find all progress files - use safer approach with pathname expansion
  local progress_files=()
  if [ -n "$SPECIFIC_PATHOGEN" ]; then
    # Look for specific pathogen progress file
    if [ -f "$MONITOR_DIR/${SPECIFIC_PATHOGEN}_progress.txt" ]; then
      progress_files=("$MONITOR_DIR/${SPECIFIC_PATHOGEN}_progress.txt")
    fi
  else
    # Find all progress files using pathname expansion instead of find
    # This is safer across different environments
    shopt -s nullglob  # Handle case where no files match
    progress_files=("$MONITOR_DIR"/*_progress.txt)
    shopt -u nullglob
    
    # Also look in subdirectories (some work dirs have nested structure)
    if [ "${#progress_files[@]}" -eq 0 ]; then
      echo "No progress files in current directory, checking subdirectories..."
      # Try progressive deeper searches until we find something
      for depth in 2 3 4 5; do
        echo "Searching with depth $depth..."
        progress_files=($(find "$MONITOR_DIR" -maxdepth $depth -name "*_progress.txt" 2>/dev/null | sort))
        if [ "${#progress_files[@]}" -gt 0 ]; then
          echo "Found ${#progress_files[@]} progress files at depth $depth"
          break
        fi
      done
      
      # If still not found, try alternate directories like /tmp and work directories
      if [ "${#progress_files[@]}" -eq 0 ]; then
        echo "Checking alternate locations..."
        for alt_dir in "$MONITOR_DIR/work" "$MONITOR_DIR/*/work" "$MONITOR_DIR/*/*/work" "/tmp"; do
          if [ -d "$alt_dir" ]; then
            echo "Checking $alt_dir..."
            alt_files=($(find "$alt_dir" -maxdepth 2 -name "*_progress.txt" 2>/dev/null | sort))
            if [ "${#alt_files[@]}" -gt 0 ]; then
              progress_files=("${alt_files[@]}")
              echo "Found ${#progress_files[@]} progress files in $alt_dir"
              break
            fi
          fi
        done
      fi
    fi
  fi
  
  # Check if we found any progress files
  if [ ${#progress_files[@]} -eq 0 ]; then
    echo -e "${YELLOW}No progress files found.${RESET}"
    if [ -n "$SPECIFIC_PATHOGEN" ]; then
      echo -e "${GRAY}Looking for: ${SPECIFIC_PATHOGEN}_progress.txt${RESET}"
    else
      echo -e "${GRAY}Looking for: *_progress.txt${RESET}"
    fi
    return
  fi
  
  # Process each progress file
  for progress_file in "${progress_files[@]}"; do
    # Extract pathogen name from filename
    local pathogen=$(basename "$progress_file" | sed 's/_progress.txt//')
    
    # Read progress data
    local pathogen_data=$(grep "PATHOGEN:" "$progress_file" 2>/dev/null | cut -d':' -f2- | tr -d ' ')
    local stage=$(grep "STAGE:" "$progress_file" 2>/dev/null | cut -d':' -f2- | tr -d ' ')
    local progress=$(grep "PROGRESS:" "$progress_file" 2>/dev/null | cut -d':' -f2- | tr -d ' ' | tr -d '%')
    local message=$(grep "MESSAGE:" "$progress_file" 2>/dev/null | cut -d':' -f2-)
    local elapsed=$(grep "ELAPSED:" "$progress_file" 2>/dev/null | cut -d':' -f2- | tr -d ' ')
    local remaining=$(grep "REMAINING:" "$progress_file" 2>/dev/null | cut -d':' -f2- | tr -d ' ')
    local timestamp=$(grep "TIMESTAMP:" "$progress_file" 2>/dev/null | cut -d':' -f2-)
    
    # Format and display
    echo -e "\n${BOLD}${BLUE}$pathogen${RESET}"
    echo -e "${GRAY}Last updated: $timestamp${RESET}"
    echo -e "${MAGENTA}Stage:${RESET} $stage"
    echo -e "$(draw_progress_bar $progress)"
    echo -e "${GRAY}Elapsed: $elapsed | Remaining: $remaining${RESET}"
    
    if [ "$message" != " -" ]; then
      echo -e "${YELLOW}$message${RESET}"
    fi
    
    # Check for log file to show recent activity
    local log_file="$MONITOR_DIR/${pathogen}_progress_log.txt"
    if [ -f "$log_file" ]; then
      echo -e "\n${GRAY}Recent activity:${RESET}"
      tail -n 5 "$log_file" | grep -v "===========" | grep -v "PATHOGEN:" | sed 's/^/  /'
    fi
    
    # Look for output files to verify progress
    local ir_files=$(find "$MONITOR_DIR" -name "${pathogen}_*.csv" 2>/dev/null | wc -l)
    local fig_files=$(find "$MONITOR_DIR" -name "${pathogen}_*.png" 2>/dev/null | wc -l)
  
    echo -e "\n${GRAY}Files generated:${RESET}"
    echo -e "  Data files: $ir_files"
    echo -e "  Figures: $fig_files"
    
    echo "────────────────────────────────────────────────────────────────────────"
  done
  
  # Also check for simple percent files (for processes not using full progress tracking)
  if [ -z "$SPECIFIC_PATHOGEN" ]; then
    local percent_files=($(find "$MONITOR_DIR" -maxdepth 1 -name "*_percent.txt" 2>/dev/null | sort))
    
    for percent_file in "${percent_files[@]}"; do
      local pathogen=$(basename "$percent_file" | sed 's/_percent.txt//')
      
      # Skip if we already processed the full progress file
      if [[ " ${progress_files[@]} " =~ " $MONITOR_DIR/${pathogen}_progress.txt " ]]; then
        continue
      fi
      
      # Read percentage
      local progress=$(cat "$percent_file" 2>/dev/null | tr -d ' \n\r\t')
      
      if [ -n "$progress" ]; then
        echo -e "\n${BOLD}${BLUE}$pathogen${RESET} (simple tracking)"
        echo -e "$(draw_progress_bar $progress)"
      fi
    done
  fi
}

# Main execution
if $FOLLOW_MODE; then
  # Continuous monitoring mode
  while true; do
    display_progress
    sleep $REFRESH_INTERVAL
  done
else
  # One-time display
  display_progress
fi