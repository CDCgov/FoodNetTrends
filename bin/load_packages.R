# List of packages
packages <- c(
  "doParallel", "iterators", "foreach",
  "HDInterval", "gtools", "haven",
  "argparse", "stringr", "dplyr",
  "purrr", "readr", "tidyr",
  "tibble", "ggplot2", "tidyverse",
  "tidybayes", "brms", "Rcpp"
)

# Install missing packages
install_if_missing <- function(packages) {
  installed <- rownames(installed.packages())
  to_install <- packages[!(packages %in% installed)]
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cran.r-project.org")
  }
}

# Run the function
install_if_missing(packages)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

cat("All packages loaded successfully.\n")

