#!/usr/bin/env Rscript

# Install GENESPACE and dependencies
# Alternative approach: install from pre-compiled binaries when possible

cat("========================================\n")
cat("Installing GENESPACE and dependencies\n")
cat("Alternative method (binary packages)\n")
cat("========================================\n\n")

# Set repository - use binary packages when available
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Function to try installing a package
try_install <- function(pkg, method = "binary") {
  cat("Installing", pkg, "...\n")
  if (method == "binary") {
    # Try binary first (faster, fewer compilation issues)
    tryCatch({
      install.packages(pkg, type = "binary", dependencies = TRUE)
      return(TRUE)
    }, error = function(e) {
      cat("  Binary not available, trying source...\n")
    })
  }
  
  # Try from source
  tryCatch({
    install.packages(pkg, type = "source", dependencies = TRUE)
    return(TRUE)
  }, error = function(e) {
    cat("  ERROR installing", pkg, ":", e$message, "\n")
    return(FALSE)
  })
}

# Install BiocManager
if (!require("BiocManager", quietly = TRUE)) {
  try_install("BiocManager")
}

# Install devtools
if (!require("devtools", quietly = TRUE)) {
  try_install("devtools")
}

# Install igraph specifically (this is the problematic one)
cat("\n=== Installing igraph (may take a few minutes) ===\n")
if (!require("igraph", quietly = TRUE)) {
  # Try with extra configuration for compilation
  Sys.setenv(MAKEFLAGS = "-j4")  # Use 4 cores for compilation
  try_install("igraph", method = "source")
}

# Install other critical dependencies
cat("\n=== Installing other dependencies ===\n")
critical_pkgs <- c("ggplot2", "data.table", "dbscan", "R.utils", "Rcpp")
for (pkg in critical_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    try_install(pkg)
  } else {
    cat("✓", pkg, "already installed\n")
  }
}

# Install Bioconductor packages
cat("\n=== Installing Bioconductor dependencies ===\n")
if (!require("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings", update = FALSE, ask = FALSE)
}
if (!require("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer", update = FALSE, ask = FALSE)
}

# Install GENESPACE from GitHub
cat("\n=== Installing GENESPACE from GitHub ===\n")
if (!require("GENESPACE", quietly = TRUE)) {
  devtools::install_github("jtlovell/GENESPACE", 
                          upgrade = "never", 
                          dependencies = TRUE,
                          force = FALSE)
}

# Verify installation
cat("\n========================================\n")
cat("Verifying installation...\n")
cat("========================================\n\n")

success <- TRUE

# Check critical packages
check_pkgs <- c("igraph", "data.table", "Biostrings", "GENESPACE")
for (pkg in check_pkgs) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "installed\n")
  } else {
    cat("✗", pkg, "FAILED\n")
    success <- FALSE
  }
}

if (success) {
  cat("\n✓ All packages successfully installed!\n")
  cat("\nYou can now run:\n")
  cat("  Rscript scripts/11_run_genespace.R results/genespace\n\n")
} else {
  cat("\n✗ Some packages failed to install.\n")
  cat("  Please check error messages above.\n\n")
}
