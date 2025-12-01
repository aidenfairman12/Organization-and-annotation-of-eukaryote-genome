#!/usr/bin/env Rscript

# Install GENESPACE and dependencies
# Run this script once before running the GENESPACE analysis

cat("========================================\n")
cat("Installing GENESPACE and dependencies\n")
cat("========================================\n\n")

# Set repository
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Check if BiocManager is installed
if (!require("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}

# Check if devtools is installed (needed for GitHub installation)
if (!require("devtools", quietly = TRUE)) {
  cat("Installing devtools...\n")
  install.packages("devtools")
}

# Install critical dependencies first (igraph often needs special handling)
cat("\nInstalling critical dependencies...\n")
critical_pkgs <- c("igraph", "ggplot2", "data.table", "dbscan", "R.utils")
for (pkg in critical_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

# Install required Bioconductor packages
cat("\nInstalling Bioconductor dependencies...\n")
BiocManager::install(c("Biostrings", "rtracklayer"), update = FALSE, ask = FALSE)

# Install GENESPACE from GitHub
cat("\nInstalling GENESPACE from GitHub...\n")
devtools::install_github("jtlovell/GENESPACE", upgrade = "never", dependencies = TRUE)

# Verify installation
cat("\n========================================\n")
cat("Verifying installation...\n")
cat("========================================\n\n")

if (require("GENESPACE", quietly = TRUE)) {
  cat("✓ GENESPACE successfully installed!\n")
  cat("  Version:", as.character(packageVersion("GENESPACE")), "\n")
} else {
  cat("✗ GENESPACE installation failed\n")
  cat("  Please check the error messages above\n")
}

cat("\nYou can now run the GENESPACE analysis:\n")
cat("  Rscript scripts/11_run_genespace.R results/genespace\n\n")
