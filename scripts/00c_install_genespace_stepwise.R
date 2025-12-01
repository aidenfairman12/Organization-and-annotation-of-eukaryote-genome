#!/usr/bin/env Rscript

# Step-by-step installation focusing on igraph compilation
cat("========================================\n")
cat("GENESPACE Installation - Step by Step\n")
cat("========================================\n\n")

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Step 1: Install basic tools
cat("Step 1: Installing basic tools...\n")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
cat("✓ BiocManager installed\n\n")

# Step 2: Install igraph dependencies first
cat("Step 2: Installing igraph dependencies...\n")
deps <- c("pkgconfig", "Matrix", "isoband", "cpp11", "cli")
for (pkg in deps) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    install.packages(pkg)
  }
}
cat("✓ Dependencies installed\n\n")

# Step 3: Try installing igraph with minimal configuration
cat("Step 3: Installing igraph (this may take 5-10 minutes)...\n")
cat("  Attempting installation with system libraries...\n")

# Set environment for compilation
Sys.setenv(MAKEFLAGS = paste0("-j", parallel::detectCores()))

if (!require("igraph", quietly = TRUE)) {
  # Try installation
  install.packages("igraph", 
                  configure.args = c("--with-external-glpk=no",
                                    "--with-external-blas=yes",
                                    "--with-external-lapack=yes"))
  
  # Check if it worked
  if (!require("igraph", quietly = TRUE)) {
    cat("\n✗ igraph installation failed\n")
    cat("\nTrying alternative approach without external libraries...\n")
    install.packages("igraph")
  }
}

if (require("igraph", quietly = TRUE)) {
  cat("✓ igraph successfully installed!\n")
  cat("  Version:", as.character(packageVersion("igraph")), "\n\n")
} else {
  cat("✗ igraph installation failed completely\n")
  cat("  Cannot proceed with GENESPACE installation\n")
  cat("\nPlease contact system administrator to install system libraries:\n")
  cat("  - libxml2-dev\n")
  cat("  - libglpk-dev\n")
  cat("  - libgmp3-dev\n")
  quit(save = "no", status = 1)
}

# Step 4: Install other CRAN dependencies
cat("Step 4: Installing other CRAN dependencies...\n")
other_deps <- c("data.table", "dbscan", "R.utils", "ggplot2", "devtools")
for (pkg in other_deps) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    install.packages(pkg)
  }
}
cat("✓ CRAN dependencies installed\n\n")

# Step 5: Install Bioconductor packages
cat("Step 5: Installing Bioconductor packages...\n")
BiocManager::install(c("Biostrings", "rtracklayer"), update = FALSE, ask = FALSE)
cat("✓ Bioconductor packages installed\n\n")

# Step 6: Install GENESPACE from GitHub
cat("Step 6: Installing GENESPACE from GitHub...\n")
if (!require("GENESPACE", quietly = TRUE)) {
  devtools::install_github("jtlovell/GENESPACE", 
                          upgrade = "never",
                          dependencies = TRUE)
}

# Final verification
cat("\n========================================\n")
cat("Installation Complete - Verification\n")
cat("========================================\n\n")

all_ok <- TRUE
packages_to_check <- c("igraph", "data.table", "Biostrings", "GENESPACE")

for (pkg in packages_to_check) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    ver <- packageVersion(pkg)
    cat(sprintf("✓ %-15s %s\n", pkg, ver))
  } else {
    cat(sprintf("✗ %-15s FAILED\n", pkg))
    all_ok <- FALSE
  }
}

cat("\n")
if (all_ok) {
  cat("SUCCESS! All packages installed.\n")
  cat("\nYou can now run:\n")
  cat("  Rscript scripts/11_run_genespace.R results/genespace\n\n")
} else {
  cat("FAILED: Some packages could not be installed.\n\n")
  quit(save = "no", status = 1)
}
