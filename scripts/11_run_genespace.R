library(GENESPACE)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get the folder where the genespace workingDirectory is located
if (length(args) == 0) {
  stop("Please provide the working directory path as an argument")
}

wd <- args[1]
cat("Working directory:", wd, "\n")

# Check that required directories exist
if (!dir.exists(file.path(wd, "bed"))) {
  stop("bed/ directory not found in working directory")
}
if (!dir.exists(file.path(wd, "peptide"))) {
  stop("peptide/ directory not found in working directory")
}

# List files for verification
cat("\nBED files found:\n")
print(list.files(file.path(wd, "bed"), pattern = "\\.bed$"))

cat("\nPeptide files found:\n")
print(list.files(file.path(wd, "peptide"), pattern = "\\.fa$"))

# Initialize GENESPACE
cat("\n=== Initializing GENESPACE ===\n")
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = "/usr/local/bin", 
  nCores = 20,
  verbose = TRUE
)

# Run GENESPACE
cat("\n=== Running GENESPACE ===\n")
out <- run_genespace(gpar, overwrite = TRUE)

# Query pangenome with TAIR10 as reference
cat("\n=== Querying pangenome ===\n")
pangenome <- query_pangenes(
  out, 
  bed = NULL, 
  refGenome = "TAIR10",
  transform = TRUE, 
  showArrayMem = TRUE, 
  showNSOrtho = TRUE,
  maxMem2Show = Inf
)

# Save pangenome object as RDS
output_file <- file.path(wd, "pangenome_matrix.rds")
cat("\n=== Saving pangenome matrix ===\n")
saveRDS(pangenome, file = output_file)
cat("Saved to:", output_file, "\n")

cat("\n=== GENESPACE analysis complete! ===\n")
