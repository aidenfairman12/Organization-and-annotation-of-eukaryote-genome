#!/usr/bin/env Rscript

# Script to count core genes and genes in accession from genespace pangenome matrix
# Usage: Rscript 20_count_genes.R [path_to_pangenome_matrix.rds]
# 
# Counts:
# - Core genes: gene families present in all accessions
# - Genes per accession: total number of genes in each accession

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set default path or use provided argument
if (length(args) == 0) {
  pangenome_path <- "/data/users/afairman/euk_org_ann/results/genespace/pangenome_matrix.rds"
} else {
  pangenome_path <- args[1]
}

cat("Loading pangenome matrix from:", pangenome_path, "\n")

# Check if file exists
if (!file.exists(pangenome_path)) {
  stop("File not found: ", pangenome_path)
}

# Load pangenome object (data.table where rows are gene families, 
# columns are accessions with character vectors of genes)
pangenome <- readRDS(pangenome_path)

cat("\n=== GENESPACE GENE COUNTS ===\n\n")

# Identify accession columns (exclude metadata columns)
metadata_cols <- c("pgID", "interpChr", "interpOrd", "og", "repGene", "genome", "chr", "start", "end")
accession_cols <- setdiff(colnames(pangenome), metadata_cols)

cat("Accessions found:", paste(accession_cols, collapse = ", "), "\n\n")

# Count core genes (gene families present in all accessions)
# A family is in an accession if the character vector is not empty
core_genes <- 0
for (i in 1:nrow(pangenome)) {
  present_in_all <- TRUE
  for (acc in accession_cols) {
    genes <- pangenome[[i, acc]]
    # Check if genes is empty (length 0, NA, or all empty strings)
    is_empty <- length(genes) == 0 || is.na(genes[1]) || all(genes == "")
    if (is_empty) {
      present_in_all <- FALSE
      break
    }
  }
  if (present_in_all) {
    core_genes <- core_genes + 1
  }
}

cat("Core genes (families present in all accessions):", core_genes, "\n\n")

# Count total genes per accession and accession-specific genes
cat("Genes per accession:\n")
genes_per_accession <- list()
accession_specific_genes <- list()

for (acc in accession_cols) {
  total_genes <- 0
  specific_genes <- 0
  
  for (i in 1:nrow(pangenome)) {
    genes <- pangenome[[i, acc]]
    # Check if genes is empty (length 0, NA, or all empty strings)
    is_empty <- length(genes) == 0 || is.na(genes[1]) || all(genes == "")
    
    if (!is_empty) {
      # Count non-empty genes
      gene_count <- length(genes)
      total_genes <- total_genes + gene_count
      
      # Check if this gene family is only in this accession
      only_in_this_acc <- TRUE
      for (other_acc in accession_cols) {
        if (other_acc != acc) {
          other_genes <- pangenome[[i, other_acc]]
          other_is_empty <- length(other_genes) == 0 || is.na(other_genes[1]) || all(other_genes == "")
          if (!other_is_empty) {
            only_in_this_acc <- FALSE
            break
          }
        }
      }
      
      if (only_in_this_acc) {
        specific_genes <- specific_genes + gene_count
      }
    }
  }
  
  genes_per_accession[[acc]] <- total_genes
  accession_specific_genes[[acc]] <- specific_genes
  cat(sprintf("  %s: %d genes (%d accession-specific)\n", acc, total_genes, specific_genes))
}

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total gene families:", nrow(pangenome), "\n")
cat("Number of accessions:", length(accession_cols), "\n")
cat("Core gene families:", core_genes, "\n")
cat("Core genes (%):", sprintf("%.2f%%", (core_genes / nrow(pangenome)) * 100), "\n")
cat("Mean genes per accession:", sprintf("%.0f", mean(unlist(genes_per_accession))), "\n")
cat("Mean accession-specific genes:", sprintf("%.0f", mean(unlist(accession_specific_genes))), "\n")
cat("Min genes per accession:", min(unlist(genes_per_accession)), "\n")
cat("Max genes per accession:", max(unlist(genes_per_accession)), "\n")
cat("Total accession-specific genes:", sum(unlist(accession_specific_genes)), "\n")

# Save results to output file
output_file <- "/data/users/afairman/euk_org_ann/results/genespace/gene_counts.txt"
cat("\n=== Saving results to:", output_file, "===\n\n")

sink(output_file)
cat("=== GENESPACE GENE COUNTS ===\n\n")
cat("Accessions:", paste(accession_cols, collapse = ", "), "\n\n")
cat("Core genes (families present in all accessions):", core_genes, "\n\n")
cat("Genes per accession:\n")
for (acc in accession_cols) {
  cat(sprintf("  %s: %d genes (%d accession-specific)\n", acc, genes_per_accession[[acc]], accession_specific_genes[[acc]]))
}
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total gene families:", nrow(pangenome), "\n")
cat("Number of accessions:", length(accession_cols), "\n")
cat("Core gene families:", core_genes, "\n")
cat("Core genes (%):", sprintf("%.2f%%", (core_genes / nrow(pangenome)) * 100), "\n")
cat("Mean genes per accession:", sprintf("%.0f", mean(unlist(genes_per_accession))), "\n")
cat("Mean accession-specific genes:", sprintf("%.0f", mean(unlist(accession_specific_genes))), "\n")
cat("Min genes per accession:", min(unlist(genes_per_accession)), "\n")
cat("Max genes per accession:", max(unlist(genes_per_accession)), "\n")
cat("Total accession-specific genes:", sum(unlist(accession_specific_genes)), "\n")
sink()

cat("Results saved successfully!\n")
