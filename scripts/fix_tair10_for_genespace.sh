#!/bin/bash

# Fix TAIR10 files to remove isoform suffixes (.1, .2, etc.)
# GENESPACE expects gene IDs without isoform notation

set -euo pipefail

GENESPACE_DIR="/data/users/afairman/euk_org_ann/results/genespace"
BED_DIR="${GENESPACE_DIR}/bed"
PEPTIDE_DIR="${GENESPACE_DIR}/peptide"

echo "=========================================="
echo "Fixing TAIR10 files for GENESPACE"
echo "=========================================="
echo ""

# Backup original files
echo "Creating backups..."
cp "${BED_DIR}/TAIR10.bed" "${BED_DIR}/TAIR10.bed.original"
cp "${PEPTIDE_DIR}/TAIR10.fa" "${PEPTIDE_DIR}/TAIR10.fa.original"
echo "  Backed up TAIR10.bed to TAIR10.bed.original"
echo "  Backed up TAIR10.fa to TAIR10.fa.original"
echo ""

# Fix BED file - remove isoform suffix from 4th column (gene ID)
echo "Fixing TAIR10.bed..."
awk 'BEGIN{OFS="\t"} {
  # Remove .1, .2, etc from gene ID (4th column)
  gsub(/\.[0-9]+$/, "", $4);
  print $0
}' "${BED_DIR}/TAIR10.bed.original" > "${BED_DIR}/TAIR10.bed"

BED_LINES=$(wc -l < "${BED_DIR}/TAIR10.bed")
echo "  Processed $BED_LINES genes"
echo ""

# Fix FASTA file - remove isoform suffix from headers and keep longest isoform per gene
echo "Fixing TAIR10.fa (keeping longest isoform per gene)..."
awk '
BEGIN { RS=">"; FS="\n" }
NR > 1 {
  # Get header and sequence
  header = $1;
  seq = "";
  for(i=2; i<=NF; i++) seq = seq $i;
  
  # Extract gene ID (remove .1, .2, etc)
  split(header, parts, /[ \t]/);
  gene_id = parts[1];
  gsub(/\.[0-9]+$/, "", gene_id);
  
  # Keep longest isoform per gene
  if(!(gene_id in best_len) || length(seq) > best_len[gene_id]) {
    best_len[gene_id] = length(seq);
    best_seq[gene_id] = seq;
  }
}
END {
  for(g in best_seq) {
    printf(">%s\n", g);
    seq = best_seq[g];
    # Wrap at 60 characters
    for(i=1; i<=length(seq); i+=60) {
      printf("%s\n", substr(seq, i, 60));
    }
  }
}
' "${PEPTIDE_DIR}/TAIR10.fa.original" > "${PEPTIDE_DIR}/TAIR10.fa"

FA_COUNT=$(grep -c "^>" "${PEPTIDE_DIR}/TAIR10.fa")
echo "  Processed $FA_COUNT genes"
echo ""

# Show sample of fixed files
echo "=========================================="
echo "Sample of fixed files:"
echo "=========================================="
echo ""
echo "TAIR10.bed (first 3 lines):"
head -3 "${BED_DIR}/TAIR10.bed"
echo ""
echo "TAIR10.fa (first header):"
head -1 "${PEPTIDE_DIR}/TAIR10.fa"
echo ""

echo "=========================================="
echo "TAIR10 files fixed!"
echo "=========================================="
echo "You can now run GENESPACE:"
echo "  sbatch scripts/11_run_genespace.sh"
echo ""
