#!/bin/bash

# Rename genespace files to remove hyphens (not allowed by GENESPACE)
# GENESPACE only allows alphanumeric characters, underscores, and periods

set -euo pipefail

GENESPACE_DIR="/data/users/afairman/euk_org_ann/results/genespace"
BED_DIR="${GENESPACE_DIR}/bed"
PEPTIDE_DIR="${GENESPACE_DIR}/peptide"

echo "=========================================="
echo "Renaming files to fix GENESPACE genome IDs"
echo "=========================================="
echo "GENESPACE requires alphanumeric + _ and . only"
echo "Converting hyphens to underscores"
echo ""

# Array of files to rename (hyphen -> underscore)
declare -A RENAME_MAP=(
  ["Altai-5"]="Altai_5"
  ["Est-0"]="Est_0"
  ["Taz-0"]="Taz_0"
)

for old_name in "${!RENAME_MAP[@]}"; do
  new_name="${RENAME_MAP[$old_name]}"
  
  echo "Renaming: $old_name -> $new_name"
  
  # Rename BED file
  if [ -f "${BED_DIR}/${old_name}.bed" ]; then
    mv "${BED_DIR}/${old_name}.bed" "${BED_DIR}/${new_name}.bed"
    echo "  ✓ ${old_name}.bed -> ${new_name}.bed"
  fi
  
  # Rename FASTA file
  if [ -f "${PEPTIDE_DIR}/${old_name}.fa" ]; then
    mv "${PEPTIDE_DIR}/${old_name}.fa" "${PEPTIDE_DIR}/${new_name}.fa"
    echo "  ✓ ${old_name}.fa -> ${new_name}.fa"
  fi
  
  echo ""
done

echo "=========================================="
echo "Final file listing:"
echo "=========================================="
echo ""
echo "BED files:"
ls -1 "${BED_DIR}"/*.bed | xargs -n1 basename
echo ""
echo "FASTA files:"
ls -1 "${PEPTIDE_DIR}"/*.fa | xargs -n1 basename
echo ""
echo "=========================================="
echo "Files renamed successfully!"
echo "=========================================="
echo "You can now run GENESPACE:"
echo "  sbatch scripts/11_run_genespace.sh"
echo ""
