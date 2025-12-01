#!/bin/bash
#SBATCH --job-name=genespace_refs
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/genespace_refs_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/genespace_refs_%j.err

# Process TAIR10 and additional Lian et al accessions for GENESPACE
# Output structure:
#   /data/users/afairman/euk_org_ann/results/genespace/
#   ├─ peptide/
#   │  ├─ TAIR10.fa
#   │  ├─ Altai-5.fa
#   │  ├─ Est-0.fa
#   │  └─ Taz-0.fa
#   └─ bed/
#      ├─ TAIR10.bed
#      ├─ Altai-5.bed
#      ├─ Est-0.bed
#      └─ Taz-0.bed

set -euo pipefail

# Output directories
OUTROOT="/data/users/afairman/euk_org_ann/results/genespace"
PEPTIDE_DIR="${OUTROOT}/peptide"
BED_DIR="${OUTROOT}/bed"
mkdir -p "$PEPTIDE_DIR" "$BED_DIR"

# Source directory for Lian et al data
LIAN_DIR="/data/courses/assembly-annotation-course/CDS_annotation/data/Lian_et_al"
TAIR_DIR="/data/courses/assembly-annotation-course/CDS_annotation"

echo "=========================================="
echo "Processing reference genomes for GENESPACE"
echo "=========================================="
echo ""

# Function to process a single accession
process_accession() {
  local ACCESSION=$1
  local GFF=$2
  local PROTEIN=$3
  
  echo "=== Processing ${ACCESSION} ==="
  
  if [ ! -f "$GFF" ]; then
    echo "ERROR: GFF not found at $GFF" >&2
    return 1
  fi
  if [ ! -f "$PROTEIN" ]; then
    echo "ERROR: Protein FASTA not found at $PROTEIN" >&2
    return 1
  fi
  
  TMPPREFIX="/tmp/${ACCESSION}.$$"
  
  # Extract gene lines and create BED
  echo "[1/2] Creating BED file..."
  grep -P "\tgene\t" "$GFF" > ${TMPPREFIX}.genes.gff3
  
  OUTBED="${BED_DIR}/${ACCESSION}.bed"
  awk -F"\t" 'BEGIN{OFS="\t"}
    {attr=$9; 
     if(match(attr,/ID=([^;]+)/,m)){ id=m[1] } 
     else if(match(attr,/Name=([^;]+)/,m2)){ id=m2[1] } 
     else { id="." }
     chr=$1; start=$4-1; end=$5; 
     print chr, start, end, id}
  ' ${TMPPREFIX}.genes.gff3 > "$OUTBED"
  
  echo "   Created: $OUTBED ($(wc -l < "$OUTBED") genes)"
  
  # Create gene list
  GENELIST="${TMPPREFIX}.genes.list"
  cut -f4 "$OUTBED" | sort -u > "$GENELIST"
  
  # Build longest-per-gene protein FASTA
  echo "[2/2] Creating longest-per-gene FASTA..."
  OUTFA="${PEPTIDE_DIR}/${ACCESSION}.fa"
  
  awk -v GENELIST="$GENELIST" -v OUTFA="$OUTFA" '
    BEGIN{ 
      while((getline < GENELIST) > 0){ 
        g=$1; 
        sub(/^[ \t]+|[ \t]+$/,"",g); 
        if(g!="") genes[g]=1 
      } 
      RS=">"; FS="\n" 
    }
    NR>1{
      header=$1; 
      seq=""; 
      for(i=2;i<=NF;i++) seq=seq $i;
      split(header,parts,/\s+/); 
      rawid=parts[1]; 
      id=rawid;
      # Extract gene ID (handle different isoform naming conventions)
      # For Lian et al: gene.1, gene.2, etc OR just gene ID
      if(index(id,".")){ 
        split(id,a,"."); 
        gid=a[1] 
      } else { 
        gid=id 
      }
      # Keep only genes in our list
      if(!(gid in genes)) next;
      # Keep longest isoform per gene
      if(!(gid in best_len) || length(seq) > best_len[gid]){ 
        best_len[gid]=length(seq); 
        best_seq[gid]=seq 
      }
    }
    END{
      for(g in best_seq){
        printf(">%s\n", g) >> OUTFA; 
        seq=best_seq[g];
        # Wrap sequences at 60 characters
        for(i=1;i<=length(seq); i+=60) {
          printf("%s\n", substr(seq,i,60)) >> OUTFA
        }
      }
      printf("   Wrote %d sequences to %s\n", length(best_seq), OUTFA)
    }
  ' "$PROTEIN"
  
  # Clean up temp files
  rm -f ${TMPPREFIX}.genes.gff3 ${TMPPREFIX}.genes.list || true
  
  echo "   Done: ${ACCESSION}"
  echo ""
}

# Process TAIR10 reference
echo "=== Processing TAIR10 ==="
if [ -f "${TAIR_DIR}/TAIR10.fa" ] && [ -f "${TAIR_DIR}/TAIR10.bed" ]; then
  cp "${TAIR_DIR}/TAIR10.fa" "$PEPTIDE_DIR/"
  cp "${TAIR_DIR}/TAIR10.bed" "$BED_DIR/"
  echo "   Copied TAIR10.fa to peptide/"
  echo "   Copied TAIR10.bed to bed/"
  echo "   Done: TAIR10"
  echo ""
else
  echo "WARNING: TAIR10 files not found in ${TAIR_DIR}"
  echo ""
fi

# Process Altai-5
process_accession "Altai-5" \
  "${LIAN_DIR}/gene_gff/selected/Altai-5.EVM.v3.5.ann.protein_coding_genes.gff" \
  "${LIAN_DIR}/protein/selected/Altai-5.protein.faa"

# Process Est-0
process_accession "Est-0" \
  "${LIAN_DIR}/gene_gff/selected/Est-0.EVM.v3.5.ann.protein_coding_genes.gff" \
  "${LIAN_DIR}/protein/selected/Est-0.protein.faa"

# Process Taz-0
process_accession "Taz-0" \
  "${LIAN_DIR}/gene_gff/selected/Taz-0.EVM.v3.5.ann.protein_coding_genes.gff" \
  "${LIAN_DIR}/protein/selected/Taz-0.protein.faa"

echo "=========================================="
echo "Summary - Directory structure"
echo "=========================================="
echo ""
echo "peptide/:"
ls -lh "$PEPTIDE_DIR" | tail -n +2 | awk '{printf "  %-20s %8s\n", $9, $5}'
echo ""
echo "bed/:"
ls -lh "$BED_DIR" | tail -n +2 | awk '{printf "  %-20s %8s\n", $9, $5}'
echo ""
echo "All files ready in: $OUTROOT"
echo "Done!"
