#!/bin/bash
#SBATCH --job-name=genespace_prep
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/genespace_prep_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/genespace_prep_%j.err

# Create BED and FASTA files for GENESPACE from filtered GFF3 and protein FASTA
# Usage: sbatch 09_make_genespace_bed_and_fasta.sh <Accession>
# Output structure:
#   /data/users/afairman/euk_org_ann/results/genespace/
#   ├─ peptide/
#   │  ├─ Accession.fa
#   │  └─ TAIR10.fa
#   └─ bed/
#      ├─ Accession.bed
#      └─ TAIR10.bed

set -euo pipefail

Accession="Co4"
WORKDIR="/data/users/afairman/euk_org_ann/results/Maker/final"
GFF="${WORKDIR}/filtered.genes.renamed.gff3"
PROT_FASTA="${WORKDIR}/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"

# Central output root for all genespace products
OUTROOT="/data/users/afairman/euk_org_ann/results/genespace"
PEPTIDE_DIR="${OUTROOT}/peptide"
BED_DIR="${OUTROOT}/bed"

# Create output directories
mkdir -p "$PEPTIDE_DIR" "$BED_DIR" "$OUTROOT"

# Validate input files
if [ ! -f "$GFF" ]; then
  echo "ERROR: GFF not found at $GFF" >&2
  exit 1
fi
if [ ! -f "$PROT_FASTA" ]; then
  echo "ERROR: protein FASTA not found at $PROT_FASTA" >&2
  exit 1
fi

echo "=== GENESPACE file preparation for ${Accession} ==="
echo "Working directory: $WORKDIR"
echo "Output root: $OUTROOT"

# Use a unique temp prefix in /tmp
TMPPREFIX="/tmp/${Accession}.$$"

# Extract gene lines and format to BED (0-based start): chr\tstart\tend\tgene_name
# GFF format columns: seqid, source, type, start, end, score, strand, phase, attributes
echo "[1/3] Extracting gene lines from GFF and creating BED file..."
grep -P "\tgene\t" "$GFF" > ${TMPPREFIX}.genes.gff3

OUTBED="${BED_DIR}/${Accession}.bed"
awk -F"\t" 'BEGIN{OFS="\t"}
  {attr=$9; 
   if(match(attr,/ID=([^;]+)/,m)){ id=m[1] } 
   else if(match(attr,/Name=([^;]+)/,m2)){ id=m2[1] } 
   else { id="." }
   chr=$1; start=$4-1; end=$5; 
   print chr, start, end, id}
' ${TMPPREFIX}.genes.gff3 > "$OUTBED"

echo "   Created BED: $OUTBED"

# Create gene list for filtering proteins
echo "[2/3] Building gene list from BED..."
GENELIST="${TMPPREFIX}.genes.list"
cut -f4 "$OUTBED" | sort -u > "$GENELIST"

# Build longest-per-gene protein FASTA
echo "[3/3] Creating longest-per-gene FASTA..."
OUTFA="${PEPTIDE_DIR}/${Accession}.fa"

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
    # Extract gene ID (remove isoform suffix if present)
    if(index(id,"-R")){ 
      split(id,a,"-R"); 
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
' "$PROT_FASTA"

# Copy TAIR10 reference files (if present)
echo "[4/4] Copying TAIR10 reference files..."
TAIR_SRC_DIR="/data/courses/assembly-annotation-course/CDS_annotation"

if [ -f "${TAIR_SRC_DIR}/TAIR10.fa" ]; then
  cp "${TAIR_SRC_DIR}/TAIR10.fa" "$PEPTIDE_DIR/"
  echo "   Copied TAIR10.fa to peptide/"
fi

if [ -f "${TAIR_SRC_DIR}/TAIR10.bed" ]; then
  cp "${TAIR_SRC_DIR}/TAIR10.bed" "$BED_DIR/"
  echo "   Copied TAIR10.bed to bed/"
fi

# Clean up temp files
rm -f ${TMPPREFIX}.genes.gff3 ${TMPPREFIX}.genes.list || true

echo ""
echo "=== Directory structure ==="
echo "peptide/:"
ls -lh "$PEPTIDE_DIR" | tail -n +2 | awk '{printf "  %s  %s\n", $9, $5}'
echo ""
echo "bed/:"
ls -lh "$BED_DIR" | tail -n +2 | awk '{printf "  %s  %s\n", $9, $5}'

echo ""
echo "Done: files for GENESPACE ready in $OUTROOT"
