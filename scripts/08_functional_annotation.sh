#!/bin/bash
# Functional annotation pipeline: BLASTP against UniProt and TAIR, keep best hits,
# then run MAKER functional annotation helpers to add names/descriptions to fasta and GFF.
# Usage: ./08_functional_annotation.sh [threads]

set -euo pipefail

THREADS=${1:-10}
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/afairman/euk_org_ann/results/Maker/final"
MAKERBIN="/data/courses/assembly-annotation-course/CDS_annotation/softwares/Maker_v3.01.03/src/bin"

# Inputs
PROTEIN_FASTA="${WORKDIR}/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"
TRANSCRIPT_FASTA="${WORKDIR}/assembly.all.maker.transcripts.fasta.renamed.filtered.fasta"

# UniProt
UNIPROT_FASTA="${COURSEDIR}/data/uniprot/uniprot_viridiplantae_reviewed.fa"
UNIPROT_DB="${UNIPROT_FASTA}"   # already formatted blast DB (user said makeblastdb done)
UNIPROT_BLAST_OUT="${WORKDIR}/blastp_uniprot.out"
UNIPROT_BEST="${UNIPROT_BLAST_OUT}.besthits"

# TAIR10
TAIR_DB="${COURSEDIR}/data/TAIR10_pep_20110103_representative_gene_model"
TAIR_BLAST_OUT="${WORKDIR}/blastp_tair.out"
TAIR_BEST="${TAIR_BLAST_OUT}.besthits"

# GFF to annotate (try common locations)
GFF_CANDIDATES=("${WORKDIR}/filtered.genes.renamed.gff3" "${WORKDIR}/filtered.genes.renamed.gff" "${WORKDIR}/assembly.all.maker.noseq.gff")
GFF_INPUT=""
for g in "${GFF_CANDIDATES[@]}"; do
  if [ -f "$g" ]; then
    GFF_INPUT="$g"
    break
  fi
done
if [ -z "$GFF_INPUT" ]; then
  echo "WARNING: could not find filtered genes GFF; you may need to edit the script to point to your GFF."
  echo "Attempting to continue without GFF-based maker_functional_gff step."
fi

# make sure MAKERBIN exists
if [ ! -x "${MAKERBIN}/maker_functional_fasta" ]; then
  echo "ERROR: maker_functional_fasta not found/executable at ${MAKERBIN}" >&2
  exit 1
fi

# Load BLAST
module load BLAST+/2.15.0-gompi-2021a || true

cd "${WORKDIR}"

# 1) BLAST against UniProt
echo "Running blastp against UniProt (threads=${THREADS})"
blastp -query "${PROTEIN_FASTA}" -db "${UNIPROT_DB}" -num_threads ${THREADS} -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out "${UNIPROT_BLAST_OUT}"

# sort to keep best hit per query (keep lowest evalue in column 12; then unique by query)
# columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# sort by query, then by evalue (col 11 or 12 depending on outfmt). Here we use evalue (field 11) and bitscore (12) as tie-breaker.
# Use generic safe sort: primary key qseqid, secondary numeric evalue ascending, tertiary bitscore descending.
awk 'BEGIN{OFS="\t"} {print $0}' "${UNIPROT_BLAST_OUT}" | sort -k1,1 -k11,11g -k12,12gr | sort -u -k1,1 --merge > "${UNIPROT_BEST}"

# 2) run maker functional fasta and gff annotation using UniProt results
cp "${PROTEIN_FASTA}" "${PROTEIN_FASTA}.Uniprot"
if [ -n "$GFF_INPUT" ]; then
  cp "$GFF_INPUT" "${GFF_INPUT}.Uniprot.gff3"
fi

echo "Running maker_functional_fasta with UniProt besthits"
"${MAKERBIN}/maker_functional_fasta" "${UNIPROT_FASTA}" "${UNIPROT_BEST}" "${PROTEIN_FASTA}" > "${PROTEIN_FASTA}.Uniprot"

if [ -n "$GFF_INPUT" ]; then
  echo "Running maker_functional_gff with UniProt besthits"
  "${MAKERBIN}/maker_functional_gff" "${UNIPROT_FASTA}" "${UNIPROT_BEST}" "$GFF_INPUT" > "${GFF_INPUT}.Uniprot.gff3"
fi

# 3) BLAST against TAIR10 (optional but requested)
echo "Running blastp against TAIR10 (threads=${THREADS})"
blastp -query "${PROTEIN_FASTA}" -db "${TAIR_DB}" -num_threads ${THREADS} -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out "${TAIR_BLAST_OUT}"

awk 'BEGIN{OFS="\t"} {print $0}' "${TAIR_BLAST_OUT}" | sort -k1,1 -k11,11g -k12,12gr | sort -u -k1,1 --merge > "${TAIR_BEST}"

echo "Functional annotation steps complete. Outputs:"
ls -lh "${UNIPROT_BLAST_OUT}" "${UNIPROT_BEST}" "${PROTEIN_FASTA}.Uniprot" "${TAIR_BLAST_OUT}" "${TAIR_BEST}" || true
if [ -n "$GFF_INPUT" ]; then ls -lh "${GFF_INPUT}.Uniprot.gff3" || true; fi

exit 0
