#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=08_functional
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/08_functional_annotation_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/08_functional_annotation_%j.err

set -euo pipefail

# Threads will match SLURM allocation
THREADS=${SLURM_CPUS_PER_TASK:-16}

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/afairman/euk_org_ann/results/Maker/final"
MAKERBIN="/data/courses/assembly-annotation-course/CDS_annotation/softwares/Maker_v3.01.03/src/bin"

# Inputs (filtered maker protein FASTA produced earlier)
PROTEIN_FASTA="${WORKDIR}/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"

# UniProt
UNIPROT_FASTA="${COURSEDIR}/data/uniprot/uniprot_viridiplantae_reviewed.fa"
UNIPROT_DB="${UNIPROT_FASTA}"
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
  echo "WARNING: could not find filtered genes GFF; maker_functional_gff step will be skipped unless you provide a GFF."
fi

cd "${WORKDIR}"

module load BLAST+/2.15.0-gompi-2021a || true

echo "[functional] Running BLASTP against UniProt (threads=${THREADS})"
blastp -query "${PROTEIN_FASTA}" -db "${UNIPROT_DB}" -num_threads ${THREADS} -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out "${UNIPROT_BLAST_OUT}"

# Keep best hit per query (evalue asc, bitscore desc)
awk 'BEGIN{OFS="\t"} {print $0}' "${UNIPROT_BLAST_OUT}" | sort -k1,1 -k11,11g -k12,12gr | sort -u -k1,1 --merge > "${UNIPROT_BEST}"

echo "[functional] Running maker_functional_fasta with UniProt besthits"
cp "${PROTEIN_FASTA}" "${PROTEIN_FASTA}.Uniprot"
"${MAKERBIN}/maker_functional_fasta" "${UNIPROT_FASTA}" "${UNIPROT_BEST}" "${PROTEIN_FASTA}" > "${PROTEIN_FASTA}.Uniprot"

if [ -n "$GFF_INPUT" ]; then
  echo "[functional] Running maker_functional_gff with UniProt besthits"
  cp "$GFF_INPUT" "${GFF_INPUT}.Uniprot.gff3"
  "${MAKERBIN}/maker_functional_gff" "${UNIPROT_FASTA}" "${UNIPROT_BEST}" "$GFF_INPUT" > "${GFF_INPUT}.Uniprot.gff3"
fi

echo "[functional] Running BLASTP against TAIR10 (threads=${THREADS})"
blastp -query "${PROTEIN_FASTA}" -db "${TAIR_DB}" -num_threads ${THREADS} -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out "${TAIR_BLAST_OUT}"
awk 'BEGIN{OFS="\t"} {print $0}' "${TAIR_BLAST_OUT}" | sort -k1,1 -k11,11g -k12,12gr | sort -u -k1,1 --merge > "${TAIR_BEST}"

echo "[functional] Done. Outputs in ${WORKDIR}:"
ls -lh "${UNIPROT_BLAST_OUT}" "${UNIPROT_BEST}" "${PROTEIN_FASTA}.Uniprot" "${TAIR_BLAST_OUT}" "${TAIR_BEST}" || true
if [ -n "$GFF_INPUT" ]; then ls -lh "${GFF_INPUT}.Uniprot.gff3" || true; fi

exit 0
