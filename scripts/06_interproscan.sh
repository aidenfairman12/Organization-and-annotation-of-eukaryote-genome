#!/bin/bash
#SBATCH --time=1-0
#SBATCH --mem=128G
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=06_interproscan
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/06_interproscan_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/06_interproscan_%j.err

set -euo pipefail

# Paths (edit if your layout differs)
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/afairman/euk_org_ann/results/Maker/final"
CONTAINER="${COURSEDIR}/containers/interproscan_latest.sif"

# Ensure SCRATCH has a sensible default if not set in environment
: "${SCRATCH:=/scratch/$USER}"

# Usage: sbatch 06_interproscan.sh [protein_basename]
# protein_basename should be the base name used earlier (e.g. "assembly.all.maker.proteins.fasta").
# The script will use ${protein}.renamed.fasta as input.
protein="assembly.all.maker.proteins.fasta"

echo "Running InterProScan on \"${protein}.renamed.fasta\" in ${WORKDIR}"

cd "${WORKDIR}"

apptainer exec \
  --bind "${COURSEDIR}/data/interproscan-5.70-102.0/data":/opt/interproscan/data \
  --bind "${WORKDIR}" \
  --bind "${COURSEDIR}" \
  --bind "${SCRATCH}":/temp \
  "${CONTAINER}" \
  /opt/interproscan/interproscan.sh \
  -appl pfam --disable-precalc -f TSV \
  --goterms --iprlookup --seqtype p \
  -i "${protein}.renamed.fasta" -o output.iprscan

echo "InterProScan finished. Output: ${WORKDIR}/output.iprscan"
