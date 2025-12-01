#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=051_map_ids
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/051_map_ids_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/051_map_ids_%j.err

set -euo pipefail

# Default paths - adjust if needed
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
MAKERBIN="${COURSEDIR}/softwares/Maker_v3.01.03/src/bin"
WORKDIR="/data/users/afairman/euk_org_ann/results/Maker"

# Usage: 051_map_ids.sh <prefix> <gff> <protein_fasta> <transcript_fasta>
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <prefix> <gff> <protein_fasta> <transcript_fasta>"
  echo "Example: $0 sample1 assembly.all.maker.gff proteins.fasta transcripts.fasta"
  exit 1
fi

prefix="$1"
gff="$2"
protein="$3"
transcript="$4"

cd "${WORKDIR}"

echo "Starting ID mapping: prefix=${prefix}, gff=${gff}, protein=${protein}, transcript=${transcript}"

echo "Checking MAKER bin: ${MAKERBIN}"
if [ ! -x "${MAKERBIN}/maker_map_ids" ]; then
  echo "ERROR: maker_map_ids not found or not executable in ${MAKERBIN}" >&2
  exit 2
fi

# 1) generate id.map
echo "Running maker_map_ids -> id.map"
"${MAKERBIN}/maker_map_ids" --prefix "${prefix}" --justify 7 "${gff}.renamed.gff" > id.map
rc=$?
if [ $rc -ne 0 ]; then
  echo "maker_map_ids failed with exit code $rc" >&2
  exit $rc
fi

# 2) map GFF ids
echo "Mapping GFF ids"
"${MAKERBIN}/map_gff_ids" id.map "${gff}.renamed.gff"
rc=$?
if [ $rc -ne 0 ]; then
  echo "map_gff_ids failed with exit code $rc" >&2
  exit $rc
fi

# 3) map protein fasta ids
echo "Mapping protein FASTA ids"
"${MAKERBIN}/map_fasta_ids" id.map "${protein}.renamed.fasta"
rc=$?
if [ $rc -ne 0 ]; then
  echo "map_fasta_ids (protein) failed with exit code $rc" >&2
  exit $rc
fi

# 4) map transcript fasta ids
echo "Mapping transcript FASTA ids"
"${MAKERBIN}/map_fasta_ids" id.map "${transcript}.renamed.fasta"
rc=$?
if [ $rc -ne 0 ]; then
  echo "map_fasta_ids (transcript) failed with exit code $rc" >&2
  exit $rc
fi

echo "ID mapping completed successfully. Outputs are in ${WORKDIR}: id.map and the renamed files."
