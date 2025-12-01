#!/bin/bash

#SBATCH --job-name=02_TEsorter
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/02_TEsorter_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/02_TEsorter_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=pibu_el8

# TEsorter clade classification for intact LTRs
# Inputs from EDTA raw outputs

WORKDIR="/data/users/afairman/euk_org_ann/results/edta_output/assembly.fasta.mod.EDTA.raw"
CONTAINER="/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif"
INPUT_FA="${WORKDIR}/assembly.fasta.mod.LTR.intact.raw.fa"
INPUT_GFF="${WORKDIR}/assembly.fasta.mod.LTR.intact.raw.gff3"
DB="rexdb-plant"

mkdir -p /data/users/afairman/euk_org_ann/logs

# sanity checks
if [ ! -f "$INPUT_FA" ]; then
    echo "ERROR: Input LTR fasta not found: $INPUT_FA" >&2
    exit 2
fi

if [ ! -f "$INPUT_GFF" ]; then
    echo "WARNING: GFF3 not found: $INPUT_GFF - proceeding without GFF" >&2
fi

# Change to workdir so outputs are written next to inputs
cd "$WORKDIR"

echo "Running TEsorter"
echo "Workdir: $WORKDIR"
echo "Container: $CONTAINER"
echo "Input fasta: $INPUT_FA"
echo "DB: $DB"
echo "Start: $(date)"

apptainer exec --bind /data/users/afairman:/data/users/afairman --bind /data/courses:/data/courses "$CONTAINER" \
    TEsorter \
    "${INPUT_FA}" -db $DB

EXIT_CODE=$?

echo "TEsorter exited with code: $EXIT_CODE"
echo "Finished: $(date)"

exit $EXIT_CODE
