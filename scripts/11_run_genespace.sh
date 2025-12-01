#!/bin/bash
#SBATCH --job-name=genespace
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/genespace_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/genespace_%j.err

# GENESPACE analysis using Apptainer container
# This script runs synteny analysis and creates a pangenome matrix

set -euo pipefail

# Directories
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/afairman/euk_org_ann"
SCRATCH="/scratch/afairman"
GENESPACE_WD="${WORKDIR}/results/genespace"

# Create scratch directory if needed
mkdir -p "$SCRATCH"

echo "=========================================="
echo "GENESPACE Analysis"
echo "=========================================="
echo "Course directory: $COURSEDIR"
echo "Working directory: $WORKDIR"
echo "GENESPACE data: $GENESPACE_WD"
echo "Scratch: $SCRATCH"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "=========================================="
echo ""

# Change to working directory
cd "$WORKDIR"

# Run GENESPACE using Apptainer container
echo "Running GENESPACE in container..."
apptainer exec \
  --bind "$COURSEDIR" \
  --bind "$WORKDIR" \
  --bind "$SCRATCH:/temp" \
  "$COURSEDIR/containers/genespace_latest.sif" \
  Rscript scripts/11_run_genespace.R "$GENESPACE_WD"

echo ""
echo "=========================================="
echo "GENESPACE analysis complete!"
echo "=========================================="
echo "Output saved to: ${GENESPACE_WD}/pangenome_matrix.rds"
