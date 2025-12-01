#!/bin/bash

#SBATCH --job-name=01_edta_annotation
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/01_edta_annotation_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/01_edta_annotation_%j.err
#SBATCH --time=50:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --partition=pibu_el8

# TE annotation using EDTA on Flye assembly
# Author: Generated script
# Date: $(date)

# Set variables
ASSEMBLY="/data/users/afairman/assembly_annotation_course/flye_results/assembly.fasta"
CONTAINER="/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif"
CDS_FILE="/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_cds_20110103_representative_gene_model_updated"
OUTPUT_DIR="/data/users/afairman/euk_org_ann/results/edta_output"
THREADS=20


# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Change to output directory
cd "$OUTPUT_DIR"

echo "Starting EDTA TE annotation..."
echo "Assembly: $ASSEMBLY"
echo "Container: $CONTAINER"
echo "CDS file: $CDS_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Start time: $(date)"

# Run EDTA using Apptainer with bind mounts
apptainer exec \
    --bind /data/users/afairman:/data/users/afairman \
    --bind /data/courses:/data/courses \
    "$CONTAINER" EDTA.pl \
    --genome "$ASSEMBLY" \
    --species others \
    --step all \
    --sensitive 1 \
    --cds "$CDS_FILE" \
    --anno 1 \
    --threads "$THREADS"

echo "EDTA annotation completed at: $(date)"

# Check if output files were generated
if [ -f "*.EDTA.TEanno.gff3" ]; then
    echo "SUCCESS: TE annotation files generated"
    ls -la *.EDTA.*
else
    echo "WARNING: Expected output files not found"
    echo "Contents of output directory:"
    ls -la
fi

echo "Script finished at: $(date)"