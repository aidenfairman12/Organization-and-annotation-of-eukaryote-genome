#!/bin/bash

#SBATCH --job-name=03_TE_circos_plots
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/03_TE_circos_plots_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/03_TE_circos_plots_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=pibu_el8

# Generate TE distribution circos plots
# This script runs the R script to create circular plots showing
# the distribution of transposable elements across the genome

# Load required modules
module load R/4.3.1

# Set working directory
SCRIPT_DIR="/data/users/afairman/euk_org_ann/scripts"
WORK_DIR="/data/users/afairman/euk_org_ann/results/edta_output"
PLOTS_DIR="/data/users/afairman/euk_org_ann/plots"

# Create plots directory if it doesn't exist
mkdir -p "$PLOTS_DIR"

echo "========================================="
echo "TE Circos Plot Generation"
echo "========================================="
echo "Start time: $(date)"
echo "Script directory: $SCRIPT_DIR"
echo "Working directory: $WORK_DIR"
echo "Output directory: $PLOTS_DIR"
echo ""

# Check if required input files exist
if [ ! -f "$WORK_DIR/assembly.fasta.mod.EDTA.TEanno.gff3" ]; then
    echo "ERROR: TE annotation GFF3 file not found!"
    echo "Expected: $WORK_DIR/assembly.fasta.mod.EDTA.TEanno.gff3"
    exit 1
fi

if [ ! -f "$WORK_DIR/assembly.fasta.fai" ]; then
    echo "ERROR: FAI index file not found!"
    echo "Expected: $WORK_DIR/assembly.fasta.fai"
    exit 1
fi

echo "Input files verified ✓"
echo ""

# Run the R script
echo "Running R script..."
Rscript "$SCRIPT_DIR/03_TE_circos_plots.R"

EXIT_CODE=$?

echo ""
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ TE circos plots generated successfully!"
    echo ""
    echo "Output files:"
    ls -lh "$PLOTS_DIR"/03_TE_*.pdf 2>/dev/null || echo "  (No PDF files found - check logs for errors)"
else
    echo "✗ Script failed with exit code: $EXIT_CODE"
    echo "Check error log for details"
fi
echo "========================================="
echo "End time: $(date)"

exit $EXIT_CODE
