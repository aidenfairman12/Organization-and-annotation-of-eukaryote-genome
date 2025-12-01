#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=60G
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=07_BUSCO
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/07_BUSCO_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/07_BUSCO_%j.err

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/afairman/euk_org_ann/results/Maker/final"
# base names (without .renamed.filtered.fasta suffix)
protein_base="assembly.all.maker.proteins.fasta.renamed.filtered"
transcript_base="assembly.all.maker.transcripts.fasta.renamed.filtered"
gff="assembly.all.maker.noseq.gff"

# full input/output paths
protein_in="${WORKDIR}/${protein_base}.fasta"
transcript_in="${WORKDIR}/${transcript_base}.fasta"
protein_long="${WORKDIR}/${protein_base}.longest.fasta"
transcript_long="${WORKDIR}/${transcript_base}.longest.fasta"

module load BUSCO/5.4.2-foss-2021a

cd "${WORKDIR}"

echo "[BUSCO] extracting longest isoform per gene for proteins and transcripts"
# Use the helper script to extract the longest sequence per gene. If the script
# is not present, fall back to a simple awk extract (shouldn't happen in this repo).
EXTRACT_SCRIPT="/data/users/afairman/euk_org_ann/scripts/07_extract_longest_per_gene.py"
if [ -x "${EXTRACT_SCRIPT}" ]; then
	echo "Running python longest-extract script"
	"${EXTRACT_SCRIPT}" "${protein_in}" "${protein_long}"
	"${EXTRACT_SCRIPT}" "${transcript_in}" "${transcript_long}"
else
	echo "Warning: extract script not found at ${EXTRACT_SCRIPT}; trying inline awk fallback"
	awk 'BEGIN{ while(getline < "list.txt"){ gsub(/^>/,"",$0); ids[$1]=1 } } /^>/{ split($0,h,"[ \t]"); id=substr(h[1],2); keep=(id in ids) } keep{ print }' "${protein_in}" > "${protein_long}"
	awk 'BEGIN{ while(getline < "list.txt"){ gsub(/^>/,"",$0); ids[$1]=1 } } /^>/{ split($0,h,"[ \t]"); id=substr(h[1],2); keep=(id in ids) } keep{ print }' "${transcript_in}" > "${transcript_long}"
fi


echo "starting busco transcripts on longest-per-gene set"
busco -i "${transcript_long}" -l brassicales_odb10 -o busco_transcripts -m transcriptome -f