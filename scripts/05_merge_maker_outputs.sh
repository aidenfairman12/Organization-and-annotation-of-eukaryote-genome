#!/bin/bash
#SBATCH --job-name=maker_merge
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/maker_merge_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/maker_merge_%j.err
#SBATCH --partition=pibu_el8
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Directories
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/afairman/euk_org_ann/results/Maker"
MAKERBIN="${COURSEDIR}/softwares/Maker_v3.01.03/src/bin"

# Move to working directory
cd "${WORKDIR}" || { echo "ERROR: cannot cd ${WORKDIR}" >&2; exit 1; }

# Check MAKER bin
if [ ! -x "${MAKERBIN}/gff3_merge" ] || [ ! -x "${MAKERBIN}/fasta_merge" ]; then
  echo "ERROR: MAKER merge binaries not found or not executable at ${MAKERBIN}" >&2
  ls -l "${MAKERBIN}" || true
  exit 2
fi

# Log file for this run
RUNLOG="merge_run.log"
echo "Merge run started: $(date -u)" > "${RUNLOG}"

# 1) gff3_merge (with sequences)
echo "Running: gff3_merge -s -d assembly.maker.output/assembly_master_datastore_index.log -> assembly.all.maker.gff" | tee -a "${RUNLOG}"
/usr/bin/time -v "${MAKERBIN}/gff3_merge" -s -d assembly.maker.output/assembly_master_datastore_index.log > assembly.all.maker.gff 2> gff3_merge_s.log
echo "gff3_merge (with seq) exit:$?" | tee -a "${RUNLOG}"

# 2) gff3_merge (no sequence)
echo "Running: gff3_merge -n -s -d assembly.maker.output/assembly_master_datastore_index.log -> assembly.all.maker.noseq.gff" | tee -a "${RUNLOG}"
/usr/bin/time -v "${MAKERBIN}/gff3_merge" -n -s -d assembly.maker.output/assembly_master_datastore_index.log > assembly.all.maker.noseq.gff 2> gff3_merge_ns.log
echo "gff3_merge (no seq) exit:$?" | tee -a "${RUNLOG}"

# 3) fasta_merge
echo "Running: fasta_merge -d assembly.maker.output/assembly_master_datastore_index.log -o assembly" | tee -a "${RUNLOG}"
/usr/bin/time -v "${MAKERBIN}/fasta_merge" -d assembly.maker.output/assembly_master_datastore_index.log -o assembly > fasta_merge.log 2>&1
echo "fasta_merge exit:$?" | tee -a "${RUNLOG}"

# List outputs
echo "Outputs (sizes):" | tee -a "${RUNLOG}"
ls -lh assembly.all.maker.gff assembly.all.maker.noseq.gff assembly.all.fasta 2>/dev/null || ls -lh assembly.*maker.* assembly.* 2>/dev/null || true | tee -a "${RUNLOG}"

echo "Merge run finished: $(date -u)" | tee -a "${RUNLOG}"

# Exit with success if the main files exist
if [ -s assembly.all.maker.gff ] && [ -s assembly.all.maker.noseq.gff ] && ( [ -s assembly.all.fasta ] || [ -s assembly.fasta ] ); then
  echo "MERGE SUCCESS" | tee -a "${RUNLOG}"
  exit 0
else
  echo "MERGE FAILED: check logs in ${WORKDIR}" >&2
  exit 3
fi
