#!/bin/bash
#SBATCH --time=7-0
#SBATCH --mem=120G
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --job-name=Maker
#SBATCH --output=/data/users/afairman/euk_org_ann/logs/Maker_gene_annotation_%j.out
#SBATCH --error=/data/users/afairman/euk_org_ann/logs/Maker_gene_annotation_%j.err

# Define directories
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/afairman/euk_org_ann/results/Maker"
REPEATMASKER_DIR="/data/courses/assembly-annotation-course/CDS_annotation/softwares/RepeatMasker"

# Add RepeatMasker to PATH
export PATH=$PATH:"${REPEATMASKER_DIR}"

# Load required modules
module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a

# Change to working directory
cd ${WORKDIR}

# Run MAKER with MPI
mpiexec --oversubscribe -n 50 apptainer exec \
--bind $SCRATCH:/TMP \
--bind $COURSEDIR \
--bind $AUGUSTUS_CONFIG_PATH \
--bind $REPEATMASKER_DIR \
--bind /data/users/afairman \
${COURSEDIR}/containers/MAKER_3.01.03.sif \
maker -mpi --ignore_nfs_tmp -TMP /TMP maker_opts.ctl maker_bopts.ctl maker_evm.ctl maker_exe.ctl
