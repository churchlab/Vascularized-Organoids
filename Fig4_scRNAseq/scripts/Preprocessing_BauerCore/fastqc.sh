#!/bin/bash
#SBATCH --time=4-12:00:00
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --partition=bos-info_priority
#SBATCH --job-name=fastqc.210120_A00794_0355_BHM3L3DSXY
#SBATCH --output=log/fastqc-%j.out
#SBATCH --error=log/fastqc-%j.err
#SBATCH --exclusive
#SBATCH --chdir=/n/boslfs02/LABS/informatics/Lab/sequencing/analysis/210120_A00794_0355_BHM3L3DSXY
ulimit -u $(ulimit -Hu)
mkdir -p /sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/QC
cd /sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/fastq/
find . -name '*.fastq.gz' ! -name 'Undetermined*' -exec env time -v fastqc -o /sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/QC --threads $SLURM_JOB_CPUS_PER_NODE {} +
