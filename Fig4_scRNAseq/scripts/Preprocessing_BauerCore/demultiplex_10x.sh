#!/bin/bash
#SBATCH --time=4-20:00:00
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --partition=bos-info_priority
#SBATCH --job-name=demultiplex_10x.210120_A00794_0355_BHM3L3DSXY
#SBATCH --output=log/demultiplex_10x-%j.out
#SBATCH --error=log/demultiplex_10x-%j.err
#SBATCH --exclusive
#SBATCH --chdir=/n/boslfs02/LABS/informatics/Lab/sequencing/analysis/210120_A00794_0355_BHM3L3DSXY
ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)
exit_code=0
mkdir -p /sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/fastq
mkdir -p /scratch/210120_A00794_0355_BHM3L3DSXY_fastq_$SLURM_JOB_ID
cd /scratch/210120_A00794_0355_BHM3L3DSXY_fastq_$SLURM_JOB_ID
env time -v cellranger mkfastq --run=/sequencing/source/210120_A00794_0355_BHM3L3DSXY --samplesheet=/sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/SampleSheet.csv --output-dir=/sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/fastq --interop-dir=/sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/InterOp --filter-single-index --localmem=$((9*$(ulimit -m)/10000000)) --loading-threads=$((SLURM_JOB_CPUS_PER_NODE/4)) --writing-threads=$((SLURM_JOB_CPUS_PER_NODE/4)) --processing-threads=$SLURM_JOB_CPUS_PER_NODE --localcores=$SLURM_JOB_CPUS_PER_NODE --barcode-mismatches=0 || exit_code=$?
cp -p */*.mri.tgz /sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/fastq || exit_code=$((exit_code | $?))
rm -rf /scratch/210120_A00794_0355_BHM3L3DSXY_fastq_$SLURM_JOB_ID
exit $exit_code
