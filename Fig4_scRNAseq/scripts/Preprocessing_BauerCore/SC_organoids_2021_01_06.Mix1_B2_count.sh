#!/bin/bash
#SBATCH --time=4-18:00:00
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --partition=bos-info_priority
#SBATCH --job-name=count_10x.210120_A00794_0355_BHM3L3DSXY_Mix1_B2
#SBATCH --output=log/count_10x.Mix1_B2-%j.out
#SBATCH --error=log/count_10x.Mix1_B2-%j.err
#SBATCH --exclusive
#SBATCH --chdir=/n/boslfs02/LABS/informatics/Lab/sequencing/analysis/210120_A00794_0355_BHM3L3DSXY
#SBATCH --nice=2147483645
ulimit -u $(ulimit -Hu)
exit_code=0
mkdir -p /sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/count/Mix1_B2
find /scratch -maxdepth 1 -name 'count_*' -user $(id -un) -exec rm -rf {} +
mkdir -p /scratch/count_210120_A00794_0355_BHM3L3DSXY_Mix1_B2_$SLURM_JOB_ID
cd /scratch/count_210120_A00794_0355_BHM3L3DSXY_Mix1_B2_$SLURM_JOB_ID
echo '{ "SC_RNA_COUNTER_CS.SC_MULTI_CS.SC_MULTI_CORE.MULTI_GEM_WELL_PROCESSOR.COUNT_GEM_WELL_PROCESSOR._BASIC_SC_RNA_COUNTER._MATRIX_COMPUTER.ALIGN_AND_COUNT": { "chunk.mem_gb": 64 } }' > count-overrides.json
env time -v cellranger count --project=SC_organoids_2021_01_06 --id=Mix1_B2 --transcriptome=/ref/refdata-gex-GRCh38-2020-A --sample=Mix1_B2 --overrides=count-overrides.json --fastqs=/sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/fastq --localmem=$((9*$(ulimit -m)/10000000)) --localcores=$SLURM_JOB_CPUS_PER_NODE || exit_code=$?

env time -v cp -Rp Mix1_B2/*.mri.tgz Mix1_B2/outs /sequencing/analysis/210120_A00794_0355_BHM3L3DSXY/count/Mix1_B2/ || exit_code=$((exit_code | $?))
rm -rf /scratch/count_210120_A00794_0355_BHM3L3DSXY_Mix1_B2_$SLURM_JOB_ID
exit $exit_code
