#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --partition=bos-info_priority
#SBATCH --job-name=checksum.210120_A00794_0355_BHM3L3DSXY
#SBATCH --output=log/checksum-%j.out
#SBATCH --error=log/checksum-%j.err
#SBATCH --exclusive
#SBATCH --chdir=/n/boslfs02/LABS/informatics/Lab/sequencing/analysis/210120_A00794_0355_BHM3L3DSXY
set -o errexit -o pipefail
cd /sequencing/analysis/210120_A00794_0355_BHM3L3DSXY
find fastq/ -name '*.fastq.gz' -print0 | env time -v xargs -0 -n 1 -P $SLURM_JOB_CPUS_PER_NODE md5sum > md5sum.txt.unsorted
sort -k 2,2 -o md5sum.txt md5sum.txt.unsorted
rm md5sum.txt.unsorted

