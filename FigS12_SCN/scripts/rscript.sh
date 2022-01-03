#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-06:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=100G                         # Memory total in MB (for all cores)
#SBATCH -o slurm_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e slurm_%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like

module load gcc/6.2.0 
module load R/4.0.1

Rscript scnet_organoid.r


