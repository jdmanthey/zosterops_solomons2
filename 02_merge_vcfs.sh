#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-41

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/02_zosterops

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# run bcftools to merge the vcf files
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/*vcf.gz > ${workdir}/04_vcf/${region_array}.vcf


