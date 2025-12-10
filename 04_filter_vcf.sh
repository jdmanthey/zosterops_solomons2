#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-41

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/02_zosterops

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# filter based on missing data for each of the subsets of data
# allow up to 5 individuals missing per site (89%+ completeness per site)

# for stats and phylogenies
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --max-missing 0.89 \
--max-alleles 2 --max-maf 0.49 --recode --recode-INFO-all \
--out ${workdir}/06_phylogenies/${region_array}

# for gene flow
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --max-missing 0.89 \
--min-alleles 2 --max-alleles 2 --mac 3 --max-maf 0.49 --recode --recode-INFO-all \
--out ${workdir}/07_dstats/${region_array}


# zip and index the stats files
bgzip ${workdir}/06_phylogenies/${region_array}.recode.vcf

tabix ${workdir}/06_phylogenies/${region_array}.recode.vcf.gz

