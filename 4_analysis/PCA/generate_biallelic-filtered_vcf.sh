#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=biallelic_LD-filter
#SBATCH --output=slurm-biallelic_LD_VCFfilter_%j.out


cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with pixy installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#set up infiles
invcf="/work/bs66/dasanthera_novaseq/VCFs/filtered_consensus_ready_no-indels.vcf.gz"
outvcf="/work/bs66/dasanthera_novaseq/VCFs/PCA-ready_biallelic_ld-filtered.vcf"

#bcftools commands to generate biallelic vcf and filter for LD
bcftools view --threads 4 --exclude 'FILTER="LowQual"' \
 --min-alleles 2 --max-alleles 2 --types snps $invcf -Ou |
 bcftools +prune --threads 4 -m 0.1 -w 50kb -Ov -o $outvcf

