#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=variants_bcftools


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with packages installed
#needs bcftools (v1.15.1 works)
conda activate mapping_etc


#hard path to the reference genome and mapped, filtered reads
#mapped filtered reads should be same directory as output from script: mapping_pipeline_array_s1a-f.sh
refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
mapped_filtered_reads="/work/bs66/dasanthera_novaseq/mapped_filtered_bams"


#path to the output directory for vcf files (called variants). Make prior to running
outdir="/work/bs66/dasanthera_novaseq/GVCF_VCF_unfiltered"


#number of cores used for mpileup and call
#change this to match number of cores allocated
numthreads=8
########################


#make list of mapped, filtered bams
bamlist=($mapped_filtered_reads/*.bam)


#Produce GT likelihoods, call variants, and normalize indels -> unfiltered .vcf
#mpileup produces genotype likelihoods from bam files. Filtering for BQ and MQ > 20
#call calls SNPs and indels from the genotype likelihoods
#norm normalizes and left-aligns indels
bcftools mpileup --threads $numthreads -Ou \
 -a FORMAT/AD,FORMAT/DP \
 --min-BQ 20 \
 --min-MQ 20 \
 --gvcf 0 -f $refgenome ${bamlist[@]} | 
 bcftools call --threads $numthreads -m --gvcf 0 -Ou | 
 bcftools norm --threads $numthreads -f $refgenome \
 -Oz -o $outdir/unfiltered_gvcf.gz


#index the gvcf
tabix $outdir/unfiltered_gvcf.gz


#convert from gvcf to a vcf (all sites, but invariants properly assessed)
bcftools convert --gvcf2vcf \
 --fasta-ref $refgenome \
 --threads $numthreads \
 $outdir/unfiltered_gvcf.gz | 
 bgzip -c --threads $numthreads > $outdir/filtered_vcf.gz

#index the vcf
tabix $outdir/unfiltered_vcf.gz

