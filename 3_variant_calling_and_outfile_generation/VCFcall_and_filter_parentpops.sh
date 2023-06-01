#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p wessinger-48core
#SBATCH --job-name=variants_bcftools


cd $SLURM_SUBMIT_DIR


##################
##### SETUP #####
##################

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with packages installed
#needs bcftools (v1.15.1 works)
conda activate mapping_etc

#load other modules
module load vcftools



#hard path to the base directory, reference genome, and mapped, filtered reads
basedir="/work/bs66/dasanthera_novaseq/"
refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
mapped_filtered_reads="/work/bs66/dasanthera_novaseq/mapped_filtered_bams_parentpops"


#change this to match number of cores allocated
numthreads=16

########################


##################
# call genotypes #
##################

#make list of mapped, filtered bams and out directory
bamlist=($mapped_filtered_reads/*.bam)
mkdir ${basedir}VCFs_parentpops
cd ${basedir}VCFs_parentpops


#Produce GT likelihoods, call variants, and normalize indels -> unfiltered .vcf
bcftools mpileup --threads $numthreads -Ou \
 -a FORMAT/AD,FORMAT/DP \
 --skip-indels \
 --min-BQ 20 \
 --min-MQ 20 \
 -f $refgenome ${bamlist[@]} | 
 bcftools call --threads $numthreads -m --gvcf 2 \
 -Oz -o unfiltered_vcf.gz


#index the vcf
tabix unfiltered_vcf.gz



##################
### filter VCF ###
##################

#We must filter variant and invariant sites separately -- applying the same filters to all sites would remove invariant sites because of the way they are coded in the file.


#only invariant sites.
 vcftools --gzvcf unfiltered_vcf.gz \
 --max-maf 0 \
 --recode --recode-INFO-all --stdout | 
 bcftools view - \
 --threads $numthreads \
 -Oz -o tmp-invariants.filtered.vcf.gz


#only variant sites. filters on QUALITY, DEPTH, and MISSINGNESS of reads.
#In addition, filter SNPs within 10 bp of indels
vcftools --gzvcf unfiltered_vcf.gz \
 --min-alleles 2 \
 --minGQ 20 \
 --min-meanDP 3 \
 --max-meanDP 60 \
 --minDP 2 \
 --max-missing 0.50 \
 --recode --recode-INFO-all --stdout | 
 bcftools filter - \
 --SnpGap 10 \
 --soft-filter LowQual \
 --exclude '%QUAL<20' \
 --set-GTs 0 \
 --threads $numthreads \
 -Oz -o tmp-variants.filtered.vcf.gz


#index both vcfs
tabix tmp-invariants.filtered.vcf.gz
tabix tmp-variants.filtered.vcf.gz


#concatenate the two vcfs and remove tmp files
bcftools concat tmp-invariants.filtered.vcf.gz tmp-variants.filtered.vcf.gz \
 --threads $numthreads --allow-overlaps -Oz -o filtered_consensus_ready.vcf.gz

#you could uncomment this to remove the tmp files, but I prefer to remove by hand, in case something goes wrong with the script.
#rm $outdir/tmp-invariants.filtered.bcf.gz*
#rm $outdir/tmp-variants.filtered.bcf.gz*


#index final merged file
tabix filtered_consensus_ready.vcf.gz




