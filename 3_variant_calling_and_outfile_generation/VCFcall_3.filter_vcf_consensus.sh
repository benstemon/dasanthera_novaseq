#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p wessinger-48core
#SBATCH --job-name=create_consensus_vcf_pipeline


##NOTE: You can either run this script, or run 3a, b, and c. 3a, b, and c are redundant with this script, but run more quickly.

cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#activate conda environment with packages installed
#need bcftools (v1.15.1 works)
conda activate mapping_etc

#load vcftools (v0.1.17 on cluster works)
module load vcftools


#Location of unfiltered vcf, from script #2 in this pipeline.
invcf="/work/bs66/dasanthera_novaseq/VCFs/unfiltered_vcf.gz"


#path to the output directory for vcfs and gvcfs. Should be made already.
outdir="/work/bs66/dasanthera_novaseq/VCFs"


#number of cores used converting from gvcf to vcf
numthreads=16


########################

#####
#VCF FILTERING#
#####

#We must filter variant and invariant sites separately -- applying the same filters to all sites would remove invariant sites because of the way they are coded in the file.


#only invariant sites.
#max-maf 0 means no sites with minor allele frequency > 0.
#we will not implement any filters here (but you could feasibly do min depth filters)
#pipe into view because the compression is much faster than bgzip
 vcftools --gzvcf $invcf \
 --max-non-ref-ac 0 \
 --recode --recode-INFO-all --stdout | 
 bcftools view - \
 --threads $numthreads \
 -Oz -o $outdir/tmp-invariants.filtered.vcf.gz


#only variant sites. filters on QUALITY, DEPTH, and MISSINGNESS of reads.
#note that exact parameter values should be thought out on analysis-by-analysis basis. This script is for filtering with the eventual goal of generating multiple sequence alignments (phylogenomics!). So, some of the population-genomic filters common at the vcf filtering stage are not helpful here (e.g., filtering on minor allele frequency, HWE proportions, etc.).

#VCFTOOLS
#--min-alleles 2 means there must be at least two alleles (variant)
#--minGQ recode genotypes below GQ [int] as missing. GQ 20 = 1% chance the call is wrong
# meanDP filters remove sites with mean allelic depths outside specified limits [int]
#--minDP recodes any genotypes with fewer than [int] reads as missing
#--minQ filters sites below [int] quality threshold. Q20 = 1% chance there is no variant
#--max-missing removes sites on basis of missing data (1 = no missing data allowed)

#BCFTOOLS
# these commands soft filter variants with QUAL <20 and change GTs to reference sequence
vcftools --gzvcf $invcf \
 --non-ref-ac 1 \
 --min-alleles 2 \
 --minGQ 20 \
 --min-meanDP 3 \
 --max-meanDP 60 \
 --minDP 2 \
 --max-missing 0.50 \
 --recode --recode-INFO-all --stdout | 
 bcftools filter - \
 --soft-filter LowQual \
 --exclude '%QUAL<20' \
 --set-GTs 0 \
 --threads $numthreads \
 -Oz -o $outdir/tmp-variants.filtered.vcf.gz


#index both vcfs
tabix $outdir/tmp-invariants.filtered.vcf.gz
tabix $outdir/tmp-variants.filtered.vcf.gz


#concatenate the two vcfs and remove tmp files
bcftools concat $outdir/tmp-invariants.filtered.vcf.gz \
 $outdir/tmp-variants.filtered.vcf.gz \
 --threads $numthreads --allow-overlaps -Oz \
 -o $outdir/filtered_consensus_ready.vcf.gz

#you could uncomment this to remove the tmp files, but I prefer to remove by hand, in case something goes wrong with the script.
#rm $outdir/tmp-invariants.filtered.bcf.gz*
#rm $outdir/tmp-variants.filtered.bcf.gz*


#index final merged file
tabix $outdir/filtered_consensus_ready.vcf.gz
