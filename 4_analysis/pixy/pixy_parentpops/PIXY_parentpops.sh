#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=pixy_fullgenome_fullspecies
#SBATCH --output=slurm-pixy_fullgenome_fullspecies_%j.out


cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with pixy installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc


#paths to (1) the input vcf, (2) the desired out-directory, and (3) the populations file
invcf="/work/bs66/dasanthera_novaseq/VCFs_parentpops/filtered_consensus_ready.vcf.gz"
outdir="/work/bs66/dasanthera_novaseq/analysis/pixyout_parentpops"
popfile="/work/bs66/dasanthera_novaseq/analysis/pixyout_parentpops/popfile_parentpops.txt"


#Run pixy 10kb
pixy --stats pi fst dxy \
 --vcf $invcf \
 --populations $popfile \
 --window_size 10000 \
 --n_cores 8 \
 --output_folder $outdir \
 --output_prefix pixy_parentpops_10kb


#Run pixy 50kb
pixy --stats pi fst dxy \
 --vcf $invcf \
 --populations $popfile \
 --window_size 50000 \
 --n_cores 8 \
 --output_folder $outdir \
 --output_prefix pixy_parentpops_50kb

