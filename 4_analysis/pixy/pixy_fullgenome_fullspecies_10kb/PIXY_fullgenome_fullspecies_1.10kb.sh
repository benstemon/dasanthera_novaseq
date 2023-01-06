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
invcf="/work/bs66/dasanthera_novaseq/VCFs/filtered_consensus_ready_no-indels.vcf.gz"
outdir="/work/bs66/dasanthera_novaseq/analysis/pixyout_fullgenome_fullspecies_10kb"
popfile="/work/bs66/dasanthera_novaseq/analysis/pixyout_fullgenome_fullspecies_10kb/popfile_allpops_fullspecies.txt"


#Run pixy
pixy --stats pi fst dxy \
 --vcf $invcf \
 --populations $popfile \
 --window_size 10000 \
 --n_cores 8 \
 --output_folder $outdir \
 --output_prefix fullgenome_fullspecies_10kb

