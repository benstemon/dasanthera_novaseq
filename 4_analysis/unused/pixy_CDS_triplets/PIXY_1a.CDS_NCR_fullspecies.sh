#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=pixy_CDS_NCR_fullspecies
#SBATCH --output=slurm-pixy_CDS_NCR_fullspecies_%j.out


cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with pixy installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc


#paths to (1) the input vcf, (2) the desired out-directory, and (3) the populations file
invcf="/work/bs66/dasanthera_novaseq/VCFs/filtered_consensus_ready_no-indels.vcf.gz"
outdir="/work/bs66/dasanthera_novaseq/analysis/pixyout_CDS_triplets"
popfile="/work/bs66/dasanthera_novaseq/analysis/pixyout_CDS_triplets/popfile_1a.NCR_fullspecies.txt"
bedfile="/work/bs66/project_compare_genomes/annot_Pdavidsonii_1mb.gffread.genes.bed"


#Run pixy
pixy --stats dxy \
 --vcf $invcf \
 --populations $popfile \
 --bed_file $bedfile \
 --n_cores 4 \
 --output_folder $outdir \
 --output_prefix NCR_fullspecies

