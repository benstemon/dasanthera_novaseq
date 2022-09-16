#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=quibl_evo2022


cd $SLURM_SUBMIT_DIR

#load conda env. gffread is in my env mapping_etc
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc



annotation="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.gff"

gffread $annotation --bed --keep-genes --sort-alpha -o annot_Pdavidsonii_gffread.genes.bed
