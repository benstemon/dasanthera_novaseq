#!/bin/sh

#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p wessinger-48core
#SBATCH --job-name=QC_s1

cd $SLURM_SUBMIT_DIR


#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with QC packages installed
conda activate QC


################
#fastqc+multiqc#
################

rawreads="/work/bs66/dasanthera_novaseq/raw_reads"

fastqc_outdir_rawreads="/work/bs66/dasanthera_novaseq/fastqc_rawreads"
rm -r $fastqc_outdir_rawreads
mkdir $fastqc_outdir_rawreads

#move to directory with reads (raw reads directory)
cd $rawreads

#store fastq files as array
files=(*.fastq.gz)

#perform fastqc -- distinction from for loop is this can process -t files simultaneously
fastqc "${files[@]}" -t 20 -o $fastqc_outdir_rawreads

#summarize results with multiqc
multiqc $fastqc_outdir_rawreads -o $fastqc_outdir_rawreads


