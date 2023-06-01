#!/bin/sh

#SBATCH -N 1
#SBATCH -n 6
#SBATCH -p wessinger-48core
#SBATCH --job-name=introtests
#SBATCH --output=slurm-Dsuite_introtests_1kb_%a.out
#SBATCH --error=slurm-Dsuite_introtests_1kb_%a.error


cd $SLURM_SUBMIT_DIR



#######
#SETUP#
#######



#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#activate conda environment with packages installed
#need python3 with numpy and other Dsuite dependencies installed
conda activate mapping_etc


#set Dsuite variable to path, and for f-branch test, the path to directory with dtools.py
dsuite="/work/bs66/software/Dsuite/Build/"
export PATH=$PATH:$dsuite


#specify invcf, outdir. Also specify the popset and triosets (same order)
invcf="/work/bs66/dasanthera_novaseq/VCFs/filtered_consensus_ready_no-indels.vcf.gz"
dstat_outdir="/work/bs66/dasanthera_novaseq/analysis/Dstats_CNRD"
popset="popset_CNRD.txt"
trioset="trioset_CNRD.txt"



#run Dsuite Dinvestigate
cd $dstat_outdir
Dsuite Dinvestigate -n 1kb -w 1000,250 $invcf $popset $trioset


