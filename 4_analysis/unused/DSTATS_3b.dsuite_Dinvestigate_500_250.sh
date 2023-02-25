#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=convert_and_Dstat
#SBATCH --output=slurm-Dsuite_Dinvestigate_%j.out
#SBATCH --error=slurm-Dsuite_Dinvestigate_%j.error


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######



#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with packages installed
#need bcftools (v1.15.1 works)
conda activate mapping_etc


#set Dsuite variable to path, and for f-branch test, the path to directory with dtools.py
dsuite="/work/bs66/software/Dsuite/Build/"
dtools="/work/bs66/software/Dsuite/utils"
export PATH=$PATH:$dsuite


#specify invcfdir, invcf, outdir, faidx file, tree file, popset, test_trios
faidxfile="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta.fai"
invcf="/work/bs66/dasanthera_novaseq/VCFs/filtered_consensus_ready_no-indels.vcf.gz"
popset="/work/bs66/dasanthera_novaseq/analysis/Dstats/popset_dtrios.txt"
test_trios="/work/bs66/dasanthera_novaseq/analysis/Dstats/test_trios_windows.txt"



##########
#ANALYSIS#
##########



#make array for specifying scaffolds -- to pull out appropriate scaffold, given array
chromlist=$(awk '{print $1}' $faidxfile)
chromlist=(${chromlist[@]})
chromname="${chromlist[$SLURM_ARRAY_TASK_ID]}"


#perform investigate analyses -- try large windows first, and on whole genome 
Dsuite Dinvestigate -w 500,250 $invcf $popset $test_trios


