#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=introtests
#SBATCH --output=slurm-Dsuite_introtests_5kb_%a.out
#SBATCH --error=slurm-Dsuite_introtests_5kb_%a.error


###NOTE: this is an array batch script. Set 0-n where n = number of popset files -1.


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
dstat_outdir="/work/bs66/dasanthera_novaseq/analysis/Dstats_introtests"

poplist=("popset_dav-new.txt" "popset_dav116-ncr.txt" "popset_rup86-new.txt" "popset_rup101-fru.txt")
triolist=("trioset_dav-new.txt" "trioset_dav116-ncr.txt" "trioset_rup86-new.txt" "trioset_rup101-fru.txt")



#specify pop and trioset from array, and make variable for run-name
poplist=(${poplist[@]});popset="${poplist[$SLURM_ARRAY_TASK_ID]}"

triolist=(${triolist[@]});trioset="${triolist[$SLURM_ARRAY_TASK_ID]}"

run_name="${popset//popset_/}";run_name="${run_name//.txt/}"


#run Dsuite Dinvestigate
cd $dstat_outdir
Dsuite Dinvestigate -n 5kb_$run_name -w 5000,1250 $invcf $popset $trioset


