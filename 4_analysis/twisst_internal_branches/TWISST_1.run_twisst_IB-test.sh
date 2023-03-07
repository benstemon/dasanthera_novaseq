#!/bin/sh

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p wessinger-48core
#SBATCH --job-name=TWISST
#SBATCH --output=slurm-TWISST_IBtest_%j.out


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate ete3_twisst


#specify invcfdir, invcf, outdir, faidx file, tree file
twisst="/work/bs66/software/twisst"
intrees="/work/bs66/dasanthera_novaseq/analysis/treemetrics/combined_10kbwindowtrees.tre"
outdir="/work/bs66/dasanthera_novaseq/analysis/twisst_IB-test"


##########
#ANALYSIS#
##########

#run twisst for rup86 tests
python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_IBL-test.weights.txt.gz \
 -g new -g car -g rup -g dav -g fru -g mon --outgroup mon \
 --method complete --groupsFile $outdir/groupsfile_IB-test.txt

