#!/bin/sh

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p wessinger-48core
#SBATCH --job-name=TWISST
#SBATCH --output=slurm-TWISST_%j.out


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
outdir="/work/bs66/dasanthera_novaseq/analysis/twisst_analyses"
groupsfile="/work/bs66/dasanthera_novaseq/analysis/twisst_analyses/twisst_groupsfile.txt"


##########
#ANALYSIS#
##########

#run twisst on the NCR clade only
python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_NCR_10kbtrees.weights.txt.gz \
 -g new -g car -g rup -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile


#run twisst and include P. davidsonii (NCRD)
python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_NCRD_10kbtrees.weights.txt.gz \
 -g new -g car -g rup -g dav -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

