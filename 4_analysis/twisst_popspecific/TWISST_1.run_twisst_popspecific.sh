#!/bin/sh

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p wessinger-48core
#SBATCH --job-name=TWISST
#SBATCH --output=slurm-TWISST_popspecific_%j.out


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
outdir="/work/bs66/dasanthera_novaseq/analysis/twisst_popspecific"


##########
#ANALYSIS#
##########

#run twisst for rup86 tests
python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_rup86_t1.weights.txt.gz \
 -g new -g car -g rup86 -g mon --outgroup mon \
 --method complete --groupsFile $outdir/groupsfile_test-rup86.txt

python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_rup86_t2.weights.txt.gz \
 -g new -g car -g rup -g mon --outgroup mon \
 --method complete --groupsFile $outdir/groupsfile_test-rup86.txt



#run twisst for dav116 tests
python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_dav116_t1.weights.txt.gz \
 -g ncr -g fru -g dav116 -g mon --outgroup mon \
 --method complete --groupsFile $outdir/groupsfile_test-dav116.txt

python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_dav116_t2.weights.txt.gz \
 -g ncr -g fru -g dav -g mon --outgroup mon \
 --method complete --groupsFile $outdir/groupsfile_test-dav116.txt


#run twisst for rup101 tests
python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_rup101_t1.weights.txt.gz \
 -g dav -g fru -g rup101 -g mon --outgroup mon \
 --method complete --groupsFile $outdir/groupsfile_test-rup101.txt

python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_rup101_t2.weights.txt.gz \
 -g dav -g fru -g rup -g mon --outgroup mon \
 --method complete --groupsFile $outdir/groupsfile_test-rup101.txt



#run twisst for new-dav tests
python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_new-dav_t1.weights.txt.gz \
 -g dav -g fru -g new -g mon --outgroup mon \
 --method complete --groupsFile $outdir/groupsfile_test-new-dav.txt

python $twisst/twisst.py -t $intrees \
 -w $outdir/twisst_new-dav_t2.weights.txt.gz \
 -g dav -g fru -g car -g mon --outgroup mon \
 --method complete --groupsFile $outdir/groupsfile_test-new-dav.txt


