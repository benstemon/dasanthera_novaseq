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
outdir="/work/bs66/dasanthera_novaseq/analysis/twisst_fullspecies"
groupsfile="/work/bs66/dasanthera_novaseq/analysis/twisst_fullspecies/twisst_fullspecies_groupsfile.txt"


##########
#ANALYSIS#
##########

#iterate over each of the ten triplets
#NCR
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_new_car_rup.txt.gz \
 -g new -g car -g rup -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

#NCF
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_new_car_fru.txt.gz \
 -g new -g car -g fru -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

#NCD
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_new_car_dav.txt.gz \
 -g new -g car -g dav -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

#NRF
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_new_rup_fru.txt.gz \
 -g new -g rup -g fru -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

#NRD
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_new_rup_dav.txt.gz \
 -g new -g rup -g dav -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

#FDN
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_fru_dav_new.txt.gz \
 -g fru -g dav -g new -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

#CRF
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_car_rup_fru.txt.gz \
 -g car -g rup -g fru -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

#CRD
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_car_rup_dav.txt.gz \
 -g car -g rup -g dav -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

#FDC
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_fru_dav_car.txt.gz \
 -g fru -g dav -g car -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile

#FDR
python $twisst/twisst.py -t $intrees \
 -w $outdir/fullspecies_fru_dav_rup.txt.gz \
 -g fru -g dav -g rup -g mon --outgroup mon \
 --method complete --groupsFile $groupsfile
