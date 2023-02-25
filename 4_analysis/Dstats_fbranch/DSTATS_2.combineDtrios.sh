#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=convert_and_Dstat
#SBATCH --output=slurm-Dsuite_combineDtrios_%j.out
#SBATCH --error=slurm-Dsuite_combineDtrios_%j.error


###NOTE: this script can only be run once the first step (DSTATS_1) has been completed, as it relies on output from those analyses


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
dtools="/work/bs66/software/Dsuite/utils"
export PATH=$PATH:$dsuite


#specify invcfdir, invcf, outdir, faidx file, tree file, popset, test_trios
dstat_outdir="/work/bs66/dasanthera_novaseq/analysis/Dstats"
faidxfile="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta.fai"
treefile="/work/bs66/dasanthera_novaseq/analysis/Dstats/intree_dtrios.tre"
popset="/work/bs66/dasanthera_novaseq/analysis/Dstats/popset_dtrios.txt"
test_trios="/work/bs66/dasanthera_novaseq/analysis/Dstats/test_trios.txt"



##########
#ANALYSIS#
##########



#make array for specifying scaffolds -- used to specify input files, assuming they were generated on a per-scaffold basis
cd $dstat_outdir
chromlist=()
chromlist=$(awk '{print $1}' $faidxfile)
chromlist=(${chromlist[@]})
acount=${#chromlist[@]}


#for loop to append appropriate directories to the list
for((i=0;i<acount;i++));
do
    chromlist[i]="${chromlist[i]}".outfiles/"${chromlist[i]}"
done


#Combine results from Dtrios runs across scaffolds
Dsuite DtriosCombine -t $treefile -o WHOLEGENOME "${chromlist[@]}"


#do fbranch test in Dsuite for wholegenome
#note that this outputs only significant results at default p=0.01
#plots from the original paper put asterisks in significant blocks, so could
#change p-thresh to 1, and then put stars in blocks that are significant.
Dsuite Fbranch $treefile $dstat_outdir/WHOLEGENOME_combined_tree.txt > $dstat_outdir/WHOLEGENOME_Fbranch.txt

Dsuite Fbranch -Z $treefile WHOLEGENOME_combined_tree.txt > WHOLEGENOME_Fbranch.txt
python $dtools/dtools.py WHOLEGENOME_Fbranch.txt $treefile

#move wholegenome files to new directories
mkdir WHOLEGENOME.outfiles
mv WHOLEGENOME_*.txt WHOLEGENOME.outfiles


#plot fbranch test output
cd WHOLEGENOME.outfiles
python $dtools/dtools.py WHOLEGENOME_Fbranch.txt $treefile
mv fbranch.png WHOLEGENOME.fbranch.png
mv fbranch.svg WHOLEGENOME.fbranch.svg


