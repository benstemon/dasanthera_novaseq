#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=convert_and_Dstat
#SBATCH --output=slurm-Dsuite_Dtrios_%j.out
#SBATCH --error=slurm-Dsuite_Dtrios_%j.error


###NOTE: this is set up as a batch array script. The array number is the number of scaffold.fasta files, -1
##0-11 for davidsonii 1mb


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
invcfdir="/work/bs66/dasanthera_novaseq/VCFs"
invcf="/work/bs66/dasanthera_novaseq/VCFs/filtered_consensus_ready_no-indels.vcf.gz"
dstat_outdir="/work/bs66/dasanthera_novaseq/analysis/Dstats"
faidxfile="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta.fai"
treefile="/work/bs66/dasanthera_novaseq/analysis/Dstats/intree_dtrios.tre"
popset="/work/bs66/dasanthera_novaseq/analysis/Dstats/popset_dtrios.txt"



##########
#ANALYSIS#
##########



#make array for specifying scaffolds -- to pull out appropriate scaffold, given array
chromlist=$(awk '{print $1}' $faidxfile)
chromlist=(${chromlist[@]})
chromname="${chromlist[$SLURM_ARRAY_TASK_ID]}"


#use bcftools to make scaffold-specific .vcf file
bcftools view $invcf --regions $chromname -o $invcfdir/$chromname.vcf.gz -Oz

#index the vcf
tabix $invcfdir/$chromname.vcf.gz


#Use Dsuite to estimate Dstatistics between all triplet pairs
#specifying the output prefix, the input vcf, and the popset
Dsuite Dtrios -o $dstat_outdir/$chromname -t $treefile $invcfdir/$chromname.vcf.gz $popset


#do fbranch test in Dsuite
#note that this outputs only significant results at default p=0.01
#plots from the original paper put asterisks in significant blocks, so could
#change p-thresh to 1, and then put stars in blocks that are significant.
#Dsuite Fbranch $treefile $dstat_outdir/"${chromname}_tree.txt" > $dstat_outdir/"${chromname}_Fbranch.txt"

#move files to new directories
cd $dstat_outdir
mkdir $chromname.outfiles
mv "${chromname}"_* $chromname.outfiles

#plot fbranch test output
#cd $chromname.outfiles
#python $dtools/dtools.py "${chromname}_Fbranch.txt" $treefile
#mv fbranch.png $chromname.fbranch.png
#mv fbranch.svg $chromname.fbranch.svg


