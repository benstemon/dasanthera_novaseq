#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --out=slurm-saguaro_%j_%a.out
#SBATCH -p wessinger-48core
#SBATCH --job-name=saguaro_convert

# This script converts alignment data for the Dasanthera species into Saguaruo binary format
# This is an array script that must be submitted using sbatch --array [0-n], where n = # input files - 1
# The array setup allows each chromosome alignment to be run independently in parallel

# INPUTS: input to the script are the rearranged full-genome fasta files extracted from the full filtered VCF. These should have one entry per sample and be separated by chromsome (e.g. LG01.rearranged.fa)

cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#path to saguaro executables
saguaro="/work/bs66/software/saguarogw-code-r44"

#use the fasta alignments for each chromosome as inputs after using the Fasta2HMMFeature to convert them to Saguaro binary features files. Then run the program separately on each chromosome in parallel as an array job.


# Make array for specifying the linkage groups 
# The desired input is a multi-fasta alignment for each scaffold of interest.
mafdir="/work/bs66/dasanthera_novaseq/consensus_alignments/scaffold_fullgenome_fastas"
chromlist=$(ls $mafdir/*.fa)
chromlist=(${chromlist[@]})
chromfile="${chromlist[$SLURM_ARRAY_TASK_ID]}"

chromname="$(basename $chromfile | sed 's/.fa//g' | sed 's/allseqs_consensus_//g')"



# Set the output files
outdir="/work/bs66/dasanthera_novaseq/analysis/saguaro_analyses"
outbin="$outdir/$chromname.saguaro.features"


##########
#ANALYSIS#
##########

# Run Saguaro to convert to features
$saguaro/./Fasta2HMMFeature -i $chromfile -o $outbin -nosame


#set up output directories for each scaffold
mkdir -p "$outdir/$chromname"
newoutdir="$outdir/$chromname"
cd $newoutdir


#now that features have been created, run Saguaro
$saguaro/./Saguaro -f $outbin -o $newoutdir -iter 40

