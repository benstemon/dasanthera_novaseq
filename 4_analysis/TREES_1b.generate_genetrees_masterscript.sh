#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=generate_genetrees

#this is a masterscript to generate sliding window trees 


cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with pixy installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#paths to important scripts and directories
pythonscript="/work/bs66/dasanthera_novaseq/analysis/TREES_1a.create_fasta_window_alignments.py"
scaffold_dir="/work/bs66/dasanthera_novaseq/consensus_alignments/scaffold_fullgenome_fastas"
outdir="/work/bs66/dasanthera_novaseq/analysis/genetree_infiles"


#list of scaffold names
scaflist=("1085" "1086" "1087" "2151" "2446" "2531" "2532" "2533" "2684" "2685" "2686" "2687")


#run for loop for python script to generate gene tree files for each scaffold
#syntax: python3 scriptname input_fasta windowsize prefix
for i in "${scaflist[@]}"
do
	python3 $pythonscript $scaffold_dir/allseqs_consensus_scaffold_$i.fa 10000 $outdir/scaffold_$i
done

