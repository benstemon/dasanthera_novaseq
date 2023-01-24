#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH --out=slurm-CFs_IQtree.%j.out
#SBATCH -p wessinger-48core
#SBATCH --job-name=CFs_IQtree


cd $SLURM_SUBMIT_DIR

#load conda env with iqtree installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#paths to important files and directories
reftree_10kb="/work/bs66/dasanthera_novaseq/analysis/astral_trees/astral_10kb_noannotations.tre"
reftree_CDS="/work/bs66/dasanthera_novaseq/analysis/astral_trees/astral_CDS_noannotations.tre"

sourcetrees_10kb="/work/bs66/dasanthera_novaseq/analysis/treemetrics/combined_10kbwindowtrees.tre"
sourcetrees_CDS="/work/bs66/dasanthera_novaseq/analysis/treemetrics/combined_CDStrees.tre"

sourcealignments_10kb="/work/bs66/dasanthera_novaseq/analysis/genetree_infiles/allscaf_alignments"
sourcealignments_CDS="/work/bs66/dasanthera_novaseq/analysis/CDS_genetree_infiles_min3"

outdir="/work/bs66/dasanthera_novaseq/analysis/treemetrics/concordance_factors"


#run IQtree gene concordance factors
cd $outdir
iqtree -t $reftree_10kb --gcf $sourcetrees_10kb --cf-verbose --scf 100 -p $sourcealignments_10kb --df-tree --prefix 10kbsource_10kbref -T 4

iqtree -t $reftree_CDS --gcf $sourcetrees_CDS --cf-verbose --scf 100 -p $sourcealignments_CDS --df-tree --prefix CDSsource_CDSref -T 4

