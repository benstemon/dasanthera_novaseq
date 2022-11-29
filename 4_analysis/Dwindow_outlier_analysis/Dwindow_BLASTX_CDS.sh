#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=blastx_Dwindow_CDS
#SBATCH --output=slurm-blastx_Dwindow_CDS.out


cd $SLURM_SUBMIT_DIR
module load blast


#specify the directory with individual CDS fastas from the reference genome:
CDSdir="/work/bs66/project_compare_genomes/davidsonii_CDS"

#specify the local blast database
blastdb="/work/bs66/project_compare_genomes/blast_database/db_Swissprot/swissprot"

#specify output directory and desired unique CDS output file
outdir="/work/bs66/dasanthera_novaseq/analysis/Dwindow_outlier_analysis"
unique_CDS="unique_CDS_1000_500_z4.txt"




#write all unique CDS and the number of times they occur to a new file:
cd $outdir
awk -v OFS='\t' '{ count[$8]++ } END { for (word in count) print word, count[word]}' CDS-hits_only_1000_500_z4.bed > $unique_CDS


#re-read this list in, generating an array to blast
blastlist=()
while read name num;
do
    blastlist+=$name" "
done < $unique_CDS


#iteratively blast each CDS in the array
for i in $blastlist;
do
    blastx -query $CDSdir/$i.fa -db $blastdb -outfmt 6 -max_target_seqs 1 -num_threads 4 -evalue 1e-50 | awk -v OFS='\t' -v CDS="$i" '{$1=CDS; print}' >> $outdir/results_blastx_1000_500_z4.txt
done


