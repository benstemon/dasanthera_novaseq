#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=blastx_fdmoutlier
#SBATCH --output=slurm-fdm_outlier_CDS_%j_%a.out


cd $SLURM_SUBMIT_DIR
module load blast


#specify the directory with individual CDS fastas from the reference genome:
CDSdir="/work/bs66/project_compare_genomes/davidsonii_CDS"

#specify the local blast database
blastdb="/work/bs66/project_compare_genomes/blast_database/db_Swissprot/swissprot"

#specify directory containing outlier .bed files
outlierdir="/work/bs66/dasanthera_novaseq/analysis/fdm_outlier_analyses/fdm_outlier_CDS_hits"

#specify output directory and desired unique CDS output file
outdir="/work/bs66/dasanthera_novaseq/analysis/fdm_outlier_analyses"


#set up the array -- list the different CDS hits here:
CDSlist=('CDS_hits_10kb.bed' 'CDS_hits_1kb.bed' 'CDS_hits_500bp.bed' 'CDS_hits_5kb.bed')
inCDS="${CDSlist[$SLURM_ARRAY_TASK_ID]}"
unique_CDS="${inCDS/.bed/_UNIQUECDS.txt}"


#write all unique CDS and the number of times they occur to a new file:
#here, gene model name occurs on the 13th column
cd $outdir
awk -v OFS='\t' '{ count[$13]++ } END { for (word in count) print word, count[word]}' $outlierdir/$inCDS > $unique_CDS


#re-read this list in, generating an array to blast
blastlist=()
while read name num;
do
    blastlist+=$name" "
done < $unique_CDS


#iteratively blast each CDS in the array
for i in $blastlist;
do
    blastx -query $CDSdir/$i.fa -db $blastdb -outfmt 6 -max_target_seqs 1 -num_threads 4 -evalue 1e-50 | awk -v OFS='\t' -v CDS="$i" '{$1=CDS; print}' >> $outdir/results_$inCDS.txt
done


