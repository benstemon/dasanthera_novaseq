#!/bin/sh

#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p wessinger-48core
#SBATCH --job-name=QC_s3-5

cd $SLURM_SUBMIT_DIR


#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with QC packages installed
conda activate QC


#######
#fastp#
#######

#path to illumina mergedreads
mergedreads="/work/bs66/dasanthera_novaseq/merged_reads"

#path to fastp outfiles (creates directory, and replaces is already present)
out_fastp="/work/bs66/dasanthera_novaseq/filtered_reads"
rm -r $out_fastp
mkdir $out_fastp


#fastp for loop
for r1in in $mergedreads/*_R1_001.fastq.gz; 
do
    r2in="${r1in/R1_001.fastq.gz/R2_001.fastq.gz}"
    r1out="${r1in##*/}"
    r2out="${r1out/R1_001.fastq.gz/R2_001.fastq.gz}"
    fastp -i "$r1in" -I "$r2in" --out1 "${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R1_001.fastq.gz}" --out2 "${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R2_001.fastq.gz}" --unpaired1 "${r1out/merged_L001_R1_001.fastq.gz/unpaired_L001_R1_001.fastq.gz}" --unpaired2 "${r1out/merged_L001_R1_001.fastq.gz/unpaired_L001_R2_001.fastq.gz}" -x -c -w 16 -h "${r1out/merged_L001_R1_001.fastq.gz/html}" -j "${r1out/merged_L001_R1_001.fastq.gz/json}"
done

#move outfiles to correct directory (probably better way to do this)
mv $mergedreads/*trimmed* $out_fastp
mv $mergedreads/*unpaired* $out_fastp
mv $mergedreads/*html $out_fastp
mv $mergedreads/*json $out_fastp



################
#fastqc+multiqc#
################

out_fastqc="/work/bs66/dasanthera_novaseq/fastqc_filtered"
rm -r $out_fastqc
mkdir $out_fastqc

#move to directory with reads (fastp outfile directory)
cd $out_fastp

#store fastq files as array
files=(*.fastq.gz)

#perform fastqc -- distinction from for loop is this can process -t files simultaneously
fastqc "${files[@]}" -t 20 -o $out_fastqc

#summarize results with multiqc
multiqc $out_fastqc -o $out_fastqc


