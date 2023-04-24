#!/bin/sh

#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p wessinger-48core
#SBATCH --job-name=QC_parentpops

cd $SLURM_SUBMIT_DIR


#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with QC packages installed
conda activate QC

#Change these: location of raw reads and location of the output directory
#(will need to make output directory)
basedir="/work/bs66/dasanthera_novaseq/"
rawreads="/work/bs66/dasanthera_novaseq/rawreads_parentpops"

#############################
#fastqc+multiqc -- raw reads#
#############################

#change to wd and make list of files as array
cd $rawreads
files=(*.fastq.gz)

#run fastqc and summarize with multiqc
mkdir ${basedir}fastqc_rawreads_parentpops
fastqc "${files[@]}" -t 20 -o ${basedir}fastqc_rawreads_parentpops
multiqc ${basedir}fastqc_rawreads_parentpops -o ${basedir}fastqc_rawreads_parentpops


#############################
### merge illumina lanes ###
#############################
cd $rawreads
for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

    do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_merged_L001_R1_001.fastq.gz

       echo "Merging R2"

cat "$i"_L00*_R2_001.fastq.gz > "$i"_merged_L001_R2_001.fastq.gz

done;

#move merged reads to new merged read directory
mkdir ${basedir}merged_reads_parentpops
mv *merged* ${basedir}merged_reads_parentpops


#############################
########### fastp ###########
#############################

#make new filtered output directory
mkdir ${basedir}filtered_reads_parentpops


#fastp for loop
cd ${basedir}merged_reads_parentpops
for r1in in *_R1_001.fastq.gz; 
do
    r2in="${r1in/R1_001.fastq.gz/R2_001.fastq.gz}"
    r1out="${r1in##*/}"
    r2out="${r1out/R1_001.fastq.gz/R2_001.fastq.gz}"
    fastp -i "$r1in" -I "$r2in" --detect_adapter_for_pe -l 30 --out1 ${basedir}filtered_reads_parentpops/"${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R1_001.fastq.gz}" --out2 ${basedir}filtered_reads_parentpops/"${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R2_001.fastq.gz}" -x -c -w 16 -h ${basedir}filtered_reads_parentpops/"${r1out/merged_L001_R1_001.fastq.gz/html}" -j ${basedir}filtered_reads_parentpops/"${r1out/merged_L001_R1_001.fastq.gz/json}"
done



#############################
#fastqc+multiqc -- filtered #
#############################

#change to wd and make list of files as array
cd ${basedir}filtered_reads_parentpops
filterfiles=(*.fastq.gz)

#run fastqc and summarize with multiqc
mkdir ${basedir}fastqc_filtered_parentpops
fastqc "${filterfiles[@]}" -t 20 -o ${basedir}fastqc_filtered_parentpops
multiqc ${basedir}fastqc_filtered_parentpops -o ${basedir}fastqc_filtered_parentpops




