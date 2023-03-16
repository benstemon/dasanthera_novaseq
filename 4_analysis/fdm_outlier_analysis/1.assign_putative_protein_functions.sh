#!/bin/sh

#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p wessinger-48core
#SBATCH --job-name=pblast_functional
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.error
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bs66@mailbox.sc.edu



#goal of the script is to add putative functions of gene models to .gff file
#run the script in the desired output directory

#load correct environments
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#load blast
module load blast


#path to maker -- add to path
makerpath="/work/bs66/software/maker/bin"
export PATH=$PATH:$makerpath


#specify blast db (uniprot database set up previously)
blastdb="/work/bs66/project_compare_genomes/blast_database/db_uniprot_sprot.fasta/uniprot_sprot.fasta"

#specify the protein fasta file (all annotated proteins), gff file, and the path to maker
proteinfasta="/work/bs66/project_compare_genomes/davidsonii_unfiltered/FUNC-ANNO_davidsonii_protein.fasta"
gff_file="/work/bs66/project_compare_genomes/davidsonii_unfiltered/annot_Pdavidsonii_genome.gff"




#run a protein blast for each of the annotated proteins. I would normally have this information from the original set of annotations, but since I don't, I need to re-do it.
blastp -num_threads 24 -query $proteinfasta -db $blastdb -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out output.blastp



#use maker_functional_gff to add functions from BLAST to GFF3
maker_functional_gff $blastdb output.blastp $gff_file > annot_Pdavidsonii_genome_putative_function.gff

#use maker_functional_fasta to add functions from BLAST to  FASTA
maker_functional_fasta $blastdb output.blastp $proteinfasta > FUNC-ANNO_davidsonii_protein_putative_function.fasta


