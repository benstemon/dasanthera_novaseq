#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=slurm_%j.out


#set up directories
refgen="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
outdir="/work/bs66/dasanthera_novaseq/ABPgenes/f3p5ph_specific"
bamdir="/work/bs66/dasanthera_novaseq/mapped_filtered_bams"


#set up variables
numthreads=1


#set up environments
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc
module load vcftools




# A. Filter bams to the region containing f3'5'h
#set up sample array
cd $outdir
bamlist=($bamdir/*.bam)

#generate filtered bams
for i in ${bamlist[@]};
do
    samplename=${i##*/}
    samtools view -b $i scaffold_1086:40697698-40699408 > $outdir/f3p5ph.$samplename
    samtools index $outdir/f3p5ph.$samplename
done




# B. call variants, including indels
#make new list of mapped, filtered bams
newbamlist=($outdir/*.bam)

#call variants, including indels, and normalize indels
bcftools mpileup --threads $numthreads -Ou \
 -a FORMAT/AD,FORMAT/DP \
 --min-BQ 20 --min-MQ 20 \
 -f $refgen ${newbamlist[@]} | 
 bcftools call --threads $numthreads -m --gvcf 2 -Ou |
 bcftools norm -f $refgen -Oz -o $outdir/unfiltered_vcf.gz

#index the vcf
tabix $outdir/unfiltered_vcf.gz




# C. filter the vcf, as done in previous iterations
# Invariant sites
 vcftools --gzvcf $outdir/unfiltered_vcf.gz \
 --max-non-ref-ac 0 \
 --recode --recode-INFO-all --stdout | 
 bcftools view - \
 --threads $numthreads \
 -Oz -o $outdir/tmp-invariants.filtered.vcf.gz

# Variant sites
vcftools --gzvcf $outdir/unfiltered_vcf.gz \
 --non-ref-ac 1 \
 --min-alleles 2 \
 --minGQ 20 \
 --min-meanDP 3 \
 --max-meanDP 60 \
 --minDP 2 \
 --max-missing 0.50 \
 --recode --recode-INFO-all --stdout | 
 bcftools filter - \
 --soft-filter LowQual \
 --exclude '%QUAL<20' \
 --set-GTs 0 \
 --threads $numthreads \
 -Oz -o $outdir/tmp-variants.filtered.vcf.gz

#index both vcfs
tabix $outdir/tmp-invariants.filtered.vcf.gz
tabix $outdir/tmp-variants.filtered.vcf.gz

#concatenate the two vcfs and remove tmp files
bcftools concat $outdir/tmp-invariants.filtered.vcf.gz \
 $outdir/tmp-variants.filtered.vcf.gz \
 --threads $numthreads --allow-overlaps -Oz \
 -o $outdir/filtered_consensus_ready.vcf.gz

#remove temporary files
rm $outdir/tmp-invariants.filtered.vcf.gz*
rm $outdir/tmp-variants.filtered.vcf.gz*

#index final merged file
tabix $outdir/filtered_consensus_ready.vcf.gz




# D. generate fasta files for each sample
# We need to pull out the exons individually, because the indels mess everything up
#exon 1: 40697698	40698069
#exon 2: 40698157	40698687
#exon 3: 40698776	40699408
for i in ${newbamlist[@]};
do
    samplename=${i##*/}
    samplename=${samplename/_mapped_filtered.bam/}
    samtools faidx $refgen scaffold_1086:40697698-40698069 |
     bcftools consensus --haplotype I --sample $i \
     --absent "N" --missing "N" --mark-del "-" \
     $outdir/filtered_consensus_ready.vcf.gz > $outdir/$samplename.exon1.fa

samtools faidx $refgen scaffold_1086:40698157-40698687 |
     bcftools consensus --haplotype I --sample $i \
     --absent "N" --missing "N" --mark-del "-" \
     $outdir/filtered_consensus_ready.vcf.gz > $outdir/$samplename.exon2.fa

samtools faidx $refgen scaffold_1086:40698776-40699408 |
     bcftools consensus --haplotype I --sample $i \
     --absent "N" --missing "N" --mark-del "-" \
     $outdir/filtered_consensus_ready.vcf.gz > $outdir/$samplename.exon3.fa

    #change the first line of each fasta to match the sample name
    sed -i "1s/.*/>$samplename exon 1/" $outdir/$samplename.exon1.fa
    sed -i "1s/.*/>$samplename exon 2/" $outdir/$samplename.exon2.fa
    sed -i "1s/.*/>$samplename exon 3/" $outdir/$samplename.exon3.fa

    #remove the first lines of the other files
    tail -n +2 $outdir/$samplename.exon2.fa > $outdir/$samplename.exon2.tmpfa
    tail -n +2 $outdir/$samplename.exon3.fa > $outdir/$samplename.exon3.tmpfa

    #concatenate exons together in single file and rename
    cat $outdir/$samplename.exon1.fa $outdir/$samplename.exon2.tmpfa $outdir/$samplename.exon3.tmpfa > tmp_catted.tmpfa
    sed -i "1s/.*/>$samplename full CDS/" tmp_catted.tmpfa

    #change to single-line fasta and remove tmpfasta
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $outdir/tmp_catted.tmpfa > $outdir/full_CDS_$samplename.fa
    rm $outdir/*.tmpfa
    
done




# E. Estimate gene tree
# The input fasta sequences are different lengths, because indels have been called.
# After we concatenate files, we must re-align them.
cat full_CDS*.fa > allseq_full_CDS_f3p5ph.fa
muscle -align allseq_full_CDS_f3p5ph.fa -output aligned_allseq_full_CDS_f3p5ph_indels.fa

#run iqtree with 100 UF bootstrap and the SH-like approximate likelihood ratio test ()
iqtree -s aligned_allseq_full_CDS_f3p5ph_indels.fa \
 --seqtype DNA \
 -m MFP \
 -T $numthreads \
 --prefix f3p5ph \
 -B 1000 \
 -alrt 1000



