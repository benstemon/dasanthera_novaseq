#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=pull_CDS_%j.out

###NOTE: this is set up as a batch array script. sbatch --array [0-17]

#set up directories
refgen="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
outdir="/work/bs66/dasanthera_novaseq/ABPgenes"
bamdir="/work/bs66/dasanthera_novaseq/mapped_filtered_bams"


#set up conda env
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc


cd $outdir

#set up sample array
bamlist=($bamdir/*.bam)
infile="${bamlist[$SLURM_ARRAY_TASK_ID]}"
insample=${infile%_mapped_filtered.bam}
insample=${insample##*/}


#set up information about the regions of interest
#"scaffold_2687:48556809-48558776:CHS1"
#"scaffold_2687:48867781-48869749:CHS2"
#"scaffold_2531:1356-3323:CHS3"
scafinfo=( "scaffold_1086:40697697-40699408:F3p5pH"
 "scaffold_2685:686427-691387:F3pH"
 "scaffold_2532:36924724-36926347:DFR"
 "scaffold_1086:38152303-38154886:FS-AFNS2"
 "scaffold_1086:12670236-12674403:FLS")


#extract the region
for i in ${scafinfo[@]};
do
    genename=$(awk -F: '{print $3}' <<< "$i")
    sinfo=$(awk -F: '{print $1 ":" $2}' <<< "$i")
    
    samtools view -b $infile $sinfo > $genename.$insample.bam
    samtools index $genename.$insample.bam
    
    #generate consensus sequence -- bamtools
    bcftools mpileup -f $refgen $genename.$insample.bam \
     --min-BQ 20 --min-MQ 20 -a FORMAT/AD,FORMAT/DP -Ou |
     bcftools call -m -Ou |
     bcftools filter --soft-filter LowQual --exclude '%QUAL<20' \
     --set-GTs . -o $genename.$insample.vcf
    
    #bgzip and index
    bgzip $genename.$insample.vcf
    tabix $genename.$insample.vcf.gz
    
    #generate consensus fasta
    samtools faidx $refgen $sinfo |
     bcftools consensus $genename.$insample.vcf.gz --haplotype I \
     --missing "N" --mark-del "-" --absent "-" -o $genename.$insample.tmpfasta
    
    #change the first line of the fasta to match species name
    sed -i "1s/.*/>$insample/" $genename.$insample.tmpfasta
    
    #change to single-line fasta and remove tmpfasta
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $genename.$insample.tmpfasta > $genename.$insample.fasta
    
    #cleanup and move to output directories
    rm $genename.$insample.tmpfasta
#    rm $genename.$insample.vcf.gz
#    rm $genename.$insample.vcf.gz.tbi
    mv $genename.$insample.* fastas_vcfs_bams
done

