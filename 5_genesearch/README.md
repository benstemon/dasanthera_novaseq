## Gene Search

### Generate focal gene sets foor further analysis

Workflow is to:
1. construct ABP genes search set (genes of interest from genbank)
2. create local blast database from davidsonii genome
3. blast the search set against the local database to identify the location of focal genes in the davidsonii genome
4. identify which CDS are within the bounds of these blast hits
5. pull those CDS locations from the resequenced genomes


#### 2. Create local blast database
```
refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"

makeblastdb -dbtype nucl -in $refgenome -out annot_Pdavidsonii_1mb
```


#### 3. Blast search set against local db
```
blastdb="/work/bs66/project_compare_genomes/davidsonii_blast_database"

#iteratively blast genes in the search set
for i in *.fasta;
do
    blastn -query $i -db $blastdb/annot_Pdavidsonii_1mb -outfmt 6 -evalue 1e-25 | awk -v OFS='\t' -v CDS="$i" '{$1=CDS; print}' >> blast_results/outresults.$i.txt
done
```


#### 4. Identify which CDS are within the bounds of blast hits

```
cd blast_results
CDSannot="/work/bs66/project_compare_genomes/annot_Pdavidsonii_1mb.gffread.genes.bed"

#pull 9th and 10th columns to produce .bed searchable file
#then use bedtools intersect to ID CDS within these bounds for each hit
#only do this for sequences with results files
for i in outresults*.txt;
do
    if test -s $i; then
    cat $i | awk '{if ($9 < $10) {print($2"\t"$9-1"\t"$10)} else {print($2"\t"$10-1"\t"$9)}}' | bedtools intersect -a stdin -b $CDSannot -wb > genomelocations.$i.bed
    fi
done
```



#### 5.1 generate file to sumarize this information

```
for i in genomelocations*;
do
    if test -s $i; then
    printname="${i/genomelocations.outresults./}";printname="${printname/.fasta.txt.bed/}"
    awk -v pn="$printname" -F'\t' '{ print $4, $5, $6, $7, pn }' $i | sort -k4 | uniq -f3 >> CDSlist_table.txt
    fi
done
```



#### 6. Pull these regions and generate new fastas
* See `pull_CDS_from_bam.sh`


#### One location had a blast hit but no annotated genes. Thoughts?

#### 







