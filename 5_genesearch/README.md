## Gene Search

### Generate focal gene sets for further analysis

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

#### Phylogenetic relationships in specific genomic regions
The first genomic region of interest is the first ~7 Mb on scaffold 2685. This region appears strongly discordant and harbors some interesting genes, including F3'H. It includes mRNA12393-mRNA12412

Estimate gene trees for concatenated CDS in the region:
1. Pull genes in the first 6MB of scaffold 2685 (f3pH_region_first6mb.bed)
2. use gffread to extract CDS for each individual -> use emboss to union fasta files
```
conda activate mapping_etc

for i in individual_fullgenome_fastas/*.fa;
do
    basename="${i##*/}"
    specname="${basename/consensus_fullgenome_bws_/}"
    
    gffread -x- -C -M -K -Y -E --sort-alpha -g $i REGION_f3ph_first6mb.bed | union -filter > f3ph_6mb/$specname.f3ph_6mb.fa
done

```
3. replace names of the fasta files and concatenate
```
cd f3ph_6mb/

for i in *.fa;
do
    sed -i "1s/.*/>$i/" $i
done

#concatenate and remove excess files
cat *.fa > combined_f3ph_6mb.fa
rm *.fa.f3ph*

#convert to single-line fasta
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' combined_f3ph_6mb.fa > nobreaks.combined_f3ph_6mb.fa
```
4. use IQtree to estimate tree from these files
```
iqtree -s nobreaks.combined_f3ph_6mb.fa --seqtype DNA -m MFP
```

5. cleanup
```
mv REGION_f3ph_first6mb.bed f3ph_6mb
mv f3ph_6mb /work/bs66/dasanthera_novaseq/genesearch
```



##### pulling useful information about gene trees from IQtree output
I was looking at the CDS gene trees from the 6mb region containing f3'h, to see if that would give me anything interesting. Most of these trees had very short branches and appears to have little information content. I wanted to quantify this, so I wrote a python script to summarize information content found in the .iqtree files.
* See [`parse_iqtree_files.py`](parse_iqtree_files.py) to pull information about number of nucleotides, invariant sites, and parsimony-informative sites for .iqtree files in a given directory. Usage:
```
python parse_iqtree_files.py --directory directory_with_.iqtree_files --outfile genetree_information.tsv
```

