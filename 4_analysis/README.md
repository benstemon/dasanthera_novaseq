## Analysis

### Trees

#### Windowed gene trees
##### Prepare fasta files for windowed gene tree inference
The fasta files generated from consensus have newlines every 50bp. Additionally, they are sorted by species. We want to have a fasta file for each scaffold, with the sequence for that scaffold from each species. There is a script, [`TREES_1.rearrange_consensus_sequences.sh`](genetrees/TREES_1.rearrange_consensus_sequences.sh), which does this. Note that this functions on Linux OS, and will need modified slightly if using on MacOS (read script for more details). Script will need edited to match variable naming scheme. Use:

`bash TREES_1.rearrange_consensus_sequences.sh`


##### Generate windowed gene tree input files
This is done in two parts, with three files:
1. Make output directories for each scaffold. 
2. Generate gene trees with [`TREES_2.ARRAY_generate_windowed_genetree_infiles_masterscript.sh`](genetrees/TREES_2.ARRAY_generate_windowed_genetree_infiles_masterscript.sh), which uses the custom python script [`TREES.create_fasta_window_alignments.py`](TREES.create_fasta_window_alignments.py). This should function without modification. Parameters are defined in the batch script, including:
	* Input fasta file
	* Window size
	* Missing data threshold (only generates windows in which all species pass missing data threshold)
	* prefix to append to outfiles (can function as outdir)


##### Estimate windowed gene trees
This pipeline uses [IQtree](http://www.iqtree.org/) for gene tree inference. Again, this is done in two parts: 
1. Make output directories for gene trees estimated along each scaffold.
2. Estimate gene trees (in array batch submission), setting outgroup (here it is P. montanus), and specifying substitution model inference and out prefix.
* See [`TREES_3.ARRAY_estimate_windowed_genetrees_iqtree.sh`](genetrees/TREES_3.ARRAY_estimate_genetrees_iqtree.sh)



#### CDS gene trees

##### Prepare CDS fasta files for gene tree inference
The fasta files generated for CDS regions are initially sorted by sample. To estimate CDS gene trees, we need a fasta for each CDS, with each sample's sequence included. First we will want to reformat our fasta files, so they are single line fasta rather than interleaved. Use this code on the CDS fastas.
```shell

for i in *.fa;
do
    awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < $i > "${i/.fa/.fixed.fa}"
    rm $i
done

```

Next, we will want to generate the infiles for gene tree inference. To generate these, see [`CDS_TREES.concat_CDS_fastas.py`](genetrees/CDS_TREES.concat_CDS_fastas.py). This python script takes input fasta, output directory, and missing data threshold, and appends to an output fasta for each CDS, naming the output after the scaffold and region the CDS corresponds to. A simple shell for loop can be run with the python script to add all samples to the output. The missing data threshold here filters for individuals, rather than windows (i.e., if a sample has more missing data than desired, that individual is not added to the output fasta, rather than the output fasta not being generated).

```shell
for i in individual_CDS_fastas/CDS_*.fa;
do
    python3 CDS_TREES.concat_CDS_fastas.py -i $i -o /work/bs66/dasanthera_novaseq/analysis/CDS_genetree_infiles -m 0.5
done
```

Here I generate aligned fasta for each CDS, excluding individuals with > 50% missing data. The script is fast, and possible to use without batch submission. But for many samples or many CDS regions, it would be wise to write a batch script so as not to overload the head node.


I also generated individual fasta files for each CDS for the reference genome. This will be useful for later inferences when looking at function of specific genomic regions.
```
#give the 1mb genome and desired output directory

refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
bedfile="/work/bs66/project_compare_genomes/annot_Pdavidsonii_1mb.gffread.genes.bed"
outdir="/work/bs66/project_compare_genomes/davidsonii_CDS"


gffread -x $outdir/concat_CDS_davidsonii_refgenome.fa -C -M -K -Y -E --sort-alpha -g $refgenome $bedfile


#change to single line per fasta
for i in *.fa;
do
    awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < $i > "${i/.fa/.fixed.fa}"
    rm $i
done


#use python script to generate individual fastas for each CDS
python3 CDS_TREES.concat_CDS_fastas.py -i concat_CDS_davidsonii_refgenome.fixed.fa -o $outdir

#remove files no longer needed
#rm CDS_TREES.concat_CDS_fastas.py
#rm concat_CDS_davidsonii_refgenome.fa
```


##### Estimate CDS gene trees
Still using IQtree for gene tree inference here. This used to be done on a scaffold-by-scaffold basis, but because that information has been lost in the header, instead we are just submitting all files for inference. They could be split up into sections to speed this part up.
* To infer gene trees, see [`CDS_TREES_1.estimate_CDS_genetrees_iqtree.sh`](genetrees/CDS_TREES_1.estimate_CDS_genetrees_iqtree.sh)




#### ASTRAL species tree: windowed analysis

First, cat trees into one large file. Change to directory with subdirectories, each containing gene trees estimated for each scaffold.

```shell
cd /work/bs66/dasanthera_novaseq/analysis/genetree_outfiles
for i in scaf_*;
do
    cat $i/*.treefile >> combined_10kbwindowtrees.tre
done
```

Also, keep a log of which trees were catted to the treefile, and in what order. This will come in handy later. It is a good idea to verify that the order is the same, but if done this way it should be correct.

```shell
for i in scaf_*;
do
    echo $i/*.treefile | tr " " "\n" >> tmpout.txt
done

cat --number tmpout.txt > numbered_10kbtreepaths.txt
rm tmpout.txt
```

I put both of these outputs in a new directory, `treemetrics`. They will be used for several downstream analyses. Next, use the combined windowtrees file to estimate the species tree in ASTRAL.

* See [`SPECIESTREE_1.astral_windowtrees.sh`](speciestrees/SPECIESTREE_1.astral_windowtrees.sh)



#### ASTRAL species tree: CDS analysis

Perform the same commands as described for the windowed analysis; cat trees to a combined tree file and keep a log of the order in which this was done. 

```shell
cd /work/bs66/dasanthera_novaseq/analysis/CDS_genetree_outfiles
for i in *.treefile;
do
    cat $i >> combined_CDStrees.tre
done

for i in *.treefile;
do
    echo $i | tr " " "\n" >> tmpout.txt
done

cat --number tmpout.txt > numbered_CDStreepaths.txt
rm tmpout.txt

mv combined_CDStrees.tre ../treemetrics
mv numbered_CDStreepaths.txt ../treemetrics
```

Then, estimate the species tree in ASTRAL.
* See [`SPECIESTREE_2.astral_CDS.sh`](speciestrees/SPECIESTREE_2.astral_CDS.sh)




#### Concatenated ML species tree
The concatenated ML species tree is estimated in IQtree. I specified the GTR+I+R model, and estimated rates for each scaffold, performing 1000 ultra-fast bootstrap replicates, and specifying _P. montanus_ as the outgroup.
* See [`SPECIESTREE_3.IQtree_concat.sh`](speciestrees/SPECIESTREE_3.IQtree_concat.sh)


#### Tree metrics
##### Robinson-Foulds distance
Plotting normalized RF distance of 10kb sliding window trees compared to the species tree. This was plotted along with genic content to visualize whether there is a relationship between RF distance and genic content along scaffolds. For the script to calculate RF distance and plot results, see [`plot_RF_x_genecontent.R`](treemetrics/plot_RF_x_genecontent.R)


#### Concordance factors
Estimated site concordance factors and gene concordance factors in IQtree. See [`TREEMETRICS_1.calculate_concordance_factors.sh`](treemetrics/TREEMETRICS_1.calculate_concordance_factors.sh). Prior to running for CDS, needed to filter out sites with fewer than 3 samples, using python script [`remove_fewsequence_CDS.py`](remove_fewsequence_CDS.py)

From this, used the .log file to match locus IDs (10kb trees) to their site concordance factors score. First, I copied the output from the log file containing identification into a new .txt file, then used `grep -v '^WARNING' input.txt > 10kbsource_cf_site_ids.txt` to generate an output text file with only this information.

* This information, was plotted (in same script as plots for TWISST internal branch plots) with [`plot_twisst_IB-test_and_CFs_IB.R`](twisst_internal_branches/plot_twisst_IB-test_and_CFs_IB.R)





### D-statistics
We used Dsuite to calculate all introgression statistics, including f-branch tests and fdM. These were performed in three main groups of analyses.

#### Dstats-fbranch (whole-genome signatures of introgression)
First, perform Dsuite Dtrios to generate D statistics for all possible triplets of taxa. These tests do not include P. lyallii and use the Astral 10kb topology as the input species tree. Samples here were treated as individuals.
* See [`DSTATS_1.ARRAY_dsuite_dtrios.sh`](Dstats_fbranch/DSTATS_1.ARRAY_dsuite_dtrios.sh) for this first batch script, as well as [`intree_dtrios.tre`](Dstats_fbranch/intree_dtrios.tre) for the input species tree and [`popset_dtrios.txt`](Dstats_fbranch/popset_dtrios.txt) for the population identification set.

Next, combine Dtrios output for each scaffold into a genome-wide analysis. This will also produce an f-branch plot for the whole genome.
* See [`DSTATS_2.combineDtrios.sh`](Dstats_fbranch/DSTATS_2.combineDtrios.sh)



#### Dstats-fullspecies-sliding (signatures of introgression across genomic regions)
We used Dsuite Dinvestigate to calculate introgression metrics in sliding windows across the genome. These scripts estimate D, fd, fdM, and Df in overlapping windows of 10kb SNPs, sliding every 2500 SNPs. Rather than as individuals, samples were assigned to species identity, including dav118 as P. fruticosus. P. lyallii was not included in these analyses, and P. montanus was used as the outgroup.

We were also interested in the potential relationship that introgression metrics have with gene density. To explore this relationship, we calculated the proportion of coding sites for each of the 10kb SNP windows using a custom script and bedtools. See script [`generate_dstat_genic_plotfiles.sh`](miscellaneous/generate_dstat_genic_plotfiles.sh) for the script. It is used as follows:
```bash
# Input for the script includes:
# -d : the d-statistics outfile produced from Dinvestigate
# -c : a .bed file with CDS coordinates
# -o : the name for the tab-separated output file


chmod +x generate_dstat_genic_plotfiles.sh

./generate_dstat_genic_plotfiles.sh -d dstat_file.txt -c CDS_regions.bed -o outfile.bed

# the script will generate a "plotting" file which can be used to compare Dinvestigate output with genic fraction.
```

Two versions of these analyses were conducted: one including all samples (except P. lyallii), and one also excluding fru106 (to see the effect removing this individual had on the relationship with genic fraction). We calculated these statistics for every possible rooted triplet of taxa (ten total rooted triplets), specifying relationships as inferred by the species tree.
* For the full sample analysis, see [`DSTATS_fullspecies_1a.sliding_10kb.sh`](Dstats_fullspecies_sliding/DSTATS_fullspecies_1a.sliding_10kb.sh) for the batch script and [`popset_fullspecies.txt`](Dstats_fullspecies_sliding/popset_fullspecies.txt) for the population set.
* For the analysis excluding fru106, see [`DSTATS_fullspecies_1b_nofru106.sliding_10kb.sh`](Dstats_fullspecies_sliding/DSTATS_fullspecies_1b_nofru106.sliding_10kb.sh) for the batch script and [`popset_fullspecies_nofru106.txt`](Dstats_fullspecies_sliding/popset_fullspecies_nofru106.txt) for the population set.
* Both analyses also implement the trioset [`trioset_fullspecies.txt`](Dstats_fullspecies_sliding/trioset_fullspecies.txt) for specifying each of the rooted triplets, and [`plot_Dstats_fullspecies_sliding_genic.R`](Dstats_fullspecies_sliding/plot_Dstats_fullspecies_sliding_genic.R) for plotting.


#### Dstats-introtests (introgression across the genome for targeted rooted triplets)
We again used Dsuite Dinvestigate to estimate D, Df, fD, and  fdM. This time, however, we generated tests to mirror those performed in the twisst-popspecific analyses. We thus have four popsets and four triosets in [`Dstats_introtests`](Dstats_introtests) that correspond to the four focal TWISST tests: (1) dav-new, (2) dav116, (3) rup86, and (4) rup101. We estimated these D statistics at four different scales:
* 10k SNP windows, sliding every 2500 SNPs [`DSTATS_1.ARRAY_introtests.sh`](Dstats_introtests/DSTATS_1.ARRAY_introtests.sh)
* 5k SNP windows, sliding every 1250 SNPs [`DSTATS_1.ARRAY_introtests_5kb.sh`](Dstats_introtests/DSTATS_1.ARRAY_introtests_5kb.sh)
* 1k SNP windows, sliding every 250 SNPs [`DSTATS_1.ARRAY_introtests_1kb.sh`](Dstats_introtests/DSTATS_1.ARRAY_introtests_1kb.sh)
* 500 SNP windows, sliding every 125 SNPs [`DSTATS_1.ARRAY_introtests_500bp.sh`](Dstats_introtests/DSTATS_1.ARRAY_introtests_500bp.sh)

We also calculated the genic fraction for each of the windows produced by these analyses, as described above (ACTUALLY, NOT YET). Finally, plots were made with [`plot_Dstats_introtests.R`](Dstats_introtests/plot_Dstats_introtests.R)





### Gene identities in significant outliers from D-window analyses

Using output from sliding window D analyses, identify the CDS in outlier windows and their function. First, use bedtools intersect to identify maker CDS in these windows.

```shell
CDSannot="/work/bs66/project_compare_genomes/annot_Pdavidsonii_1mb.gffread.genes.bed"
fullannot="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.gff"

#bedtools command
bedtools intersect -a significant_windows_1000_500_z4.bed -b $CDSannot -wb > CDS-hits_only_1000_500_z4.bed

```

Next we will want to make a local blastx database against which we will blast the CDS. To do this we use the swissprot database. Set this up like so:

```shell
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
gunzip swissprot.gz

module load blast
makeblastdb -dbtype prot -in swissprot

mkdir db_Swissprot
mv swissprot* db_Swissprot/
```

Finally, see [`Dwindow_BLASTX_CDS.sh`](Dwindow_outlier_analysis/Dwindow_BLASTX_CDS.sh) to run the blast search, and [`explore_CDS.R`](Dwindow_outlier_analysis/explore_CDS.R) to filter results and generate final tables of CDS function.



### pixy fullgenome
#### fullspecies (10kb and 50kb)
For this analysis, we estimated pi, dxy, and fst genome-wide for all cross-species comparisons. These were estimated in 10kb and 50kb stretches. Species were assigned to their taxonomic affinity. "dav118" was assigned to P. fruticosus.

The goal of these analyses are to visualize how these diversity and differentiation metrics vary across the genome, and whether there are shared patterns across species and/or with particular genomic features.

* See [`PIXY_fullgenome_fullspecies_1.10kb.sh`](pixy/pixy_fullgenome_fullspecies_10kb) and [`PIXY_fullgenome_fullspecies_1.50kb.sh`](pixy/pixy_fullgenome_fullspecies_50kb) for batch scripts used. popfiles for each analysis are within their respective directories.
* See [`plot_pixy_fullgenome_fullspecies.R`](pixy/plot_pixy_fullgenome_fullspecies.R) for plotting. Several plots were made, including:
	- genic fraction vs. average dxy, pi, and fst
	- average pi, dxy, and fst plotted along scaffolds for each species/interaction





### Genomic PCA
Using the filtered .vcf file, generate a genome-wide PCA.
* See script [`generate_biallelic-filtered_vcf.sh`](PCA/generate_biallelic-filtered_vcf.sh)

Then use following code in plink to generate PCA
```
#generate files for PCA
plink --vcf PCA-ready_biallelic_ld-filtered.vcf --double-id --allow-extra-chr \
 --set-missing-var-ids @:# --make-bed --pca --out PCA
```
I noticed that this produces eigenvalues that don't look right (negative values). So, I performed a second PCA using a different method to compare to the plink version.

* See [`plot_plinkPCA_and_generate_adegenetPCA.R`](PCA/plot_plinkPCA_and_generate_adegenetPCA.R) For the R script which plots the plink output and has a second section to perform PCA in adegenet.






### TWISST
We ran TWISST using the 10kb sliding window trees as input. There are three main TWISST modules that we generated. 

#### TWISST fullspecies
This module employed all samples except for P. lyallii. We ran twisst for every possible combination of three taxa (rooted with P. montanus -- 10 total tests). We then smoothed topology weights in 2Mb windows, spacing every 5kb, to produce smoothed topology weights. We then summarized TWISST results by creating two new metrics based on the smoothed topology weights, where topo1 is always the topology concordant with the species tree: (1) total discordance, which is the topology weight for the two topologies discordant from the species tree topology (w(topo2) + w(topo3)), and (2) discordance imbalance, which is the difference in topology weights of the two discordant topologies (abs(w(topo2)-w(topo3))). These were plotted in combination with genic fraction information.
* See [`TWISST_1.run_twisst_fullspecies.sh`](twisst_fullspecies/TWISST_1.run_twisst_fullspecies.sh) for the batch script and [`twisst_fullspecies_groupsfile.txt`](twisst_fullspecies/twisst_fullspecies_groupsfile.txt) for the groupsfile as input for twisst fullspecies analyses.

* See [`plot_twisst_fullspecies.R`](twisst_fullspecies/plot_twisst_fullspecies.R) for plotting script.

#### TWISST popspecific
This module examined more intensively some of the focal introgressed individuals identified from the f-branch statistics. Each focal introgressed individual (P2) has an introgression partner (P3), and is part of a broader species identification (P2.species). For each P2 + P3 combination targeted, we ran two tests of rooted triplets: (1) P2, P3, and sister taxon (sister to either P2 or P3), and (2) P2.species (non-introgressed), P3, and sister taxon. This way, differences in topology weights between test 1 and 2 can be attributed directly to the evolutionary history unique to P2. These were then plotted using the functions made available by the TWISST developers, smoothing values over 2Mb in 5kb increments.
* See [`TWISST_1.run_twisst_popspecific_v2.sh`](twisst_popspecific_v2/TWISST_1.run_twisst_popspecific_v2.sh) for the batch script, and [`plot_twisst_popspecific_v2.R`](twisst_popspecific_v2/plot_twisst_popspecific_v2.R) for the plotting script.


#### TWISST internal branch test
This last module considered all possible tree topologies consistent wth the three main internal branches of the species tree: 1: new+car+rup, 2: new+car, and 3: dav+fru. 
* See [`groupsfile_IB-test.txt`](twisst_internal_branches/groupsfile_IB-test.txt) for the groupsfile, and [`TWISST_1.run_twisst_IB-test.sh`](twisst_internal_branches/TWISST_1.run_twisst_IB-test.sh) for the batch script to run these analyses.
* We plotted the topology weights for each 5 taxon tree consistent with each of these three internal branches with [`plot_twisst_IB-test_and_CFs_IB.R`](twisst_internal_branches/plot_twisst_IB-test_and_CFs_IB.R). Note that this R script also includes code for plotting concordance factors for these same internal branches.


### Miscellaneous
This also falls under the Dstats category, but see the script [`generate_dstat_genic_plotfiles.sh`](miscellaneous/generate_dstat_genic_plotfiles.sh) to turn the output from Dinvestigate into a plottable format that includes information about the genic content in that region of the genome. Use of this script requires a the Dinvestigate output file and a .bed file with CDS coordinates (or other genomic intervals of interest). Use of the script is as follows:

```shell
chmod +x generate_dstat_genic_plotfiles.sh

./generate_dstat_genic_plotfiles.sh -d dstat_file.txt -c CDS_regions.bed -o outfile.bed
```



To obtain information about gene content in sliding windows, see [`calculate_percentage_CDS.sh`](miscellaneous/calculate_percentage_CDS.sh). Using this script will require a .txt file with the names and sizes of each chromosome (see [`genomesize_scaffolds_davidsonii_1mb.txt`](miscellaneous/genomesize_scaffolds_davidsonii_1mb.txt)) and a .bed file with the genomic coordinates of CDS.
* To enable use of the script as intended, it must be made executable with `chmod +x calculate_percentage_CDS.sh` and run with `./calculate_percentage_CDS.sh`




### Gene Density plots
* plotted gene density vs. RF distance
* plotted gene density vs. pixy
