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
Using the conda install version of ete3 here. I am having trouble getting this to work through batch submission, so used an interactive node instead.

```
idev

#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#activate conda environment with packages installed
#need ete3 installed
conda activate quibl_etc


#specify reference (species) tree, and the combined window tree file
#reference tree should not have annotations so I made a new, cleaned up treefile
reftree="/work/bs66/dasanthera_novaseq/analysis/astral_trees/astral_10kb_noannotations.tre"
treelist="/work/bs66/dasanthera_novaseq/analysis/treemetrics/combined_10kbwindowtrees.tre"
outdir="/work/bs66/dasanthera_novaseq/analysis/treemetrics"


#run ete3
ete3 compare --src_tree_list $treelist -r $reftree --unrooted --taboutput > $outdir/RFdistance_10kbwindow_astralref.txt
```

##### Finding discordant topologies given specified triplets
- conda activate quibl_etc

- manually remove all annotations and branch length information -- just use topology of CDS astral species tree

Use R script to determine whether pruned triplets are concordant, discordant, etc.





### D-statistics and related metrics
* Using Dsuite to calculate various metrics related to introgression
* Using an input tree to additionally test for significance of D with respect to tree topology. Tree for input here is the MLE tree for concatenated species tree, rooted with P. montanus. P. lyalli has been removed from the tree, since we only need one outgroup and it complicates things slightly to keep it.


First. perform Dsuite dtrios to generate D statistics for all possible triplets. Include tree topology to generate tree.txt files, which are D metrics with triplets arranged as they are in the tree. This is done on a scaffold-by-scaffold basis.
* See [`DSTATS_1.ARRAY_dsuite_dtrios.sh`](Dstats/DSTATS_1.ARRAY_dsuite_dtrios.sh)

Next, combine Dtrios output for each scaffold into a genome-wide analysis. These first two scripts will also generate plots of significant f-branch results.
* See [`DSTATS_2.combineDtrios.sh`](Dstats/DSTATS_2.combineDtrios.sh). Will also need [`intree_NOLYALLII_dtrios.tre`](Dstats/intree_NOLYALLII_dtrios.tre) and [`popset_dtrios.txt`](Dstats/popset_dtrios.txt)

Finally, investigate targeted triplets of interest for sliding window introgression metrics. I tested a few different window sizes. As an example:
* For conducting this analysis, see [`DSTATS_3a.dsuite_Dinvestigate_1000_500.sh.sh`](Dstats/DSTATS_3a.dsuite_Dinvestigate_1000_500.sh). Will also need [`popset_dtrios.txt`](Dstats/popset_dtrios.txt)
* For plotting output of this analysis and generating outfile for outlier windows, see, for example: [`plot_Dinvestigate_Dwindow_1000_500.R`](Dstats/plot_Dinvestigate_Dwindow_1000_500.R)




#### Gene identities in significant outliers from D-window analyses

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






### pixy CDS

#### pixy NCR triplets
We obtained unbiased estimates of Dxy from pixy between species in the NCR clade (newberryi, cardwellii, rupicola) to assess how inferred relationships differ depending on the gene of interest. First, we compared Dxy between P. rupicola, P. newberryi, and P. cardwellii in all CDS regions. CDS were required to have > 100 sites and missing count/count comparisons was required to be < 50%. We generated plots for scenarios grouping all individuals by species (fullspecies) as well as specific individual-by-individual comparisons. As expected, Dxy in CDS regions between P. newberryi and P. cardwellii were on average lower than average Dxy between P. newberryi and P. rupicola, and between P. rupicola and P. cardwellii. This pattern held whether we grouped individuals by species (fullspecies) or whether we made individual-by-individual comparisons (popspecific).
* For pixy analyses, see [`PIXY_1a.CDS_NCR_fullspecies.sh`](pixy/pixy_CDS/PIXY_1a.CDS_NCR_fullspecies.sh) and [`popfile_1a.NCR_fullspecies.txt`](pixy/pixy_CDS/popfile_1a.NCR_fullspecies.txt)
* For plotting, see [`plot_pixy_CDS_dxy.R`](pixy/pixy_CDS/plot_pixy_CDS_dxy.R)

#### pixy allpops
We also obtained estimates of Dxy from all possible samples. This data set is redundant with pixy NCR, so at some point I should reconcile that.
* See [`PIXY_allpops_1.CDS_NCR_popspecific.sh`](pixy/pixy_CDS_allpops/PIXY_allpops_1.CDS_NCR_popspecific.sh) for the bash script and [`popfile_allpops.NCR_popspecific.txt`](pixy/pixy_CDS_allpops/popfile_allpops.NCR_popspecific.txt) for the popfile.



### geodist_gendist

We examined the relationship between genome-wide dxy and geographic distance amongst all possible population comparisons. Plots were made on a per-species basis. We also used these dxy comparisons to make MDS plots (multidimensional scaling).
* See [`geodist_gendist_MDS_plotting.R`](geodist_gendist_MDS/geodist_gendist_MDS_plotting.R) for the R script to filter and plot data, and [`dasanthera_coords.txt`](geodist_gendist_MDS/dasanthera_coords.txt) for coordinate data for each population.







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
Using the 10kb gene trees, we ran TWISST. We focused on the NCR and NCRD clades for simplification. Plots were made by smoothing weightings in 2Mb increments.
* See [`TWISST_1.run_twisst.sh`](twisst/TWISST_1.run_twisst.sh) for the shell script to run TWISST. This requires [`twisst_groupsfile.txt`](twisst/twisst_groupsfile.txt) to function.
* See [`plot_twisst_output.R`] for plotting TWISST output.

