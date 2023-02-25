In this readme is information pertaining to analyses that were performed but ultimately not used. Saving temporarily in case they are needed again, and for reference.

### pixy CDS (unused)
Note: the pixy information here is a bit disorganized because I ended up doing a lot of different analyses, many of which were not ultimately used to produce anything meaningful. I will need to clean up this information and the scripts within to make clear what actually happened.

#### pixy CDS triplets (unused)
We obtained unbiased estimates of Dxy from pixy between species in the NCR clade (newberryi, cardwellii, rupicola) to assess how inferred relationships differ depending on the gene of interest. First, we compared Dxy between P. rupicola, P. newberryi, and P. cardwellii in all CDS regions. CDS were required to have > 100 sites and missing count/count comparisons was required to be < 50%. We generated plots for scenarios grouping all individuals by species (fullspecies) as well as specific individual-by-individual comparisons. As expected, Dxy in CDS regions between P. newberryi and P. cardwellii were on average lower than average Dxy between P. newberryi and P. rupicola, and between P. rupicola and P. cardwellii. This pattern held whether we grouped individuals by species (fullspecies) or whether we made individual-by-individual comparisons (popspecific).
* For pixy analyses, see [`PIXY_1a.CDS_NCR_fullspecies.sh`](pixy/pixy_CDS/PIXY_1a.CDS_NCR_fullspecies.sh) and [`popfile_1a.NCR_fullspecies.txt`](pixy/pixy_CDS/popfile_1a.NCR_fullspecies.txt)
* For plotting, see [`plot_pixy_CDS_dxy.R`](pixy/pixy_CDS/plot_pixy_CDS_dxy.R)

#### pixy CDS allpops (unused)
Obtained estimates of dxy from all possible samples, but in CDS regions.
* See [`PIXY_allpops_1.CDS_NCR_popspecific.sh`](pixy/pixy_CDS_allpops/PIXY_allpops_1.CDS_NCR_popspecific.sh) for the bash script and [`popfile_allpops.NCR_popspecific.txt`](pixy/pixy_CDS_allpops/popfile_allpops.NCR_popspecific.txt) for the popfile.


### geodist_gendist (unused)

We examined the relationship between genome-wide dxy and geographic distance amongst all possible population comparisons. Plots were made on a per-species basis. We also used these dxy comparisons to make MDS plots (multidimensional scaling).
* See [`geodist_gendist_MDS_plotting.R`](geodist_gendist_MDS/geodist_gendist_MDS_plotting.R) for the R script to filter and plot data, and [`dasanthera_coords.txt`](geodist_gendist_MDS/dasanthera_coords.txt) for coordinate data for each population.


### Old Dinvestigate scripts
Finally, investigate targeted triplets of interest for sliding window introgression metrics. I tested a few different window sizes. As an example:
* For conducting this analysis, see [`DSTATS_3a.dsuite_Dinvestigate_1000_500.sh.sh`](Dstats/DSTATS_3a.dsuite_Dinvestigate_1000_500.sh). Will also need [`popset_dtrios.txt`](Dstats/popset_dtrios.txt)
* For plotting output of this analysis and generating outfile for outlier windows, see, for example: [`plot_Dinvestigate_Dwindow_1000_500.R`](Dstats/plot_Dinvestigate_Dwindow_1000_500.R)

