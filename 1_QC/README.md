## NovaSeq data for 18 accessions of *Penstemon* subgenus *Dasanthera*

### Note: this pipeline requires the use of conda installs of various softwares. For info on navigating installation of these packages with conda on the cluster at South Carolina, see [this readme file](/conda_info/README.md)

### Quality control pipeline
1. Run QC of raw sequencing reads (fastqc) and summarize (multiqc)
2. Merge Illumina lanes by forward and reverse reads
3. 
    a. Trim adapters, quality filter, and enable base correction in overlapped regions (fastp)
    b. Run QC on trimmed, filtered data (fastqc)
    c. Summarize results (multiqc)



#### 1. QC on raw reads with fastqc, summarized with multiqc
* see [`QC_s1_rawqc-summarize.sh`](QC_s1_rawqc-summarize.sh)


#### 2. Merge data across lanes
* See [`QC_s2_merge-illumina-lanes.sh`](QC_s2_merge-illumina-lanes.sh)
* After using this script, move the merged reads to a new directory





#### 3a-c. Trimming and quality filtering (fastp), check QC on filtered reads (fastqc) and summarize (multiqc)
#### Trial 1:
Quality filtering with fastp, and check QC of trimmed reads
* Default options for quality filtering
* Base correction for overlapping reads enabled
* poly-x trimming on 3' ends enabled
* Unmerged reads that pass filter are included in separate files (default is to delete. reads failing to pass filter are not included)
This produces trimmed reads, that have been paired (in *trimmed*.fastq.gz) and unpaired (in *unpaired*.fastq.gz).
Moving forward with mapping, these unpaired reads can be mapped in the same way as paired reads, and then downstream BAMs can be merged with Picard's MergeSamFiles (or some other approach) 

**However, the initial trial highlighted some remaining issues, which were remediated in a second trial:**
* Still evidence of adapter contamination (didn't realize default auto-detection of adapters is off for PE data)
* Some reads probably too short (I now limit read length to 30 bp vs 15)
* The unpaired reads appeared to be low quality (I now do not keep unpaired reads)

**For the finalized version of steps 3a-c:**
* See [`QC_s3_fastp-filterqc-summarize.sh`](QC_s3_fastp-filterqc-summarize.sh)


### Alternative QC pipeline
All of the above steps can be performed in a single script, which is included.
* See [`QC_parentpops.sh`](QC_parentpops.sh)

