### Mapping, deduplication, overlap clipping, and summary

* Reads are mapped to the *P. davidsonii* reference genome with bwa mem
* Duplicates marked with samtools markdup
* Overlapping paired end reads clipped with bamutil clipOverlap
* Some basic summary statistics generated with samtools coverage and samtools stats

All of the main commands are piped to avoid the creation of many intermediate files. See [`array_mapping_pipe.sh`](https://github.com/benstemon/dasanthera_novaseq/blob/main/mapping_etc/array_mapping_pipe.sh)

### Filter out reads with MQ <20 and that were marked as duplicates
* See [`array_filter_mapped_bams.sh`](https://github.com/benstemon/dasanthera_novaseq/blob/main/mapping_etc/array_filter_mapped_bams.sh)