# Testing peak calling by SEACR in bulk CUT&Tag data

I am here testing the optimal way to call peaks by SEACR on BAM files from PE reads.

SEACR is developed to call peaks from bedgraph files. Hence, after read mapping to the reference genome, the BAM file needs to be converted to bedgraph. This sounds like an easy task, however, while doing so, one should take into account of the SE or PE nature of the sequencing data.

I have the feeling that the way nf_core CUT&Run pipeline is developed does not take into account information from PE, but considers the alignemnts in the BAM file as they derive from SE and not PE data.

To investigate this I am comparing the method implemented in the nf_core CUT&Run pipeline (https://github.com/nf-core/cutandrun) with what is reported on SEACR GitHub page (https://github.com/FredHutch/SEACR).


## Rationale behind nf_core CUT&Run pipeline and custom code from SEACR GitHub
Both approaches convert the BAM file to bedgraph by ```bedtools genomecov```. However, while nf_core CUT&Run pipeline generates the bedgraph directly from the BAM file, SEACR GitHub suggests to first get the fragment coordinates from the BAM and then to convert such a file to bedgraph.

The main difference is that the approach reccomended by SEACR GitHub ensures that the bedgraph files reflect the density across read *pairs*, whereas by using nf_core CUT&Run pipeline the bedgraph files reflect the density of individual reads. 


## nf_core CUT&Run pipeline
This follows the code developed in the bedtools genomecov module of the nf_core CUT&Run pipeline (https://github.com/nf-core/cutandrun/blob/master/modules/nf-core/bedtools/genomecov/main.nf).

```

# Define args
sample=chr1
genome=Mus_musculus_GRCm39.chrom.sizes

# Convert BAM to bedgraph

# Call peaks

```
