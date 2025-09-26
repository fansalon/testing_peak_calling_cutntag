# testing_peak_calling_cutntag

I am here testing the optimal way to call peaks by SEACR on BAM files from PE reads.

SEACR is developed to call peaks from bedgraph files. Hence, after read mapping to the reference genome, the BAM file needs to be converted to bedgraph. This sounds like an easy task, however, while doing so, one should take into account of the SE or PE nature of the sequencing data.

I have the feeling that the way nf_core CUT&Run pipeline is developed does not take into account information from PE, but considers the alignemnts in the BAM file as they derive from SE and not PE data.

To investigate this I am comparing the method implemented in the nf_core CUT&Run pipeline (https://github.com/nf-core/cutandrun) with what is reported on SEACR GitHub page (https://github.com/FredHutch/SEACR).
