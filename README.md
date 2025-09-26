# Testing peak calling by SEACR in bulk CUT&Tag data

I am here testing the optimal way to call peaks by SEACR on BAM files from PE reads.

SEACR is developed to call peaks from bedgraph files. Hence, after read mapping to the reference genome, the BAM file needs to be converted to bedgraph. This sounds like an easy task, however, while doing so, one should take into account of the SE or PE nature of the sequencing data.

I have the feeling that the way nf_core CUT&Run pipeline is developed does not take into account information from PE, but considers the alignemnts in the BAM file as they derive from SE and not PE data.

To investigate this I am comparing the method implemented in the nf_core CUT&Run pipeline (https://github.com/nf-core/cutandrun) with what is reported on SEACR GitHub page (https://github.com/FredHutch/SEACR).


## Rationale behind nf_core CUT&Run pipeline and custom code from SEACR GitHub
Both approaches convert the BAM file to bedgraph by ```bedtools genomecov```. However, while nf_core CUT&Run pipeline generates the bedgraph directly from the BAM file, SEACR GitHub suggests to first get the fragment coordinates from the BAM and then to convert such a file to bedgraph.

The main difference is that the approach reccomended by SEACR GitHub ensures that the bedgraph files reflect the density across read *pairs*, whereas by using nf_core CUT&Run pipeline the bedgraph files reflect the density of individual reads. 

Before starting, get the repository and get into the directory:
```
git clone https://github.com/fansalon/testing_peak_calling_cutntag.git
cd testing_peak_calling_cutntag
```


## nf_core CUT&Run pipeline
This follows the code developed in the bedtools genomecov module of the nf_core CUT&Run pipeline (https://github.com/nf-core/cutandrun/blob/master/modules/nf-core/bedtools/genomecov/main.nf).

```
# Define args
sample=chr19

# Prepare directory
mkdir nfcore_approach
cd nfcore_approach

# Convert BAM to bedgraph
bedtools genomecov -bga -ibam ../input/${sample}.bam > ${sample}.bedgraph

# Call peaks
SEACR_1.3.sh ${sample}.bedgraph 0.01 norm stringent seacr_out

cd ..
```


## SEACR GitHub
This follows the code described in SEACR GitHub (https://github.com/FredHutch/SEACR?tab=readme-ov-file#preparing-input-bedgraph-files)

```
# Define args
sample=chr19
genome=Mus_musculus_GRCm39.chrom.sizes

# Prepare directory
mkdir seacr_approach
cd seacr_approach

# sort by query	
samtools sort -@ 10 -n -o ${sample}_sorted.bam ../input/${sample}.bam

# convert to bed
bedtools bamtobed -bedpe -i ${sample}_sorted.bam > ${sample}.bed

# extract fragments shorter than 1kb
cat ${sample}.bed | awk '$1==$4 && $6-$2 < 1000 {print $0}' > $sample.clean.bed
cut -f 1,2,6 $sample.clean.bed | sort -k1,1 -k2,2n -k3,3n > $sample.fragments.bed

# convert to bedgraph
bedtools genomecov -bg -i $sample.fragments.bed -g ../input/$genome > $sample.fragments.bedgraph

# call peaks
SEACR_1.3.sh ${sample}.fragments.bedgraph 0.01 norm stringent seacr_out

cd ../
```


## Compare peaks

### Jaccard statistic

A first, generic, way to compare bed files is to use Jaccard statistic.
We can do that, by using ```bedtools jaccard``` as follow:
```
bedtools jaccard -a nfcore_approach/seacr_out.stringent.bed -b seacr_approach/seacr_out.stringent.bed 
```

Which returns:
```
intersection	union	jaccard	n_intersections
375126	606045	0.618974	289
```

This means that of the 403 peaks identified in ```nfcore_approach/seacr_out.stringent.bed``` only 289 (~70%) were also found in ```seacr_approach/seacr_out.stringent.bed ```.\
In addition, the Jaccard statistics is 0.62 - which I wouldn't call "great".


### Peak signal

A second, probably better, way to compare the peaks is to compare the SEACR peak signal in the peaks identified by the one or the other method.

To this end, peaks identified by both methods are intersected and the signal from the ones overlapping for at least 1 bp extracted:
```
bedtools intersect -wo -a nfcore_approach/seacr_out.stringent.bed -b seacr_approach/seacr_out.stringent.bed | awk '{ print $1"_"$2"_"$3"@"$7"_"$8"_"$9 "\t" $4 "\t" $10 }' > common_peaks.bed
```

Now we need to add the peaks identified by one tool only. These peaks will be assigned a signal of 0 in the method that failed in identifying them:
```
# identified only by nf_core
bedtools intersect -c -a nfcore_approach/seacr_out.stringent.bed -b seacr_approach/seacr_out.stringent.bed | awk '{ if ($NF==0) print $1"_"$2"_"$3"@NA\t" $4 "\t" 0 }' > nfcore_specific.bed
bedtools intersect -c -a seacr_approach/seacr_out.stringent.bed -b nfcore_approach/seacr_out.stringent.bed | awk '{ if ($NF==0) print "NA@"$1"_"$2"_"$3"\t0\t" $4 }' > seacr_specific.bed
```

Combine them all, and sort:
```
cat common_peaks.bed seacr_specific.bed nfcore_specific.bed | sort -k1,1 > intersected_peaks.bed
```

Clean
```
rm common_peaks.bed *spec*bed
```
