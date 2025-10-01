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

### Plot results - this needs to be executed in an R/Rstudio env
```
################################################################################
library(ggplot2)

# Read file containing peaks from both methods and their signal
res <- read.table("intersected_peaks.bed",sep='\t',header=F)
names(res) <- c("peak","nfcore","seacr")

################################################################################
# Get metrics

# number of peaks identified only by nfcore approach
n1 <- nrow(subset(res,nfcore>0 & seacr==0))
n2 <- nrow(subset(res,nfcore==0 & seacr>0))

# total number of peaks
t1 <- nrow(subset(res,nfcore>0))
t2 <- nrow(subset(res,seacr>0))

# shared
s1 <- nrow(subset(res,nfcore>0 & seacr>0))
################################################################################

ggplot(res,aes(x=nfcore,y=seacr)) +
  theme_classic() +
  geom_point(size=2,shape=21,fill="indianred",color='black') +
  # x axis
  xlab("nfcore approach") +
  theme(axis.title.x = element_text(size = 12,color = "black")) +
  theme(axis.text.x = element_text(size = 10)) +
  # y axis
  ylab("seacr GitHub approach") +
  theme(axis.title.y = element_text(size = 12,color = "black")) +
  theme(axis.text.y = element_text(size = 10)) +
  coord_cartesian(xlim = c(0,max(max(res$nfcore),max(res$seacr))),
                  ylim = c(0,max(max(res$nfcore),max(res$seacr)))) +
  # add line
  geom_abline(intercept = 1, slope = 1) +
  # add cor test results
  annotate("text",x=20000,y=41000,
           label=paste0("R = ",
                        round(cor.test(res$nfcore,res$seacr)$estimate,2),
                        "; p = ",
                        formatC(cor.test(res$nfcore,res$seacr)$p.value,digits = 2,format = 'e')),
           color='red') +
  # add metrics
  annotate("text",x=4000,y=30000,label=paste0("SEACR\nTot peaks = ",t2,
                                              "\nSpecific = ",n2,
                                              "\nShared = ",s1),size=3) +
  annotate("text",x=35000,y=3000,label=paste0("nfcore\nTot peaks = ",t1,
                                              "\nSpecific = ",n1,
                                              "\nShared = ",s1),size=3)
```

<img src="https://github.com/fansalon/testing_peak_calling_cutntag/blob/main/results/cor_res.png" width="750" height="750"/>

The figure above shows a discrete degree of correlation. Although the correlation is statistically significant (p=1.15e-89), the correlation score is good, but not super high (R=0.78) - considering we are just testing two different methods on the very same sample!

In particular, what worries me the most is that there are 114 peaks that are identified only by the nfcore approach, and not by the other. Whereas, on the other hand, there are only 36 peaks identified only by the SEACR GitHub approach. This means that almost all the peaks identified by SEACR GitHub are also identified by the nfcore approach, while the vice-versa is not true. This suggests that the nfcore approach idenitifies more peaks: either nfcore approach is more sensitive and these peaks are real peaks, or it is less specific and these peaks are false positives.

On the positive side, all the tool-specific peaks are associated to lower CUT&Tag signal (<1000).


### IGV screenshots
Defining which of the two approaches is correct is hard, if not impossible. But a perhaps good way to test this is to visually inspect some IGV screenshots trying to define why and where some peaks are called only when following the nfcore approach.

To do so, I have uploaded H3K27ac data for the sample in analysis (chr19 only): bedgraph files from both SEACR GitHub and nfcore methods, bed files containing the identified peaks as well as the bam file.

Then, I started looking for some peaks identified by the nfcore approach only, and found one at Ndufv1 promoter:

<img src="https://github.com/fansalon/testing_peak_calling_cutntag/blob/main/results/screenshot_Ndufv1.png" width="1000" height="750"/>

The bedgraph tracks from the two approaches look quite similar, however, the nfcore-generated bigwig (light blue) display higher signal than the SEACR GitHub one (light red). The reason for this is a consequence of how the tracks were generated. nfcore approach worked at the level of each single read (*i.e.*, each pair was counted twice) while the approach suggested by SEACR GitHub worked at the level of pairs (*i.e.*, each pair was counted once).

We can have a direct proof of this when looking at the peak summit (chr19:4,062,652-4,062,672).

```

# define peak summit
echo -e "19\t4062652\t4062672" > peak_summit.bed

# extract signal from bedgraphs relative to the peak summit
bedtools intersect -wao -a peak_summit.bed -b nfcore_approach/chr1*bedgraph
# --> 19	4062652	4062672	19	4062652	4062672	11	20

bedtools intersect -wao -a peak_summit.bed -b seacr_approach/chr1*bedgraph
# --> 19	4062652	4062672	19	4062637	4062672	8	20

```

From the code above, it is clear how the nfcore approach counted 11 reads in the peak summit region, whereas the SEACR GitHub approach counted only 8. When checking this on IGV by looking at the reads in the BAM file mapping to the peak summit, it appears evident how there are 8 reads in total, 3 of which are mapped as pairs (pink AND violet together/overlapping) and 5 mapped as single reads (pink OR violet).

<img src="https://github.com/fansalon/testing_peak_calling_cutntag/blob/main/results/screenshot_Ndufv1_summit.png" width="1000" height="750"/>

These were counted 8 by the SEACR GitHub approach (each pair is counted once [=3], each single read is counted once [=5]) and 11 by the nfcore approach (each pair is counted twice [=6], each single read is counted once [=5]).


The way the two different approaches are implemented justifies the different number of reads associated to the peak summit, and to the whole peak in general. However, this also affects the way SEACR calls peaks: probably SEACR considered the signal from the nfocre approach enough to call a peak in that region and the one from the SECAR GitHub approach one not enough to call a peak.



