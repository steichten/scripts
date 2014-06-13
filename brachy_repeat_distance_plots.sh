#!/usr/bin/bash

#plotting DNA methylation over Brachypodium repeat models
#12-5-14
#SRE
#######################

filename=$1
####################### bedtools to map per-cytosine levels to repeats
closestBed -D "ref" -a ${filename}_CHG.bed -b ../brachy_annotations/Bd21Control_repeatmasker_v1_200bp.bed > ${filename}_CHG_repeat.bed
closestBed -D "ref" -a ${filename}_CHH.bed -b ../brachy_annotations/Bd21Control_repeatmasker_v1_200bp.bed > ${filename}_CHH_repeat.bed
closestBed -D "ref" -a ${filename}_CpG.bed -b ../brachy_annotations/Bd21Control_repeatmasker_v1_200bp.bed > ${filename}_CpG_repeat.bed

#subset to the regions within 100bp of a repeat (make the files more manageable for R)
awk -F$'\t' '$10<1000 && $10>-1000' ${filename}_CHG_repeat.bed > ${filename}_CHG_repeat.1k.bed
awk -F$'\t' '$10<1000 && $10>-1000' ${filename}_CHH_repeat.bed > ${filename}_CHH_repeat.1k.bed
awk -F$'\t' '$10<1000 && $10>-1000' ${filename}_CpG_repeat.bed > ${filename}_CpG_repeat.1k.bed
#######################

#initiate the R script to create the plots
Rscript repeat_methylation_plots.r ${filename}
