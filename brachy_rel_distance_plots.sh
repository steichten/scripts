#!/usr/bin/bash

#plotting DNA methylation over Brachypodium gene models
#12-5-14
#SRE
#######################

filename=$1
#######################
closestBed -D "ref" -a ${filename}_CHG.bed -b ../brachy_annotations/Bdistachyon.MIPS_1_2.predicted.genes.bed > ${filename}_CHG_gene.bed
closestBed -D "ref" -a ${filename}_CHH.bed -b ../brachy_annotations/Bdistachyon.MIPS_1_2.predicted.genes.bed > ${filename}_CHH_gene.bed
closestBed -D "ref" -a ${filename}_CpG.bed -b ../brachy_annotations/Bdistachyon.MIPS_1_2.predicted.genes.bed > ${filename}_CpG_gene.bed

#subset to the regions within 100bp of a gene (make the files more manageable for R)
awk -F$'\t' '$10<1000 && $10>-1000' ${filename}_CHG_gene.bed > ${filename}_CHG_gene.1k.bed
awk -F$'\t' '$10<1000 && $10>-1000' ${filename}_CHH_gene.bed > ${filename}_CHH_gene.1k.bed
awk -F$'\t' '$10<1000 && $10>-1000' ${filename}_CpG_gene.bed > ${filename}_CpG_gene.1k.bed
#######################

#initiate the R script to create the plots
Rscript gene_methylation_plots.r ${filename}