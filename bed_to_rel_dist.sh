#!/bin/bash

#plotting DNA methylation over Brachypodium gene models
#12-5-14
#SRE
#######################
# REQUIREMENTS
# bedtools
# awk
# R with libraries: fields
#######################

if [ "$#" -ne 3 ]; then
echo "USAGE: bed_to_rel_dist.sh <input path to bed file you want to map to> <filename prefix> <output map name>"
echo "EXAMPLE: bed_to_rel_dist.sh /home/steve/brachy_annotation/genes.bed epignome.500k gene"
exit 1
fi

bedpath=$1
filename=$2
outname=$3

####################### some bedtools stuff
echo "Performing closestBed of CHG methylation..."
closestBed -D "ref" -a ${filename}_CHG.bed -b $bedpath > ${filename}_CHG_${outname}.bed
echo "Performing closestBed of CHH methylation..."
closestBed -D "ref" -a ${filename}_CHH.bed -b $bedpath > ${filename}_CHH_${outname}.bed
echo "Performing closestBed of CpG methylation..."
closestBed -D "ref" -a ${filename}_CpG.bed -b $bedpath > ${filename}_CpG_${outname}.bed

#subset to the regions within 100bp of a gene (make the files more manageable for R)
echo "subsetting files to within 1kb..."
awk -F$'\t' '$NF<1000 && $NF>-1000' ${filename}_CHG_${outname}.bed > ${filename}_CHG_${outname}.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' ${filename}_CHH_${outname}.bed > ${filename}_CHH_${outname}.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' ${filename}_CpG_${outname}.bed > ${filename}_CpG_${outname}.1k.bed
#######################
echo "Performing R plots..."
#initiate the R script to create the plots
Rscript rel_methylation_plots.r ${filename} ${outname}
