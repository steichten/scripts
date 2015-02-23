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

if [ "$#" -ne 4 ]; then
echo "USAGE: bed_to_rel_dist.sh <-wig or -bed> <input path to bed file you want to map to> <filename prefix> <output map name>"
echo "EXAMPLE: bed_to_rel_dist.sh -bed /home/steve/brachy_annotation/genes.bed epignome.500k gene"
exit 1
fi
type=$1
bedpath=$2
filename=$3
outname=$4

if [ "$1" == "-wig" ];then

sed -e "1d" ${filename}_CpG_100bp.wig > ${filename}_CpG_100bp.bed
sed -e "1d" ${filename}_CHG_100bp.wig > ${filename}_CHG_100bp.bed
sed -e "1d" ${filename}_CHH_100bp.wig > ${filename}_CHH_100bp.bed

fi

sort -k1,1 -k2,2n ${filename}_CpG_100bp.bed -o ${filename}_CpG_100bp.bed
sort -k1,1 -k2,2n ${filename}_CHG_100bp.bed -o ${filename}_CHG_100bp.bed
sort -k1,1 -k2,2n ${filename}_CHH_100bp.bed -o ${filename}_CHH_100bp.bed

#get total number of columns for both input files
l1="$(cat ${filename}_CpG_100bp.bed | awk 'BEGIN{FS="\t"};{print NF}' | head -n 1)"
l2="$(cat ${bedpath} | awk 'BEGIN{FS="\t"};{print NF}' | head -n 1)"


#convert the wigs to bed
####################### some bedtools stuff
echo "Performing closestBed of CHG methylation..."
closestBed -D "ref" -a ${filename}_CHG_100bp.bed -b $bedpath > ${filename}_CHG_${outname}.bed
echo "Performing closestBed of CHH methylation..."
closestBed -D "ref" -a ${filename}_CHH_100bp.bed -b $bedpath > ${filename}_CHH_${outname}.bed
echo "Performing closestBed of CpG methylation..."
closestBed -D "ref" -a ${filename}_CpG_100bp.bed -b $bedpath > ${filename}_CpG_${outname}.bed

#subset to the regions within 100bp of a gene (make the files more manageable for R)
echo "subsetting files to within 1kb..."
awk -F$'\t' '$NF<1000 && $NF>-1000' ${filename}_CHG_${outname}.bed > ${filename}_CHG_${outname}.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' ${filename}_CHH_${outname}.bed > ${filename}_CHH_${outname}.1k.bed
awk -F$'\t' '$NF<1000 && $NF>-1000' ${filename}_CpG_${outname}.bed > ${filename}_CpG_${outname}.1k.bed
#######################


echo "Performing R plots..."
#initiate the R script to create the plots
Rscript $HOME/scripts/rel_methylation_plots.r ${filename} ${outname} ${l1} ${l2}
