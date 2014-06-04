#!/bin/bash
set -u

export PATH=$PATH:/home/steve/bin

# MethylKit-edmr DMR calling from bismark sam files v0.2
# SRE
# Updated 3-6-2014
###################
#
#
#
###################
#
#usage:
if [ "$#" -ne 3 ]; then
echo "USAGE: wgbs_dmr_pipeline.sh <input path to sam file dir> <CpG/CHG/CHH> <minimum coverage>"
echo "EXAMPLE: wgbs_dmr_pipeline.sh sorted_sam CpG 1"
exit 1
fi


#gather input variables
sam_dir=$1; #the input fastq file
dmr_context=$2 #context to pull cytosines for DMR analysis (CpG / CHG / CHH)
min_cov=$3 #minimum coverage required to analyze the cytosine
dow=$(date +"%F-%H-%m-%S") #timestamp

#develop directory structure
cd $sam_dir
#mkdir ${dow}_${dmr_context}_mincov_${min_cov}
#mv *.sam ${dow}_${dmr_context}_mincov_${min_cov}
cd ${dow}_${dmr_context}_mincov_${min_cov}
mkdir 1_MethylKit_results
mkdir 2_eDMR
mkdir 3_eDMR_filtered
mkdir 4_DMRmet_results
#call R script for methylkit and eDMR work
Rscript /home/steve/scripts/methylkit_edmr_pipeline_v0.2.r $dmr_context $min_cov

	
#bedops to merge all of the filtered eDMR files into one. Merge to make the largest boundary for any overlapping DMRs
cd 3_eDMR_filtered
sort-bed *.bed | bedops -m - > merged.edmr.regions.bed
cp merged.edmr.regions.bed ../4_DMRmet_results

cd ../4_DMRmet_results
cp ../../covfiles/*${dmr_context}*.cov .

#bedtools to grab methylation coverage for each samples (match to your sam file set)
for file in *.cov
do
	bedtools intersect -wa -wb -a merged.edmr.regions.bed -b "$file" | bedtools groupby -i stdin -g 1,2,3 -c 7,8,9 -o mean,sum,sum > "${file}.dmr"
done

############
#R to merge all the .dmr files together into the final output

Rscript /home/steve/scripts/dmr_merge.r