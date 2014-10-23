#!/bin/bash
set -u

export PATH=$PATH:/home/steve/bin
export PATH=$PATH:/home/steve/bin/
# 100bp_dmrs.sh
# This script identifies 100bp windows (from bismark alignment pipeline) that display 
# a difference (of choice) in methylation (context of choice) with a required coverage
# level (of choice). These windows are selected across every pairwise comparison of 100bp 
# window wig files and aggregated into a single (bed) file. All adjacent windows are 
# collapsed into a single DMR. This file is then used to grab all individual met/unmet 
# reads for each DMR (from .cov files in bismark pipeline) for all samples.


#REQUIRES THAT sample names do not contain '_', as this will screw up the final steps
######################

#execute from directory containing all the wig and cov files from all samples
#usage:
if [ "$#" -ne 3 ]; then
echo "USAGE: 100bp_dmrs.v0.1.sh <context> <diffmeth> <coverage_req>"
echo "EXAMPLE: 100bp_dmrs.v0.1.sh CpG 80 10"
echo "Look at CG context, difference of 80% methylation with 10 reads minimum over window"
exit 1
fi
######################
#define arguments
context=$1
difference=$2
coverage=$3
######################
Rscript /home/steve/scripts/100bp_wig_to_dmrs.r ${context} ${difference} ${coverage}


#bedtools to intersect the bed file w. the coverage files
for file in *${context}*.cov
do
	bedtools intersect -wa -wb -a 100bp_${context}_${difference}diff_${coverage}collapsed.bed -b "$file" | bedtools groupby -i stdin -g 4 -c 1,2,3,5,9,10,11 -o first,mean,mean,mean,mean,sum,sum > "${file}.${context}_${difference}diff_${coverage}.dmr"
done


#file structure cleanup
mkdir 100bp_${context}_${difference}diff_${coverage}_out
mv 100bp_${context}_${difference}diff_${coverage}* 100bp_${context}_${difference}diff_${coverage}_out/
mv *.${context}_${difference}diff_${coverage}* 100bp_${context}_${difference}diff_${coverage}_out/
cd 100bp_${context}_${difference}diff_${coverage}_out/

#collapse them all into a summary table
Rscript /home/steve/scripts/100bp_dmr_merge.r ${context} ${difference} ${coverage}
