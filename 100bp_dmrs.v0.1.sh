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
#if [ "$#" -ne 3 ]; then
#echo "USAGE: 100bp_dmrs.v0.1.sh <context> <diffmeth> <coverage_req>"
#echo "EXAMPLE: 100bp_dmrs.v0.1.sh CpG 80 10"
#echo "Look at CG context, difference of 80% methylation with 10 reads minimum over window"
#exit 1
#fi

usage() { 
echo "############################################################"
echo
echo "Usage: $0 [-c <CpG | CHG | CHH>] [-m <0|100>] [-d <integer>]" 1>&2
echo
echo "EXAMPLE: $0 -c CpG -m 80 -d 10"
echo "Look at CpG context, difference of 80% methylation with 10 reads minimum over window"
echo
echo "This script will create DMRs from 100bp window wig files"
echo "Execute the script from a directory containing wig files and cov files"
echo
echo "This script utilizes R and the package 'fields'"
echo
echo "############################################################"
exit 1
}

flag1=0
flag2=0
flag3=0
while getopts ":c:m:d:" opt; do
	case $opt in
     c)  context=$OPTARG; flag1=1;;
     m)  difference=$OPTARG; flag2=1;;
     d)  coverage=$OPTARG; flag3=1;;
    \?)  usage;;
     :)  echo "Option -$OPTARG requires an argument." >&2; usage;;	
     *)  usage
	esac
done

if [ $flag1 == 0 ]; then
	echo "############################################################"
	echo "context argument ( -c ) required!"
fi
if [ $flag2 == 0 ]; then
	echo "############################################################"
	echo "methylation difference argument ( -m ) required!"
fi
if [ $flag3 == 0 ]; then
	echo "############################################################"
	echo "coverage argument ( -d ) required!"
fi

if [ $((flag1+flag2+flag3)) != 3 ]; then
	usage
fi

echo "c = ${context}"
echo "m = ${difference}"
echo "d = ${coverage}"

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
