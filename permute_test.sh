#!/bin/bash
set -u

usage() { 
echo "############################################################"
echo
echo "Help & usage TBD" 1>&2
echo "############################################################"
exit 1
}

flag1=0
flag2=0
flag3=0
while getopts ":n:i:c:" opt; do
	case $opt in
     n)  perm=$OPTARG; flag1=1;;
     i)  file=$OPTARG; flag2=1;;
     c)  compare=$OPTARG; flag3=1;;
    \?)  usage;;
     :)  echo "Option -$OPTARG requires an argument." >&2; usage;;	
     *)  usage
	esac
done

if [ $flag1 == 0 ]; then
	echo "############################################################"
	echo "permutation number argument ( -n ) required!"
fi
if [ $flag2 == 0 ]; then
	echo "############################################################"
	echo "input file argument ( -i ) required!"
fi
if [ $flag3 == 0 ]; then
	echo "############################################################"
	echo "comparison file argument ( -c ) required!"
fi
if [ $((flag1+flag2+flag3)) != 3 ]; then
	usage
fi


#####################
output=NULL

bedtools closest -D "ref" -a ${file} -b ${compare} > nonrandom.temp.to.genes.txt
output=$(awk -F$'\t' '$NF<500 && $NF>-500' nonrandom.temp.to.genes.txt | wc | awk '{print $1}')
echo "0	${output}"

for i in `seq 1 ${perm}`;
do
	bedtools shuffle -i ${file} -seed ${i} -g /home/steve/genomes/Bd21Control/chromosome_sizes.txt > random.temp.bed
	bedtools closest -D "ref" -a random.temp.bed -b ${compare} > random.temp.to.genes.txt
	output=$(awk -F$'\t' '$NF<500 && $NF>-500' random.temp.to.genes.txt | wc | awk '{print $1}')
	echo "${i}	${output}"
done

rm nonrandom.temp.to.genes.txt
rm random.temp.bed
rm random.temp.to.genes.txt