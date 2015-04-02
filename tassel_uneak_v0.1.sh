#!/bin/bash
set -u

#############################
# TASSEL UNEAK Pipeline Script
# v0.1
# Borevitz Lab - SRE
###############

usage() { 
echo "############################################################"
echo
echo "USAGE: tassel_uneak_vX.X.sh "
echo "		-f <path to fastq directory> "
echo "		-k <path to keyfile> "
echo "		-e <enzyme used for GBS digestion>"
echo
echo "This script will perform the TASSEL SNP Discovery pipeline"
echo
echo "REQUIREMENTS:"
echo "TASSEL v3"
echo "BWA or Bowtie2"
echo "Indexed Reference genome to align to (fasta format)"
echo "Keyfile formatted for TASSEL sample demuxing (ends in '_key.txt')"
echo
echo "############################################################"
exit 1
}

flag1=0
flag2=0
flag3=0
flag4=0
flag5=0
flag6=0
while getopts ":f:k:e:o:" opt; do
	case $opt in
     f)  fastqdir=$OPTARG; flag2=1;;
     k)  keyfile=$OPTARG; flag3=1;;
     e)  enzyme=$OPTARG; flag6=1;;
     o)  prefix=$OPTARG; flag1=1;;
    \?)  usage;;
     :)  echo "Option -$OPTARG requires an argument." >&2; usage;;	
     *)  usage
	esac
done

if [ $flag1 == 0 ]; then
	echo "############################################################"
	echo "output prefix ( -o ) required!"
fi
if [ $flag2 == 0 ]; then
	echo "############################################################"
	echo "path to directory containing fastq files ( -f ) required!"
fi
if [ $flag3 == 0 ]; then
	echo "############################################################"
	echo "tassel keyfile ( -k ) required!"
fi
if [ $flag6 == 0 ]; then
	echo "############################################################"
	echo "GBS enzyme digestion ( -e ) required!"
fi

if [ $((flag1+flag2+flag3+flag6)) != 4 ]; then
	usage
fi

timestamp=$(date +'%d-%m-%Y-%H%M')

echo "###################################"
echo "# TASSEL UNEAK Pipeline Script"
echo "# v0.1"
echo "################# Arguments provided:"
echo "time=${timestamp}"
echo "f = ${fastqdir}"
echo "k = ${keyfile}"
echo "e = ${enzyme}"
echo "o = ${prefix}"

mkdir ${prefix}_UNEAK_${timestamp}
cd ${prefix}_UNEAK_${timestamp}

#make directory structure
run_pipeline.pl -fork1 -UCreatWorkingDirPlugin -w . -endPlugin -runfork1

#move in fastq to /illumina
cp ../${fastqdir}/* Illumina

#move in key to /key
cp ../${keyfile} key

#UFastqToTagCountPlugin
run_pipeline.pl -Xmx32g -fork1 -UFastqToTagCountPlugin -w . -e ApeKI -endPlugin -runfork1

#UMerge Taxa Tag Count
run_pipeline.pl -Xmx32g -fork1 -UMergeTaxaTagCountPlugin -w . -c 5 -endPlugin -runfork1


#UTag Count To Tag Pair
run_pipeline.pl -Xmx32g -fork1 -UTagCountToTagPairPlugin -w . -e 0.03 -endPlugin -runfork1

#UTag Pair To TBT
run_pipeline.pl -Xmx32g -fork1 -UTagPairToTBTPlugin -w . -endPlugin -runfork1

#U TBT To Map Info
run_pipeline.pl -Xmx32g -fork1 -UTBTToMapInfoPlugin -w . -endPlugin -runfork1

#U Map Info To HapMap
run_pipeline.pl -Xmx32g -fork1 -UMapInfoToHapMapPlugin -w . -mnMAF 0.05 -mxMAF 0.5 -mnC 0 -mxC 1 -endPlugin -runfork1

