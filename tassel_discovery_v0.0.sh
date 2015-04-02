#!/bin/bash
set -u

#############################
# TASSEL Discovery Pipeline Script
# v0.1
# Borevitz Lab - SRE
###############

usage() { 
echo "############################################################"
echo
echo "USAGE: tassel_discovery_vX.X.sh "
echo "		-a <aligner to use bwa or bt2> "
echo "		-f <path to fastq directory> "
echo "		-k <path to keyfile> "
echo "		-e <enzyme used for GBS digestion>"
echo "		-r <path to reference> "
echo "		-o <output name prefix>"
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
while getopts ":a:f:k:e:r:o:" opt; do
	case $opt in
     a)  aligner=$OPTARG; flag1=1;;
     f)  fastqdir=$OPTARG; flag2=1;;
     k)  keyfile=$OPTARG; flag3=1;;
     e)  enzyme=$OPTARG; flag6=1;;
     r)  reference=$OPTARG; flag4=1;;
     o)  prefix=$OPTARG; flag5=1;;
    \?)  usage;;
     :)  echo "Option -$OPTARG requires an argument." >&2; usage;;	
     *)  usage
	esac
done

if [ $flag1 == 0 ]; then
	echo "############################################################"
	echo "read aligner method ( -a ) required!"
fi

if [ $flag2 == 0 ]; then
	echo "############################################################"
	echo "path to directory containing fastq files ( -f ) required!"
fi
if [ $flag3 == 0 ]; then
	echo "############################################################"
	echo "tassel keyfile ( -k ) required!"
fi
if [ $flag4 == 0 ]; then
	echo "############################################################"
	echo "reference genome ( -r ) required!"
fi
if [ $flag5 == 0 ]; then
	echo "############################################################"
	echo "output naming prefix ( -o ) required!"
fi
if [ $flag6 == 0 ]; then
	echo "############################################################"
	echo "GBS enzyme digestion ( -e ) required!"
fi

if [ $((flag1+flag2+flag3+flag4+flag5+flag6)) != 6 ]; then
	usage
fi

timestamp=$(date +'%d-%m-%Y-%H%M')

echo "###################################"
echo "# TASSEL Discovery Pipeline Script"
echo "# v0.1"
echo "################# Arguments provided:"
echo "time=${timestamp}"
echo "a = ${aligner}"
echo "f = ${fastqdir}"
echo "k = ${keyfile}"
echo "e = ${enzyme}"
echo "r = ${reference}"
echo "o = ${prefix}"

#make directory structure
mkdir ${prefix}_discovery_${timestamp}
cd ${prefix}_discovery_${timestamp}
mkdir 1_fastqtotagcount
mkdir 2_mergedtagcount
mkdir 3_tagcounttofastq
mkdir 4_mergedtagcounts
mkdir 5_topm
mkdir 6_tbt
mkdir 7_hapmap
cd 7_hapmap
mkdir raw
mkdir mergedSNPs
mkdir filt
cd ../


#FastQ To Tag Count
run_pipeline.pl -Xmx32g -fork1 -FastqToTagCountPlugin -i ../${fastqdir} -k ../${keyfile} -e ${enzyme} -o 1_fastqtotagcount -endPlugin -runfork1

#Merge Multiple Tag Count
run_pipeline.pl -Xmx32g -fork1 -MergeMultipleTagCountPlugin -i 1_fastqtotagcount -o 2_mergedtagcount/${prefix}_masterGBStags.cnt -c 5 -endPlugin -runfork1

#Tag Count To Fastq
run_pipeline.pl -Xmx32g -fork1 -TagCountToFastqPlugin -i 2_mergedtagcount/${prefix}_masterGBStags.cnt -o 3_tagcounttofastq/${prefix}_masterGBStags.fq -c 5 -endPlugin -runfork1

if [ "$aligner" == "bwa" ]; then

#BWA Alignment
bwa aln -t 4 ../${reference} 3_tagcounttofastq/${prefix}_masterGBStags.fq > 4_mergedtagcounts/${prefix}_AlignedMasterTags.sai

#BWA sai to sam file conversion
bwa samse ../${reference} 4_mergedtagcounts/${prefix}_AlignedMasterTags.sai 3_tagcounttofastq/${prefix}_masterGBStags.fq > 4_mergedtagcounts/${prefix}_AlignedMasterTags.sam

fi

if [ "$aligner" == "bt2" ]; then

#bowtie2 alignment
bowtie2 -M 4 -p 15 --very-sensitive-local -x ../${reference%%.fasta} -U 3_tagcounttofastq/${prefix}_masterGBStags.fq -S 4_mergedtagcounts/${prefix}_AlignedMasterTags.sam

fi

#Sam Convertor
run_pipeline.pl -Xmx32g -fork1 -SAMConverterPlugin -i 4_mergedtagcounts/${prefix}_AlignedMasterTags.sam -o 5_topm/${prefix}_MasterTags.topm -endPlugin -runfork1

#Fastq To TBT
run_pipeline.pl -Xmx32g -fork1 -FastqToTBTPlugin -i ../${fastqdir} -k ../${keyfile} -e ${enzyme} -o 6_tbt -y -t 2_mergedtagcount/${prefix}_masterGBStags.cnt -endPlugin -runfork1

#Merge Tags By Taxa
run_pipeline.pl -Xmx32g -fork1 -MergeTagsByTaxaFilesPlugin -i 6_tbt -o 6_tbt/${prefix}_merged.tbt.byte -endPlugin -runfork1

#Tags To SNP By Alignment
run_pipeline.pl -Xmx32g -fork1 -TagsToSNPByAlignmentPlugin -i 6_tbt/${prefix}_merged.tbt.byte -y -m 5_topm/${prefix}_MasterTags.topm -mUpd 5_topm/${prefix}_MasterTagsWithVariants.topm -o 7_hapmap/raw/${prefix}_GBSGenos_chr+.hmp.txt.gz -mnF 0.8 -mnMAF 0.02 -mnMAC 10000 -ref ../${reference} -sC 1 -eC 10 -endPlugin -runfork1

#Merge Duplicate SNPs
run_pipeline.pl -Xmx32g -fork1 -MergeDuplicateSNPsPlugin -hmp 7_hapmap/raw/${prefix}_GBSGenos_chr+.hmp.txt.gz -o 7_hapmap/mergedSNPs/${prefix}_GBSGenos_mergedSNPs_chr+.hmp.txt.gz -misMat 0.1 -callHets -sC 1 -eC 10 -endPlugin -runfork1

#GBS HapMap Filters
#TBD - largely dependent on experimental goal



