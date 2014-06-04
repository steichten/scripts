#!/bin/bash

# Bisulfite sequence analysis pipeline v0.2
# SRE
# Updated 6-5-2014
###################
# This script it designed to take a single end fastq file and process it through the 
# bismark aligner call methylated cytosines, and develop per-c bed files and 100bp 
# window wig files for CpG, CHG, and CHH methylation levels
#
#
###################
#
#usage:
#bisulfite_se_mapping.sh <input fastq file> <path to genome folder> <file ID for output files>



#gather input variables
fq_file=$1;
genome_path=$2;
fileID=$3;
dow=$(date +"%F")
echo "$dow"

echo $fq_file
echo $genome_path
echo $fileID
echo newfile_named_${fileID}

#develop directory tree
mkdir ${fileID}_wgbspipeline_${dow}
mv $fq_file ${fileID}_wgbspipeline_${dow}
cd ${fileID}_wgbspipeline_${dow}

#fastqc
mkdir 1_fastqc
fastqc $fq_file
mv ${fq_file%%.fastq}_fastqc* 1_fastqc

#trim_galore
mkdir 2_trimgalore
cd 2_trimgalore
trim_galore ../$fq_file
cd ../

#fastqc_again
mkdir 3_trimmed_fastqc
fastqc 2_trimgalore/${fq_file%%.fastq}_trimmed.fq
mv 2_trimgalore/${fq_file%%.fastq}_trimmed.fq_* 3_trimmed_fastqc
mkdir 0_rawfastq
mv $fq_file 0_rawfastq

#bismark
mkdir 4_bismark_alignment
cd 4_bismark_alignment
bismark -n 2 -l 20 ../../$genome_path ../2_trimgalore/${fq_file%%.fastq}_trimmed.fq
#sam to bam
samtools view -bS ${fq_file%%.fastq}_trimmed.fq_bismark.sam > ${fq_file%%.fastq}_trimmed.fq_bismark.bam
samtools sort ${fq_file%%.fastq}_trimmed.fq_bismark.bam ${fq_file%%.fastq}_trimmed.fq_bismark.sorted
samtools index ${fq_file%%.fastq}_trimmed.fq_bismark.sorted.bam

#methylation extraction
bismark_methylation_extractor --comprehensive --report --buffer_size 8G -s ${fq_file%%.fastq}_trimmed.fq_bismark.sam

#bedgraph creation
bismark2bedGraph --CX CpG* -o ${fileID}_CpG.bed
bismark2bedGraph --CX CHG* -o ${fileID}_CHG.bed
bismark2bedGraph --CX CHH* -o ${fileID}_CHH.bed
cd ../
mkdir 5_output_files
mv 4_bismark_alignment/*.bed* 5_output_files

#100bp window creation

perl ../C_context_window.pl 4_bismark_alignment/CpG* 100 0 ${fileID}_CpG | tee -a ../${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CpG*.wig 5_output_files/${fileID}_CpG_100bp.wig
perl ../C_context_window.pl 4_bismark_alignment/CHG* 100 0 ${fileID}_CHG | tee -a ../${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CHG*.wig 5_output_files/${fileID}_CHG_100bp.wig
perl ../C_context_window.pl 4_bismark_alignment/CHH* 100 0 ${fileID}_CHH | tee -a ../${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CHH*.wig 5_output_files/${fileID}_CHH_100bp.wig
