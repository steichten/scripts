#!/bin/bash
set -u
# Bisulfite sequence analysis pipeline v0.1
# SRE
# Updated 29-5-2014
###################
# This script it designed to take a single end fastq file and process it through the 
# bismark aligner call methylated cytosines, and develop per-c bed files and 100bp 
# window wig files for CpG, CHG, and CHH methylation levels
#
#
###################
#
#usage:
if [ "$#" -ne 3 ]; then
echo "Missing required arguments!"
echo "USAGE: wgbs_se_pipelinev0.2.sh <input fastq file> <path to bismark genome folder> <fileID for output files>"
exit 1
fi


#gather input variables
fq_file=$1; #the input fastq file
genome_path=$2; #the path to the genome to be used (bismark genome prepped)
fileID=$3;
dow=$(date +"%F-%H-%m-%S")

echo "##################"
echo "Performing Bismark single-end alignment with the following parameters:"
echo "Input File: $fq_file"
echo "Path to bismark genome folder: $genome_path"
echo "Output ID: $fileID"
echo "Time of analysis: $dow"
echo ""
echo "Full logfile of steps: ${fileID}_logs_${dow}.log"
echo "##################"

#develop directory tree
mkdir ${fileID}_wgbspipeline_${dow}
mv $fq_file ${fileID}_wgbspipeline_${dow}
cd ${fileID}_wgbspipeline_${dow}

#fastqc
mkdir 1_fastqc
fastqc $fq_file 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv ${fq_file%%.fastq}_fastqc* 1_fastqc

raw_reads=$(wc -l < $fq_file)
raw_reads=$(($raw_reads / 4 ))
echo $raw_reads

#trim_galore
mkdir 2_trimgalore
cd 2_trimgalore
trim_galore ../$fq_file 2>&1 | tee -a ../${fileID}_logs_${dow}.log
cd ../

#fastqc_again
mkdir 3_trimmed_fastqc
fastqc 2_trimgalore/${fq_file%%.fastq}_trimmed.fq 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 2_trimgalore/${fq_file%%.fastq}_trimmed.fq_* 3_trimmed_fastqc
flt_reads=$(wc -l < 2_trimgalore/${fq_file%%.fastq}_trimmed.fq)
flt_reads=$(($flt_reads / 4))
echo $flt_reads

mkdir 0_rawfastq
mv $fq_file 0_rawfastq

#bismark
mkdir 4_bismark_alignment
cd 4_bismark_alignment
bismark -n 2 -l 20 ../../$genome_path ../2_trimgalore/${fq_file%%.fastq}_trimmed.fq 2>&1 | tee -a ../${fileID}_logs_${dow}.log
#sam to bam
samtools view -bSh ${fq_file%%.fastq}_trimmed.fq_bismark.sam > ${fq_file%%.fastq}_trimmed.fq_bismark.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
samtools sort ${fq_file%%.fastq}_trimmed.fq_bismark.bam ${fq_file%%.fastq}_trimmed.fq_bismark.sorted 2>&1 | tee -a ../${fileID}_logs_${dow}.log
samtools index ${fq_file%%.fastq}_trimmed.fq_bismark.sorted.bam 2>&1 | tee -a ../${fileID}_logs_${dow}.log
#sam sort for MethylKit
grep -v '^[[:space:]]*@' ${fq_file%%.fastq}_trimmed.fq_bismark.sam | sort -k3,3 -k4,4n > ${fq_file%%.fastq}_trimmed.fq_bismark.sorted.sam
#methylation extraction
bismark_methylation_extractor --comprehensive --report --buffer_size 8G -s ${fq_file%%.fastq}_trimmed.fq_bismark.sam 2>&1 | tee -a ../${fileID}_logs_${dow}.log

#bedgraph creation
bismark2bedGraph --CX CpG* -o ${fileID}_CpG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark2bedGraph --CX CHG* -o ${fileID}_CHG.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
bismark2bedGraph --CX CHH* -o ${fileID}_CHH.bed 2>&1 | tee -a ../${fileID}_logs_${dow}.log
cd ../
mkdir 5_output_files
mv 4_bismark_alignment/*.bed* 5_output_files

#100bp window creation

perl ../C_context_window_SREedits.pl 4_bismark_alignment/CpG* 100 0 ${fileID}_CpG 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CpG*.wig 5_output_files/${fileID}_CpG_100bp.wig
perl ../C_context_window_SREedits.pl 4_bismark_alignment/CHG* 100 0 ${fileID}_CHG 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CHG*.wig 5_output_files/${fileID}_CHG_100bp.wig
perl ../C_context_window_SREedits.pl 4_bismark_alignment/CHH* 100 0 ${fileID}_CHH 2>&1 | tee -a ${fileID}_logs_${dow}.log
mv 4_bismark_alignment/CHH*.wig 5_output_files/${fileID}_CHH_100bp.wig

echo "#####################"
echo "providing pipeline metrics to wgbs pipeline logfile..."
echo "#####################"

#get all relevant numbers for final log summary
bismark_version=$(bismark --version | grep "Bismark Version:" | cut -d":" -f2 | tr -d ' ')
samtools_version=$(samtools 3>&1 1>&2 2>&3 | grep "Version:" | cut -d' ' -f2 | tr -d ' ')

map_ef=$(grep 'Mapping efficiency:' 4_bismark_alignment/${fq_file%%.fastq}_trimmed.fq_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
unique_aln=$(grep 'Number of alignments with a unique best hit from the different alignments:' 4_bismark_alignment/${fq_file%%.fastq}_trimmed.fq_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
no_aln=$(grep 'Sequences with no alignments under any condition:' 4_bismark_alignment/${fq_file%%.fastq}_trimmed.fq_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
multi_aln=$(grep 'Sequences did not map uniquely:' 4_bismark_alignment/${fq_file%%.fastq}_trimmed.fq_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t')
cpg_per=$(grep 'C methylated in CpG context:' 4_bismark_alignment/${fq_file%%.fastq}_trimmed.fq_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chg_per=$(grep 'C methylated in CHG context:' 4_bismark_alignment/${fq_file%%.fastq}_trimmed.fq_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)
chh_per=$(grep 'C methylated in CHH context:' 4_bismark_alignment/${fq_file%%.fastq}_trimmed.fq_bismark_SE_report.txt  | cut -d: -f2 | tr -d '\t' | cut -d'%' -f1)

#add it to the full pipeline logfile
printf "${dow}\t${fq_file}\t${fileID}\t${genome_path##../}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n"
printf "${dow}\t${fq_file}\t${fileID}\t${genome_path##../}\t${bismark_version}\t${samtools_version}\t${raw_reads}\t${flt_reads}\t${map_ef}\t${unique_aln}\t${no_aln}\t${multi_aln}\t${cpg_per}\t${chg_per}\t${chh_per}\n" >> /home/steve/wgbs_se_pipeline_analysis_record.log






