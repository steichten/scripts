Readme file


This is a updating set of scripts that I've pulled together to attempt and make analysis easier.

Methylome work:

fastq > bismark outputs - use wgbs_se_pipeline_v0.2.sh:

USAGE: wgbs_se_pipelinev0.2.sh <input fastq file> <rel. path to bismark genome folder> <fileID for output files>
Example: wgbs_st_pipelinev0.2.sh Bd21.fastq ../brachy_genomes/Bd21Control Bd21

This will create a folder structure in the directory you are in that contains the fastq files, fastqc reports, trimmed + filtered reads, secondary fastqc reports, bismark outputs, and windowed/coverage/bed files for visualization and analysis.



How to setup:

execute in a directory containg:
your fastq file
C_context_window_SREedits.pl

