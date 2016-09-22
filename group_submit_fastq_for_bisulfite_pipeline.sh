
#grab each fastq.gz file and run them through the alignment pipeline. Output filenames are the SampleID number for BVZ0049

N=12
for FILE in *.fastq.gz
	do
	((i=i%N)); ((i++==0)) && wait
	../scripts/wgbs_pipelinev0.4.sh -se_epi $FILE ../genomes_temp/Bd21Control $(echo $FILE | cut -d "_" -f1) &
done

mkdir 100bp_tiles
find -name "*.wig" | xargs mv -t 100bp_tiles

#copy tiles to local machine for IGV check