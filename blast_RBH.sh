

#make blast databases
makeblastdb -in Bdistachyon_192_cds_primaryOnly.fa -dbtype nucl -parse_seqids
makeblastdb -in Bstacei_316_v1.1.cds_primaryTranscriptOnly.fa -dbtype nucl -parse_seqids


#run blast with output format 6
blastn -db Bstacei_316_v1.1.cds_primaryTranscriptOnly.fa -query Bdistachyon_192_cds_primaryOnly.fa -outfmt 6 -out dis_to_stadb.out2

blastn -db Bdistachyon_192_cds_primaryOnly.fa -query Bstacei_316_v1.1.cds_primaryTranscriptOnly.fa -outfmt 6 -out sta_to_disdb.out2


chmod 770 reciprocal_blast_hits.py

#run script to find RBHs

./reciprocal_blast_hits.py dis_to_stadb.out2 sta_to_disdb.out2 1 2 12 high dis_sta.hits.out