#!/bin/bash

#run this in the directory containing all of the pipeline output folders for all samples (i.e. the ./reads folder)

#decompress them
find -name "*bismark.txt.gz" > list.to.uncompress.txt
parallel -a list.to.uncompress.txt gzip -d {}

find -name "*bismark.txt" > list.to.rewindow.txt
sed 's/.*C/C/' list.to.rewindow.txt | sed 's/_trimmed.*//' | sed 's/context_//' > wig.name.list.txt

parallel --xapply -a list.to.rewindow.txt -a wig.name.list.txt perl $HOME/scripts/C_context_window_SREedits.pl {1} 100 0 {2}

find -name "*bismark.txt" | xargs pigz