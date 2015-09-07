#!/bin/bash

cd /home/steve/conviron_log_plot/

#go grab all the log csv files from the large-date
for FILE in /network/phenocam-largedatasets/a_data/_webroot/NCRIS-Chambers/chamber_logs/*.csv
do
tail -n 1250 $FILE > ${FILE}_recent.chamberlog
echo "1250 lines from ${FILE} written to chamberlog file"
cp ${FILE}_recent.chamberlog .
done

#create html file
R -e "library(rmarkdown); render('status.Rmd')"

#move html file to large-data webspace
mv status.html /network/phenocam-largedatasets/a_data/_webroot/NCRIS-Chambers/
