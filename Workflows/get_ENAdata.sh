#!/bin/bash
# GET ENA DATA USIGN FTP
DATA=${1}
for sample in `awk '{print $1}' $DATA`
  do
    echo -e "Going to download "$sample""
    fastq-dump "$sample"
done
