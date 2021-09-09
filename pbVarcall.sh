#!/bin/bash

INFILE=$1
# REFSEQ="LPL_cds.fa"
REFSEQ=$2
OUTFILE=$(echo "$INFILE"|sed -r "s/\\.txt.gz$/_varcalls.txt/")

mkdir -p chunks
zcat "$INFILE"|split -l 1000 - chunks/chunk
for CHUNK in $(ls chunks/*); do
    ID=$(basename $CHUNK)
    submitjob.sh -t 00:30:00 -n "$ID" -- ./pbVarcall_worker.sh "$CHUNK" "$REFSEQ"
done

waitForJobs.sh

cat chunks/*_varcall.txt>$OUTFILE&&rm chunks/*

submitjob.sh -t 00:30:00 -n translate -- pbVarcall_translate.R $OUTFILE parameters.json

waitForJobs.sh
