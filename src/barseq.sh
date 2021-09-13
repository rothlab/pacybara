#!/bin/bash
# Copyright (C) 2021  Jochen Weile, Roth Lab
#
# This file is part of BarseqPro.
#
# BarseqPro is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BarseqPro is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with BarseqPro.  If not, see <https://www.gnu.org/licenses/>.

#This script depends on submitjob.sh and waitForJobs.sh !

#helper function to print usage information
usage () {
  cat << EOF

barseq.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

Processes BarSeq data ona SLURM HPC cluster.
Usage: barseq.sh <INDIR> <PARAMS>

<INDIR>     : The input directory containing the fastq.gz files
<PARAMS>    : A barseq parameter sheet file

EOF
 exit $1
}

#Parse Arguments
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      usage 0
      shift
      ;;
    --) # end of options indicates that the main command follows
      shift
      PARAMS="$PARAMS $@"
      eval set -- ""
      ;;
    -*|--*=) # unsupported flags
      echo "ERROR: Unsupported flag $1" >&2
      usage 1
      ;;
    *) # positional parameter
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
#reset command arguments as only positional parameters
eval set -- "$PARAMS"

#folder containing barseq fastq files
INPUTFOLDER="$1"
if [[ -z "$INPUTFOLDER" ]]; then
  echo "No INDIR provided!"
elif ! [[ -d "$INPUTFOLDER" ]]; then
  echo "$INPUTFOLDER is not a valid directory!">&2
  exit 1
fi
#parameter file
PARAMETERS=${2:-parameters.txt}
if ! [[ -r "$PARAMETERS" ]]; then
  echo "$PARAMETERS not found or unreadable!">&2
  exit 1
fi

#blacklist of nodes with broken R installations
# BLACKLIST=galen1,galen2,galen4,galen5,galen10,galen24,galen25,galen27,galen31,galen35,galen40,galen41,galen43,galen54,galen57,galen58,galen59,galen60,galen61,galen62,galen63,galen64,galen65

#helper function to extract relevant sections from a parameter file
extractParamSection() {
  INFILE="$1"
  SECTION="$2"
  mkdir -p tmp
  case "$SECTION" in
    ARGUMENTS) TMPFILE=$(mktemp -p tmp/);;
    *SEQUENCE*) TMPFILE=$(mktemp -p tmp/ --suffix=.fasta);;
    SAMPLE) TMPFILE=$(mktemp -p tmp/ --suffix=.tsv);;
    *) echo "Unrecognized section selected!"&&exit 2;;
  esac
  RANGE=($(grep -n "$SECTION" "$INFILE"|cut -f1 -d:))
  sed -n "$((${RANGE[0]}+1)),$((${RANGE[1]}-1))p;$((${RANGE[1]}))q" "$INFILE">"$TMPFILE"
  echo "$TMPFILE"
}

#Load parameters
source $(extractParamSection $PARAMETERS ARGUMENTS)
FLANKING=$(extractParamSection $PARAMETERS 'FLANKING SEQUENCES')
SAMPLES=$(extractParamSection $PARAMETERS SAMPLE)
CDS=$(extractParamSection $PARAMETERS 'CODING SEQUENCE')

#set the correct argument for reverse complement
if [[ $REVCOMP == 1 ]]; then
  RCARG="--rc"
else
  RCARG=""
fi

#Create output directory for this run
WORKSPACE=${TITLE}_$(date +%Y%m%d_%H%M%S)/
mkdir "$WORKSPACE"
cp "$PARAMETERS" "$WORKSPACE"

#create folders for logs and temporary chunk files
mkdir -p ${WORKSPACE}logs
mkdir -p ${WORKSPACE}chunks
mkdir -p ${WORKSPACE}counts
mkdir -p ${WORKSPACE}scores

for INPUTFQ in $(ls $INPUTFOLDER/*fastq.gz); do

  #prefix for input chunks
  PREFIX=$(basename "$INPUTFQ"|sed -r "s/\\.fastq\\.gz$/_/")
  #output file
  OUTFILE=${WORKSPACE}counts/$(basename "$INPUTFQ"|sed -r "s/\\.fastq\\.gz$/_counts.txt/")


  echo "Splitting FASTQ file into chunks"
  zcat "$INPUTFQ"|split -a 3 -l 200000 --additional-suffix .fastq - "${WORKSPACE}chunks/$PREFIX"
  CHUNKS=$(ls ${WORKSPACE}chunks/${PREFIX}*.fastq)


  JOBS=""
  echo "Processing chunks on HPC cluster."
  for CHUNK in $CHUNKS; do
    TAG=barseq$(echo $CHUNK|sed -r "s/.*_|\\.fastq//g")
    #start barseq.R job and capture the job-id number
      # --blacklist $BLACKLIST \
    RETVAL=$(submitjob.sh -n $TAG -l ${WORKSPACE}logs/${TAG}.log -t 24:00:00 \
      barseq_caller.R $RCARG --flanking $FLANKING --bcLen $BCLEN \
      --maxErr $BCMAXERR $CHUNK $LIBRARY)
    JOBID=${RETVAL##* }
    if [ -z "$JOBS" ]; then
      #if jobs is empty, set it to the new ID
      JOBS=$JOBID
    else
      #otherwise append the id to the list
      JOBS=${JOBS},$JOBID
    fi
  done

  waitForJobs.sh -v "$JOBS"

  #Validate outputs
  for CHUNK in $CHUNKS; do
    TAG=$(echo $CHUNK|sed -r "s/.*_|\\.fastq//g")
    LOG=${WORKSPACE}logs/barseq${TAG}.log
    STATUS=$(tail -1 $LOG)
    if [ "$STATUS" != "Done!" ]; then
      echo "$TAG failed!"
    fi
  done

  echo "Consolidating results..."
  RESULTS=$(ls ${WORKSPACE}chunks/${PREFIX}*_hits.csv.gz)
  zcat $RESULTS|cut -f 1 -d,|sort|uniq -c>$OUTFILE

  #clean up 
  rm $CHUNKS $RESULTS
  rm ${WORKSPACE}logs/*

done

#consolidates all OUTFILEs (*_counts.txt) into "allCounts.csv"
barseq_consolidator.R "${WORKSPACE}counts/" "$SAMPLES" "$LIBRARY"

barseq_enrichment.R "${WORKSPACE}counts/allCounts.csv" "$SAMPLES" "${WORKSPACE}scores/"

#cleanup temp files
rm tmp/*
