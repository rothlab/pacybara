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

extractRX() {
  if [[ $1 =~ $2 ]]; then
    echo ${BASH_REMATCH[1]}
  fi
}

#Load parameters
source $(extractParamSection $PARAMETERS ARGUMENTS)
FLANKING=$(extractParamSection $PARAMETERS 'FLANKING SEQUENCES')
SAMPLES=$(extractParamSection $PARAMETERS SAMPLE)
CDS=$(extractParamSection $PARAMETERS 'CODING SEQUENCE')

#validate parameters
if [[ ! -r $LIBRARY ]]; then
  echo "Library at $LIBRARY could not be found or read!">&2
  exit 1
fi

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

#if we're in paired-end mode, look for R1 and R2 files separately
if [[ $PAIREDEND == 1 ]]; then
  R1FQS=$(ls $INPUTFOLDER/*_R1_*fastq.gz)
  #Do a preliminary scan to see if all files are accounted for
  for R1FQ in $R1FQS; do
    R2FQ=$(echo "$R1FQ"|sed -r "s/_R1_/_R2_/")
    if  [[ ! -r "$R2FQ" ]]; then
      echo "ERROR: Unable to find or read R2 file $R2FQ !">&2
      exit 1
    fi
  done
else
  #otherwise just use all fastq files directly
  R1FQS=$(ls $INPUTFOLDER/*fastq.gz)
fi


for R1FQ in $R1FQS; do

  #output file
  OUTFILE=${WORKSPACE}counts/$(basename "$R1FQ"|sed -r "s/\\.fastq\\.gz$/_counts.txt/")

  #prefix for input chunks
  R1PREFIX=$(basename "$R1FQ"|sed -r "s/\\.fastq\\.gz$/_/")

  echo "Splitting FASTQ file into chunks"
  zcat "$R1FQ"|split -a 3 -l 200000 --additional-suffix .fastq - \
    "${WORKSPACE}chunks/$R1PREFIX"
  CHUNKS=$(ls ${WORKSPACE}chunks/${R1PREFIX}*.fastq)

  #deal with R2 reads
  if [[ $PAIREDEND == 1 ]]; then
    #infer name of R2 file again (we already checked for its existence above)
    R2FQ=$(echo "$R1FQ"|sed -r "s/_R1_/_R2_/")
    R2PREFIX=$(basename "$R2FQ"|sed -r "s/\\.fastq\\.gz$/_/")
    zcat "$R2FQ"|split -a 3 -l 200000 --additional-suffix .fastq - \
      "${WORKSPACE}chunks/$R2PREFIX"
  fi


  JOBS=""
  echo "Processing chunks on HPC cluster."
  for CHUNK in $CHUNKS; do
    TAG=barseq$(echo $CHUNK|sed -r "s/.*_|\\.fastq//g")

    if [[ $PAIREDEND == 1 ]]; then
      R2CHUNK=$(echo "$CHUNK"|sed -r "s/_R1_/_R2_/")
      R2FLAG="--r2 $R2CHUNK"
    else 
      R2FLAG=""
    fi

    #start barseq.R job and capture the job-id number
    RETVAL=$(submitjob.sh -n $TAG -l ${WORKSPACE}logs/${TAG}.log -t 24:00:00 \
      barseq_caller.R $RCARG --flanking $FLANKING --bcLen $BCLEN \
      --maxErr $BCMAXERR --r1 $CHUNK $R2FLAG $LIBRARY)
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
  echo "Validating jobs..."
  for CHUNK in $CHUNKS; do
    TAG=$(echo $CHUNK|sed -r "s/.*_|\\.fastq//g")
    LOG=${WORKSPACE}logs/barseq${TAG}.log
    STATUS=$(tail -1 $LOG)
    if [ "$STATUS" != "Done!" ]; then
      echo "$TAG failed!"
    fi
  done

  #Consolidate exception counts
  for CHUNK in $CHUNKS; do
    TAG=$(echo $CHUNK|sed -r "s/.*_|\\.fastq//g")
    LOG=${WORKSPACE}logs/barseq${TAG}.log
    EXCLINE=$(grep failedExtraction $LOG)
    ((FAILEDEXTRACTION += $(extractRX $EXCLINE "failedExtraction=([0-9]+)") ))
    ((NOMATCH += $(extractRX $EXCLINE "noMatch=([0-9]+)") ))
    ((AMBIGUOUS += $(extractRX $EXCLINE "ambiguous=([0-9]+)") ))
  done
  echo "failedExtraction=$FAILEDEXTRACTION">${WORKSPACE}/counts/exceptions.txt
  echo "noMatch=$NOMATCH">>${WORKSPACE}/counts/exceptions.txt
  echo "ambiguous=$AMBIGUOUS">>${WORKSPACE}/counts/exceptions.txt

  #Consolidate individual result chunks
  echo "Consolidating results..."
  RESULTS=$(ls ${WORKSPACE}chunks/${R1PREFIX}*_hits.csv.gz)
  zcat $RESULTS|cut -f 1 -d,|sort|uniq -c>$OUTFILE

  #clean up 
  rm $CHUNKS $RESULTS
  tar czf ${WORKSPACE}logs/countLogs.tgz ${WORKSPACE}logs/*log
  rm ${WORKSPACE}logs/*log

done

#consolidates all OUTFILEs (*_counts.txt) into "allCounts.csv"
barseq_consolidator.R "${WORKSPACE}counts/" "$SAMPLES" "$LIBRARY"

barseq_enrichment.R "${WORKSPACE}counts/allCounts.csv" "$SAMPLES" "${WORKSPACE}scores/"

#cleanup temp files
rm tmp/*&&rmdir tmp
