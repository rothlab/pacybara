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
Usage: barseq.sh [-b|--blacklist <BLACKLIST>] <INDIR> <PARAMS>

<INDIR>        : The input directory containing the fastq.gz files
<PARAMS>       : A barseq parameter sheet file
-b|--blacklist : An optional comma-separated blacklist of nodes to avoid

EOF
 exit $1
}

#Parse Arguments
PARAMS=""
BLACKLIST=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      usage 0
      shift
      ;;
    -b|--blacklist)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        BLACKLIST=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
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

if [[ -z $BLACKLIST ]]; then
  BLARG=""
else
  BLARG="--blacklist $BLACKLIST"
fi

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

#helper function to extract string using regex
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

#TODO:
#Check that sample names in parameter sheet match fastq names!

#Create output directory for this run
WORKSPACE=${TITLE}_$(date +%Y%m%d_%H%M%S)/
mkdir "$WORKSPACE"
cp "$PARAMETERS" "$WORKSPACE"

#create folders for logs and temporary chunk files
mkdir -p ${WORKSPACE}logs
mkdir -p ${WORKSPACE}chunks
mkdir -p ${WORKSPACE}counts
mkdir -p ${WORKSPACE}scores

#helper function to check that FASTQ files exist for each sample
#listed in the parameter sheet
function validateFASTQs() {
  SAMPLENAMES=$(tail -n +2 ${SAMPLES}|cut -f 1)
  R1FQS=""
  for SAMPLENAME in $SAMPLENAMES; do
    FILEMATCH=$(echo ${INPUTFOLDER}/${SAMPLENAME}*R1*.fastq.gz)
    if [[ -r $FILEMATCH ]]; then
      R1FQ=$(ls ${INPUTFOLDER}/${SAMPLENAME}*R1*.fastq.gz)
      R1FQS="${R1FQS} ${R1FQ}"
    else 
      echo "ERROR: Unable to find or read FASTQ file for sample $SAMPLENAME">&2
      exit 2
    fi
    #in case of paired-end mode, also check for R2 file
    if [[ $PAIREDEND == 1 ]]; then
      R2FQ=$(echo "$R1FQ"|sed -r "s/_R1_/_R2_/")
      if  [[ ! -r "$R2FQ" ]]; then
        echo "ERROR: Unable to find or read R2 FASTQ file $R2FQ !">&2
        exit 1
      fi
    fi
  done
  echo "$R1FQS"
}
R1FQS="$(validateFASTQs)"
#the code block below was superceded by the new function above
# #if we're in paired-end mode, look for R1 and R2 files separately
# if [[ $PAIREDEND == 1 ]]; then
#   R1FQS=$(ls $INPUTFOLDER/*_R1_*fastq.gz)
#   #Do a preliminary scan to see if all files are accounted for
#   for R1FQ in $R1FQS; do
#     R2FQ=$(echo "$R1FQ"|sed -r "s/_R1_/_R2_/")
#     if  [[ ! -r "$R2FQ" ]]; then
#       echo "ERROR: Unable to find or read R2 file $R2FQ !">&2
#       exit 1
#     fi
#   done
# else
#   #otherwise just use all fastq files directly
#   R1FQS=$(ls $INPUTFOLDER/*fastq.gz)
# fi

#helper function to process a list of FASTQ chunks
processChunks() {

  CHUNKS=$1

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
    RETVAL=$(submitjob.sh -n $TAG -l ${WORKSPACE}logs/${TAG}.log \
      -e ${WORKSPACE}logs/${TAG}.log -t 24:00:00 $BLARG \
      -m 12G -c 4 \
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

}

#helper function to check for failed jobs.
#returns a list of any failed jobs
checkForFailedJobs() {

  CHUNKS=$1

  FAILEDCHUNKS=""
  for CHUNK in $CHUNKS; do
    TAG=$(echo $CHUNK|sed -r "s/.*_|\\.fastq//g")
    LOG=${WORKSPACE}logs/barseq${TAG}.log
    STATUS=$(tail -1 $LOG)
    if [ "$STATUS" != "Done!" ]; then
      echo "Process $TAG failed!">&2
      if [[ -z $FAILEDCHUNKS ]]; then
        FAILEDCHUNKS="$CHUNK"
      else
        FAILEDCHUNKS="$FAILEDCHUNKS $CHUNK"
      fi
    fi
  done
  echo "$FAILEDCHUNKS"
}

for R1FQ in $R1FQS; do

  echo ""
  echo "Processing $R1FQ"

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

  #process chunks and wait for completion
  processChunks "$CHUNKS"
  echo "Validating jobs..."
  FAILEDCHUNKS=$(checkForFailedJobs "$CHUNKS")

  #If any jobs failed, try 2 more times
  TRIES=1
  while ! [[ -z $FAILEDCHUNKS ]]; do
    if [[ $TRIES < 3 ]]; then
      echo "Attempting to re-run failed jobs."
      processChunks "$FAILEDCHUNKS"
      echo "Validating jobs..."
      FAILEDCHUNKS=$(checkForFailedJobs "$CHUNKS")
      ((TRIES++))
    else
      echo "ERROR: Exhausted 3 attempts at re-running failed jobs!">&2
      exit 1
    fi
  done


  #Consolidate exception counts
  FAILEDEXTRACTION=0
  NOMATCH=0
  AMBIGUOUS=0
  for CHUNK in $CHUNKS; do
    TAG=$(echo $CHUNK|sed -r "s/.*_|\\.fastq//g")
    LOG=${WORKSPACE}logs/barseq${TAG}.log
    EXCLINE=$(grep failedExtraction $LOG)
    ((FAILEDEXTRACTION += $(extractRX $EXCLINE "failedExtraction=([0-9]+)") ))
    ((NOMATCH += $(extractRX $EXCLINE "noMatch=([0-9]+)") ))
    ((AMBIGUOUS += $(extractRX $EXCLINE "ambiguous=([0-9]+)") ))
  done
  echo "$R1PREFIX">>${WORKSPACE}/counts/exceptions.txt
  echo "failedExtraction=$FAILEDEXTRACTION">>${WORKSPACE}/counts/exceptions.txt
  echo "noMatch=$NOMATCH">>${WORKSPACE}/counts/exceptions.txt
  echo "ambiguous=$AMBIGUOUS">>${WORKSPACE}/counts/exceptions.txt

  #Consolidate individual result chunks
  echo "Consolidating results..."
  RESULTS=$(ls ${WORKSPACE}chunks/${R1PREFIX}*_hits.csv.gz)
  zcat $RESULTS|cut -f 1 -d,|sort|uniq -c>$OUTFILE

  #clean up 
  rm $CHUNKS $RESULTS
  tar czf "${WORKSPACE}logs/countLogs_${R1PREFIX}.tgz" ${WORKSPACE}logs/*log
  rm ${WORKSPACE}logs/*log

done

#consolidates all OUTFILEs (*_counts.txt) into "allCounts.csv"
echo "Consolidating counts from all samples..."
barseq_consolidator.R "${WORKSPACE}counts/" "$SAMPLES" "$LIBRARY"

echo "Calculating enrichment ratios and scores..."
barseq_enrichment.R "${WORKSPACE}counts/allCounts.csv" "$SAMPLES" "${WORKSPACE}scores/"

#cleanup temp files
rm tmp/*&&rmdir tmp

echo "Done!"
