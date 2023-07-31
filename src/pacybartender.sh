#!/bin/bash
# Copyright (C) 2021, 2022  Jochen Weile, The Roth Lab
#
# This file is part of Pacybara.
#
# Pacybara is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Pacybara is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Pacybara.  If not, see <https://www.gnu.org/licenses/>.

#This script depends on submitjob.sh and waitForJobs.sh !

#helper function to print usage information
usage () {
  cat << EOF

pacybartender.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

A convenient wrapper for Bartender that performs lookups against Pacybara
libraries and performs downstream analysis
Usage: pacybartender.sh [-b|--blacklist <BLACKLIST>] [-c|--conda <ENV>] <INDIR> <PARAMS>

<INDIR>        : The input directory containing the fastq.gz files
<PARAMS>       : A barseq parameter sheet file
-b|--blacklist : An optional comma-separated blacklist of nodes to avoid
-c|--conda     : Conda environment to activate for bartender jobs (for python2.7)

EOF
 exit $1
}

#Parse Arguments
PARAMS=""
BLACKLIST=""
CONDAARG=""
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
    -c|--conda)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        CONDAARG="--conda $2"
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

TMPDIR=$(mktemp -d -p ./)

#helper function to extract relevant sections from a parameter file
extractParamSection() {
  INFILE="$1"
  SECTION="$2"
  TMPDIR="$3"
  # mkdir -p tmp
  case "$SECTION" in
    ARGUMENTS) TMPFILE=$(mktemp -p "$TMPDIR");;
    *SEQUENCE*) TMPFILE=$(mktemp -p "$TMPDIR" --suffix=.fasta);;
    SAMPLE) TMPFILE=$(mktemp -p "$TMPDIR" --suffix=.tsv);;
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
source $(extractParamSection $PARAMETERS ARGUMENTS "$TMPDIR")
FLANKING=$(extractParamSection $PARAMETERS 'FLANKING SEQUENCES' "$TMPDIR")
SAMPLES=$(extractParamSection $PARAMETERS SAMPLE "$TMPDIR")
CDS=$(extractParamSection $PARAMETERS 'CODING SEQUENCE' "$TMPDIR")

#validate parameters
if [[ ! -r $LIBRARY ]]; then
  echo "Library at $LIBRARY could not be found or read!">&2
  exit 1
fi

#bartender doesn't accept flanking sequences longer than 5bp
#so we trunkate them accordingly here
FLINES=(`cat $FLANKING`)
UPSTREAM=${FLINES[1]:(-5):5}
DOWNSTREAM=${FLINES[3]:0:5}
#and use them to construct the bartender extraction pattern
PATTERN="${UPSTREAM}[$((BCLEN-BCMAXERR))-$((BCLEN+BCMAXERR))]${DOWNSTREAM}"

#Create output directory for this run
WORKSPACE=${TITLE}_$(date +%Y%m%d_%H%M%S)/
mkdir "$WORKSPACE"
#make convenience copies of the parameter sheet and sample tables
cp "$PARAMETERS" "$WORKSPACE"
cp "$SAMPLES" "$WORKSPACE/samples.txt"

#create folders for logs and results
mkdir -p ${WORKSPACE}logs
mkdir -p ${WORKSPACE}extract
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
  done
  echo "$R1FQS"
}

processChunks() {

  CHUNKS="$1"

  JOBS=""
  echo "Processing chunks on HPC cluster."
  for CHUNK in $CHUNKS; do
    # TAG=barseq$(echo $CHUNK|sed -r "s/.*_|\\.fastq//g")
    TAG="bt$(basename ${CHUNK%.fastq*})"

    #set the correct argument for reverse complement
    if [[ $REVCOMP == 1 ]]; then
      RCARG="rc"
    else
      RCARG="f"
    fi

    #start barseq.R job and capture the job-id number
    RETVAL=$(submitjob.sh -n $TAG -l ${WORKSPACE}logs/${TAG}.log \
      -e ${WORKSPACE}logs/${TAG}.log -t 24:00:00 $BLARG \
      -m 16G -c 8 $CONDAARG \
      bartender_wrapper.sh -d $RCARG -p "$PATTERN" -m $BCMAXERR \
      -q "?" -t 8 -w ${WORKSPACE}/ $CHUNK )
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
    # TAG=$(echo $CHUNK|sed -r "s/.*_|\\.fastq//g")
    TAG="bt$(basename ${CHUNK%.fastq*})"
    LOG=${WORKSPACE}logs/${TAG}.log
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


R1FQS="$(validateFASTQs)"
processChunks "$R1FQS"
echo "Validating jobs..."
FAILEDCHUNKS=$(checkForFailedJobs "$R1FQS")

#If any jobs failed, try 2 more times
TRIES=1
while ! [[ -z $FAILEDCHUNKS ]]; do
  if [[ $TRIES < 3 ]]; then
    echo "Attempting to re-run failed jobs."
    processChunks "$FAILEDCHUNKS"
    echo "Validating jobs..."
    FAILEDCHUNKS=$(checkForFailedJobs "$R1FQS")
    ((TRIES++))
  else
    echo "ERROR: Exhausted 3 attempts at re-running failed jobs!">&2
    exit 1
  fi
done

echo "Cleaning up intermediate files..."
# mv "${WORKSPACE}/*_barcode.*.gz" ${WORKSPACE}/extract/
# mv "${WORKSPACE}/*_cluster.csv" ${WORKSPACE}/counts/
tar czf ${WORKSPACE}/counts/qualityMatrices.tgz ${WORKSPACE}/counts/*_quality.csv&&rm ${WORKSPACE}/counts/*_quality.csv

echo "Consolidating counts from all samples..."
bartender_consolidator.R "${WORKSPACE}counts/" "$SAMPLES" "$LIBRARY"

#compress intermediate files
echo "Compressing files..."
gzip ${WORKSPACE}/counts/*_cluster.csv


#assemble parameter list for bartender_combiner
# COMBOLIST=""
# for R1FQ in $R1FQS; do
#   PREFIX=${WORKSPACE}/$(basename ${R1FQ%.fastq.gz})
#   if [[ -z "$COMBOLIST" ]]; then
#     COMBOLIST="${PREFIX}_cluster.csv,${PREFIX}_quality.csv"
#   else
#     COMBOLIST="${COMBOLIST},${PREFIX}_cluster.csv,${PREFIX}_quality.csv"
#   fi
# done
# #run bartender combiner to consolidate results
# bartender_combiner_com -f "$COMBOLIST" -c 1 -o "${WORKSPACE}/counts/rawCounts"

if [[ -n $BNFILTER ]]; then
  BNARG="--bnFilter $BNFILTER"
else
  BNARG=""
fi
if [[ -n $FREQFILTER ]]; then
  FFARG="--freqFilter $FREQFILTER"
else 
  FFARG=""
fi

echo "Calculating enrichment ratios and scores..."
barseq_enrichment.R "${WORKSPACE}counts/allCounts.csv" "$SAMPLES" "${WORKSPACE}scores/" $FFARG $BNARG

echo "Running QC"
#FIXME: Need to process different QC
bartender_qc.R "${WORKSPACE}scores/allLRs.csv" "${WORKSPACE}counts/allCounts.csv" "$SAMPLES" "${WORKSPACE}qc/" --logfolder "${WORKSPACE}logs/" $FFARG $BNARG

#cleanup temp files
rm -r "$TMPDIR"

echo "Done!"
