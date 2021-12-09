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

NCHUNKS=1000
BLACKLIST=""

#helper function to print usage information
usage () {
  cat << EOF

pbVarcall.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

Runs parallel variant calling on Pacbio consensus sequences using needle  
SLURM HPC cluster.
Usage: pacbioCCS.sh [-c|--chunks <CHUNKS>] [-b|--blacklist <BLACKLIST>] 
    <INFILE> <PARAMETERS>

-c|--chunks    : The number of chunks to process in parallel.
                 Defaults to $NCHUNKS .
-b|--blacklist : An optional comma-separated blacklist of nodes to avoid
<INFILE>       : The compressed, tab-delimited .txt.gz file 
                 containing Pacbio sequences
<PARAMETERS>   : A barseq parameter sheet file

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
    -c|--chunks)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        NCHUNKS=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
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

INFILE="$1"
if [ -z "$INFILE" ]; then
  echo "No INFILE provided!">&2
  usage 1
elif [[ "$INFILE" != *.txt.gz ]] || ! [[ $(file "$INFILE") =~ "gzip compressed data" ]]; then
  echo "$INFILE is not a gzip-compressed text file!">&2
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

REFSEQ=$(extractParamSection $PARAMETERS 'CODING SEQUENCE')

# REFSEQ="LPL_cds.fa"
# REFSEQ="$2"
# if [ -z "$REFSEQ" ]; then
#   echo "No REFERENCE FASTA provided!">&2
#   usage 1
# elif [[ "$REFSEQ" != *.fasta ]] && [[ "!REFSEQ" != *.fa ]]; then
#   echo "$REFSEQ must be a FASTA file!"
#   exit 1
# fi

#set outfile based on infile
OUTFILE=$(echo "$INFILE"|sed -r "s/\\.txt.gz$/_varcalls.txt/")

#split file into chunks and process using worker scripts
echo "Splitting file into chunks..."
mkdir -p chunks
zcat "$INFILE"|split -l $NCHUNKS - chunks/chunk

echo "Calling mutations in chunks..."
JOBS=""
for CHUNK in $(ls chunks/*); do
    ID=$(basename $CHUNK)
    RETVAL=$(submitjob.sh -t "00:30:00" -n "$ID" $BLARG -- \
      pbVarcall_worker.sh "$CHUNK" "$REFSEQ")
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

#consolidate outputs
echo "Consolidating results..."
cat chunks/*_varcall.txt>$OUTFILE&&rm chunks/*

#translate variant calls to amino acid level
echo "Translating results to amino acid levels..."
RETVAL=$(submitjob.sh -t "02:00:00" -c 8 -n translate $BLARG -- \
pbVarcall_translate.R "$OUTFILE" "$REFSEQ")
JOBID=${RETVAL##* }

waitForJobs.sh -v "$JOBID"

#cleanup 
rmdir chunks alignments

echo "Done!"
