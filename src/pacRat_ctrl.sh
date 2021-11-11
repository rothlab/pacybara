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

#fail on error, even within pipes; disallow use of unset variables
# set -euo pipefail 

BARCODE=SWSWSWSWSWSWSWSWSWSWSWSWS

# BCPOS=153
ORFSTART=207
ORFEND=2789
THREADS=24
QARG=""

#helper function to print usage information
usage () {
  cat << EOF

pacRat_ctrl.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

Runs PacbioSubassembly and PacRat
Usage: pacRat_ctrl.sh [-b|--barcode <BARCODE>] [-s|--orfStart <ORFSTART>] [-q|--queue <QUEUE>]
   [-e|--orfEnd <ORFEND>] [-c|--cpus <NUMBER>] <INDIR> <FASTA> [<WORKSPACE>]

-b|--barcode   : The barcode degeneracy code sequence, defaults to
                 $BARCODE
-c|--cpus      : Number of CPUs to use, defaults to $THREADS
-s|--orfStart  : The ORF start position, defaults to $ORFSTART
-e|--orfEnd    : The ORF end position, defaults to $ORFEND
-q|--queue     : The queue (slurm partition) to use
<INDIR>        : The input directory containing fastq.gz files to process
<FASTA>        : The raw reference fasta file 
<WORKSPACE>    : The workspace directory. Defaults to 'workspace/'

EOF
 exit $1
}


#Parse Arguments
NUMRX='^[0-9]+$'
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      usage 0
      shift
      ;;
    -b|--barcode)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        BARCODE=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -c|--cpus)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        if ! [[ $2 =~ $NUMRX ]] ; then
           echo "ERROR: cpus must be a positive integer number" >&2
           usage 1
        fi
        THREADS=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -s|--orfStart)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        if ! [[ $2 =~ $NUMRX ]] ; then
           echo "ERROR: orfStart must be a positive integer number" >&2
           usage 1
        fi
        ORFSTART=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -e|--orfEnd)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        if ! [[ $2 =~ $NUMRX ]] ; then
           echo "ERROR: orfEnd must be a positive integer number" >&2
           usage 1
        fi
        ORFEND=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -q|--queue)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        QARG="-q $2"
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    --) # end of options indicates that only positional arguments follow
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

#validate positional arguments
INDIR=$1
if ! [[ -d "$INDIR" ]]; then
  echo "ERROR: $INDIR is not a directory!">&2
  exit 1
fi
#add trailing "/" to directory name, if missing
RX='/$'
if ! [[ "$INDIR" =~ $RX ]]; then
  INDIR="${INDIR}/"
fi
REFFASTA=$2
if ! [[ -r "$REFFASTA" ]]; then
  echo "ERROR: Unable to read FASTA file $REFFASTA">&2
  exit 1
fi
RX='\.fasta|\.fa$'
if ! [[ "$REFFASTA" =~ $RX ]]; then
   echo 'ERROR: First parameter must be a *.fasta file' >&2
   exit 1
fi
WORKSPACE=${3:-workspace/}
if ! [[ -d "$WORKSPACE" ]]; then
  echo "ERROR: $WORKSPACE is not a directory!">&2
  exit 1
fi
RX='/$'
if ! [[ "$WORKSPACE" =~ $RX ]]; then
  WORKSPACE="${WORKSPACE}/"
fi


function removeBarcode() {
  INFASTA=$1
  BARCODE=${2:-SWSWSWSWSWSWSWSWSWSWSWSWS}
  #read FASTA file
  mapfile -t FLINES<"$INFASTA"
  FHEADER=${FLINES[0]}
  #ensure all bases are upper case and there are no whitespaces
  FBODY=$(echo ${FLINES[*]:1}|tr [a-z] [A-Z]|tr -d ' ')
  #remove BARCODES
  FBODY2=$(echo $FBODY|sed -r "s/${BARCODE}//g")
  #write to output file
  OUTFASTA=$(echo $INFASTA|sed -r "s/.fa$/_noBC.fa/")
  echo $FHEADER>$OUTFASTA
  printf "$FBODY2\n"|fold -w 80 >>$OUTFASTA
  #return the name of the new FASTAFILE
  echo $OUTFASTA
}

function findBarcodPos() {
  INFASTA=$1
  BARCODE=${2:-SWSWSWSWSWSWSWSWSWSWSWSWS}
  #read FASTA file
  mapfile -t FLINES<"$INFASTA"
  FHEADER=${FLINES[0]}
  #ensure all bases are upper case and there are no whitespaces
  FBODY=$(echo ${FLINES[*]:1}|tr [a-z] [A-Z]|tr -d ' ')
  #cut off the barcode and anything following it
  PREFIX=${FBODY%%$BARCODE*}
  #then return the length of the remaining prefix
  echo $((${#PREFIX}+1))
}

#find barcode position in reference
BCPOS=$(findBarcodPos $REFFASTA $BARCODE)
#remove barcode from reference
REFFASTANOBC=$(removeBarcode $REFFASTA $BARCODE)
#and create an index for the barcodeless reference
samtools faidx "$REFFASTANOBC"
bwa index -a is "$REFFASTANOBC"


mkdir -p ${WORKSPACE}alignments
mkdir -p ${WORKSPACE}consensus
mkdir -p ${WORKSPACE}logs

#iterate over all input FASTQ files
i=0
INFQS=$(ls ${INDIR}*fastq.gz)
for INFQ in $INFQS; do
  
  #set log file
  LOGFILE=${WORKSPACE}logs/$(basename $INFQ|sed -r "s/.fastq.gz$/.log/")
  #And start worker jobs
  ((i++))
  echo "Scheduling job #$i for $INFQ"
  RETVAL=$(submitjob.sh -n "pacrat$i" -t 24:00:00 $QARG\
    -c $THREADS -m $(($THREADS*2))G -l "$LOGFILE" -e "${LOGFILE}err" -- \
    scripts/pacRat_worker.sh --barcode $BARCODE --barcodePos $BCPOS \
    --orfStart $ORFSTART --orfEnd $ORFEND \
    "$INFQ" "$REFFASTANOBC" "$WORKSPACE")
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

echo "Done!"