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

#fail on error, even within pipes; disallow use of unset variables, enable history tracking
set -euo pipefail -o history -o histexpand

BARCODE=SWSWSWSWSWSWSWSWSWSWSWSWS
BCPOS="153,2812"
ORFSTART=207
ORFEND=2789
THREADS=24


#helper function to print usage information
usage () {
  cat << EOF

pacifica_worker.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

Runs Pacifica worker
Usage: pacifica_worker.sh [-b|--barcode <BARCODE>] [-p|--barcodePos <BCPOS>]
   [-s|--orfStart <ORFSTART>] [-e|--orfEnd <ORFEND>] [-c|--cpus <NUMBER>]
   <FASTQ> <FASTA> <WORKSPACE>

-b|--barcode   : The barcode degeneracy code sequence, defaults to
                 $BARCODE
-p|--barcodePos: Comma-separated list of barcode position, defaults to $BCPOS
-s|--orfStart  : The ORF start position, defaults to $ORFSTART
-e|--orfEnd    : The ORF end position, defaults to $ORFEND
-c|--cpus      : Number of CPU cores, defaults to $THREADS
<FASTQ>        : The input fastq.gz file
<FASTA>        : The reference fasta file (with removed barcodes)
<WORKSPACE>    : The workspace directory

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
    -p|--barcodePos)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        # if ! [[ $2 =~ $NUMRX ]] ; then
        #    echo "ERROR: barcode position must be a positive integer number" >&2
        #    usage 1
        # fi
        BCPOS=$2
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


#derive barcode length from barcode string
BLEN=${#BARCODE}

#validate positional arguments
INFQ="$1"
if ! [[ -r "$INFQ" ]]; then
  echo "ERROR: Unable to read FASTQ file $INFQ">&2
  exit 1
fi
RX='\.fastq$'
if ! [[ "$INFQ" =~ $RX ]]; then
   echo 'ERROR: First parameter must be a *.fastq file' >&2
   exit 1
fi

REFFASTANOBC="$2"
if ! [[ -r "$REFFASTANOBC" ]]; then
  echo "ERROR: Unable to read FASTA file $REFFASTANOBC">&2
  exit 1
fi
RX='\.fasta|\.fa$'
if ! [[ "$REFFASTANOBC" =~ $RX ]]; then
   echo 'ERROR: Second parameter must be a *.fasta file' >&2
   exit 1
fi

WORKSPACE="$3"
if ! [[ -d "$WORKSPACE" ]]; then
  echo "ERROR: $WORKSPACE is not a directory">&2
  exit 1
fi
RX='/$'
if ! [[ "$WORKSPACE" =~ $RX ]]; then
  WORKSPACE="${WORKSPACE}/"
fi
echo $WORKSPACE



#define intermediate files
CHUNKPREFIX=$(basename $INFQ|sed -r "s/.fastq$//")
#alignment output
ALNFILE="${WORKSPACE}/${CHUNKPREFIX}_aligned.bam"
#extracted barcodes file
EXTRACTDIR="${WORKSPACE}/${CHUNKPREFIX}_extract/"

#create the directory
mkdir -p $EXTRACTDIR

#align to template
if [[ ! -s "$ALNFILE" ]]; then
  echo "Running alignment..."
  # bwa mem -t $THREADS -C -M -L 80 "$REFFASTANOBC" $INFQ | samtools view -u - \
  #   | samtools sort -o "$ALNFILE" - 
  bwa mem -t $THREADS -C -M -L 80 "$REFFASTANOBC" $INFQ | samtools view -b -o "$ALNFILE" - 
  # dieOnError $? !!
else
  echo "Using existing alignment"
fi

# #generate alignment stats
# if [[ ! -s "${ALNFILE}_stats" ]]; then
#   echo "Generating alignment stats..."
#   samtools flagstat "$ALNFILE" > "${ALNFILE}_stats"
#   #inspect alignment quality
#   samtools view "$ALNFILE" |awk '{ if ($4 == 1) print $6 }' \
#     |sort |uniq -c |sort -nr |head >"${ALNFILE}_topCIGARs" || true
# else
#   echo "Using existing stats"
# fi

#extract barcodes
if [[ ! -s "${EXTRACTDIR}/bcExtract_1.fastq.gz" ]]; then
  echo "Extracting barcodes..."
  pacifica_extractBCORF.R <(samtools view "$ALNFILE") "$REFFASTANOBC" "$EXTRACTDIR"\
    --bcLen $BLEN --bcPos $BCPOS --orfStart $ORFSTART --orfEnd $ORFEND
  #calculate barcode length distributions
  # zcat "${EXTRACTDIR}/bcExtract_1.fastq.gz"|grep len=|cut -f 3,3 -d'='|\
  #   sort -n|uniq -c>"${EXTRACTDIR}/bc1len_distr.txt"
  # zcat "${EXTRACTDIR}/bcExtract_combo.fastq.gz"|grep len=|cut -f 3,3 -d'='|\
  #   sort -n|uniq -c>"${EXTRACTDIR}/bccombolen_distr.txt"
else 
  echo "Using existing extracted barcodes"
fi

echo "Done!"