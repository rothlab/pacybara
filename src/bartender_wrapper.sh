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

#helper function to print usage information
DIRECTION=f
MAXERR=2
QUALITY=?
THREADS=8
WORKSPACE=./

usage () {
  cat << EOF

bartender_wrapper.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

A convenient wrapper for Bartender
Usage: bartender_wrapper.sh [-d|--direction {f|rc}] [-p|--pattern <PATTERN>] 
       [-m|--maxErr <INT>] [-q|--quality <CHAR>] [-t|--threads <INT>] 
       [-w|--workspace <DIR>] <FASTQ>

<FASTQ>        : The fastq.gz file to process
-d|--direction : Either f for forward or rc for reverse-complement. Default: $DIRECTION
-p|--pattern   : The extraction pattern (see bartender help for details), 
                 e.g. ACGGG[23-27]CTAAC
-m|--maxErr    : The maximum number of errors allowed (for barcode extraction 
                 and cluster merging). Default: $MAXERR
-q|--quality   : The minimum PHRED quality ASCII code. Default: $QUALITY
-t|--threads   : The number of parallel threads to run. Default: $THREADS
-w|--workspace : The workspace directory. Default: $WORKSPACE

EOF
 exit $1
}

#Parse Arguments
PARAMS=""
BLACKLIST=""
DEBUG=0
while (( "$#" )); do
  case "$1" in
    -h|--help)
      usage 0
      shift
      ;;
    -d|--direction)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        DIRECTION=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -p|--pattern)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        PATTERN=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -m|--maxErr)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        MAXERR=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -q|--quality)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        QUALITY=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -t|--threads)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        THREADS=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -w|--workspace)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        WORKSPACE=$2
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

FQGZ=$1
#prefix for intermediate and output files
PREFIX=${WORKSPACE}/$(basename ${FQGZ%.fastq.gz})

#unfortunately bartender can't deal with the <(gzip -dc $FQGZ) file handle
#so we'll have to extract the gz file first for bartender to be able to use it
FQ=$(mktemp --suffix=.fq)
gzip -dc $FQGZ>$FQ
#extract barcodes
bartender_extractor_com -f $FQ -o "$PREFIX" \
    -q $QUALITY -p "${PATTERN}" -m $MAXERR -d $DIRECTION
rm $FQ

#run clustering
bartender_single_com -f ${PREFIX}_barcode.txt -o ${PREFIX} \
    -d $MAXERR -t $THREADS

#compress intermediate results
gzip "${PREFIX}_barcode."*
# gzip ${PREFIX}_barcode.csv
mv "${PREFIX}_barcode."* "${WORKSPACE}/extract/"
mv "${PREFIX}_cluster.csv" "${PREFIX}_quality.csv" "${WORKSPACE}/counts/"

echo "Done!"
