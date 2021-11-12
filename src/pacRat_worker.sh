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
BCPOS=153
ORFSTART=207
ORFEND=2789
THREADS=24
MINREADS=2
MINFREQ=0.6

#helper function to print usage information
usage () {
  cat << EOF

pacRat_worker.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

Runs PacbioSubassembly and PacRat
Usage: pacRat_worker.sh [-b|--barcode <BARCODE>] [-p|--barcodePos <BCPOS>]
   [-s|--orfStart <ORFSTART>] [-e|--orfEnd <ORFEND>] [-c|--cpus <NUMBER>]
   <FASTQ> <FASTA> <WORKSPACE>

-b|--barcode   : The barcode degeneracy code sequence, defaults to
                 $BARCODE
-p|--barcodePos: The barcode position, defaults to $BCPOS
-s|--orfStart  : The ORF start position, defaults to $ORFSTART
-e|--orfEnd    : The ORF end position, defaults to $ORFEND
-c|--cpus      : Number of CPU cores, defaults to $THREADS
-r|--minReads  : Minimum number of reads required to accept barcode, 
                 defaults to $MINREADS
-f|--minFreq   : Minimum frequency of alignment base for consensus,
                 defaults to $MINFREQ
<FASTQ>        : The input fastq.gz file
<FASTA>        : The reference fasta file (with removed barcodes)
<WORKSPACE>    : The workspace directory

EOF
 exit $1
}

# #helper function for more helpful error messages
# dieOnError() {
#     EXITCODE=$1
#     LASTCMD=${@:2}
#     if [ $EXITCODE -ne 0 ]; then
#         >&2 echo "ERROR: \"${LASTCMD}\" failed with exit code ${EXITCODE}."
#         exit $EXITCODE
#     fi
# }


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
        if ! [[ $2 =~ $NUMRX ]] ; then
           echo "ERROR: barcode position must be a positive integer number" >&2
           usage 1
        fi
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
    -r|--minReads)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        if ! [[ $2 =~ $NUMRX ]] ; then
           echo "ERROR: minReads must be a positive integer number" >&2
           usage 1
        fi
        MINREADS=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -f|--minFreq)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        if ! [[ $2 > 0 && $2 < 1 ]] ; then
           echo "ERROR: minFreq must be a number between 0 and 1" >&2
           usage 1
        fi
        MINFREQ=$2
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
RX='\.fastq\.gz$'
if ! [[ "$INFQ" =~ $RX ]]; then
   echo 'ERROR: First parameter must be a *.fastq.gz file' >&2
   exit 1
fi

REFFASTANOBC="$2"
if ! [[ -r "$REFFASTANOBC" ]]; then
  echo "ERROR: Unable to read FASTA file $REFFASTANOBC">&2
  exit 1
fi
RX='\.fasta|\.fa$'
if ! [[ "$REFFASTANOBC" =~ $RX ]]; then
   echo 'ERROR: First parameter must be a *.fasta file' >&2
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
OUTPREFIX=$(basename $INFQ|sed -r "s/.fastq.gz$//")
#alignment output
ALNFILE="${WORKSPACE}alignments/${OUTPREFIX}_aligned.bam"
#re-formatted BAM of alignment
ALNFILE2="${WORKSPACE}alignments/${OUTPREFIX}_aligned_cropped.bam"
#extracted barcodes file
EXTRACTFILE="${WORKSPACE}alignments/${OUTPREFIX}_extracted.tsv.gz"
#output of PacbioSubassembly
HIQFILE="${WORKSPACE}alignments/${OUTPREFIX}_HiQ.tsv.gz"
#final output folder
OUTFOLDER="${WORKSPACE}consensus/${OUTPREFIX}/"
mkdir -p "$OUTFOLDER"
#final output file
OUTFILE="${OUTPREFIX}_subassembly.txt"


#align to template
if [[ ! -s "$ALNFILE" ]]; then
  echo "Running alignment..."
  bwa mem -t $THREADS -C -M -L 80 "$REFFASTANOBC" $INFQ | samtools view -u - \
    | samtools sort -o "$ALNFILE" - 
  # dieOnError $? !!
else
  echo "Using existing alignment"
fi

#generate alignment stats
if [[ ! -s "${ALNFILE}_stats" ]]; then
  echo "Generating alignment stats..."
  samtools flagstat "$ALNFILE" > "${ALNFILE}_stats"
  #inspect alignment quality
  samtools view "$ALNFILE" |awk '{ if ($4 == 1) print $6 }' \
    |sort |uniq -c |sort -nr |head >"${ALNFILE}_topCIGARs"
else
  echo "Using existing stats"
fi

#remove non-standard BAM fields for compatibility with pysam and re-index
if [[ ! -s "${ALNFILE2}" ]]; then
  echo "Removing non-standard BAM fields..."
  samtools view -h "$ALNFILE"|cut -f -14 |samtools view -Sb - > "$ALNFILE2" \
    && samtools index -b "$ALNFILE2"
  # dieOnError $? !!
else
  echo "Using existing cropped alignment"
fi

#extract barcodes
if [[ ! -s $EXTRACTFILE ]]; then
  echo "Extracting barcodes..."
  python AssemblyByPacBio/extractBarcodeInsertPairs_moreQC.py "$ALNFILE2" \
    --length=$BLEN --position=$BCPOS --start=$(($ORFSTART-$BLEN-1)) --end=$(($ORFEND-$BLEN)) \
    | gzip -c >"$EXTRACTFILE"
  # dieOnError $? !!
else 
  echo "Using existing extracted barcodes"
fi

#generate "hiqual" file
if [[ ! -s $HIQFILE ]]; then
  echo "Running subassembly preprocessor..."
  python AssemblyByPacBio/unifyAssignment.py "$EXTRACTFILE" | gzip -c > "$HIQFILE"
  # dieOnError $? !!
else 
  echo "Using existing HiQ-file"
fi

#activate conda environment
# CONDA_PREFIX=$(conda info | grep -i 'base environment'|awk '{print $4}')
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate msa_ccs

#start PacRat
echo "Running PacRat..."
python3.6 PacRAT/msa_pacbio.py -d "$OUTFOLDER" -o "$OUTFILE" \
  --highQual <(gzip -cd "$HIQFILE") --inputSeqs <(gzip -cd "$EXTRACTFILE") \
  -c $MINREADS -t $MINFREQ -s -r -m $(which muscle) -n $(which needle)
# dieOnError $? !!

#compressing result file
gzip "${OUTFOLDER}$OUTFILE"

#run pbvarcall in a subshell within the output folder
(cd $OUTFOLDER && pbVarcall.sh "$OUTFILE" $PARAMETERS)

#deactivate environment
conda deactivate

echo "Done!"
