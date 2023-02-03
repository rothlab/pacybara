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

MAXERR=2
#helper function to print usage information
usage () {
  cat << EOF

barseq2_worker.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2022

Processes BarSeq data ona SLURM HPC cluster.
Usage: barseq2_worker.sh -r|--reference <FASTA> -b|--bcLen <INTEGER> 
      [-m|--maxErr <INTEGER>] -1|--r1 <STRING> [-2|--r2 <STRING>] 
      [-l|--logfile <FILE>] -y|--library <FILE>

-r|--reference : reference FASTA file with removed barcode sequence
-b|--bcLen     : barcode length
-m|--maxErr    : maxiumum allowed errors in barcode. Defaults to $MAXERR
-1|--r1        : forward read FASTQ file
-2|--r2        : reverse read FASTQ file (optional)
-l|--logfile   : log file name
-y|--library   : barcode library CSV file

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
    -r|--reference)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        REFFASTANOBC=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -b|--bcLen)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        BCLEN=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -b|--maxERR)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        MAXERR=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -1|--r1)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        R1FILE=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -2|--r2)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        R2FILE=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -l|--logfile)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        LOGFILE=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -y|--library)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        LIBRARYFILE=$2
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

CHUNKPREFIX=$(basename $R1FILE|sed -r "s/.fastq$//")
CHUNKDIR=$(dirname $R1FILE)
ALNFILE=${CHUNKDIR}/${CHUNKPREFIX}_refaln.bam

#STEP1: Align to reference
echo "Running alignment..."
if [[ -z $R2FILE ]]; then
  #Use SED to correct SAM format violation by BWA (barcode tag format)
  #BWA options:
  # -C            append FASTQ comments ; 
  # -M            mark shorter hits as secondary ; 
  # -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
  # -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
  # -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
  # -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
  # -B INT        penalty for a mismatch [4]
  # bwa mem -t $THREADS -C -M -L 80 "$REFFASTANOBC" $R1FILE \
  #  | sed -e 's/1:N:0:/BC:Z:/g' \
  #  | samtools view -b -o "$ALNFILE" - 
  bwa mem -t $THREADS -M -A 3 -O 6 -E 1 -L 80 "$REFFASTANOBC" $R1FILE \
   | samtools view -b -o "$ALNFILE" - 
else
  bwa mem -t $THREADS -M -L 80 "$REFFASTANOBC" $R1FILE $R2FILE\
   | samtools view -b -o "$ALNFILE" - 
fi
#STEP2: Extract barcode from alignment
pacybara_extractBCORF.R <(samtools view "$ALNFILE") "$REFFASTANOBC" "$EXTRACTDIR" \
  --bcLen $BLEN --bcPos $BCPOS --barseqMode 1

#STEP3: Align barcodes to library
bowtie2 --no-head --norc --very-sensitive --all \
  -x "${LIBRARY%.csv}" \
  -U "${EXTRACTDIR}/bcExtract_1.fastq.gz" 2>/dev/null\
  |gzip -c> "${CLUSTERDIR}/bcMatches.sam.gz"