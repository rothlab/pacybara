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

# BCPOS=153
ORFSTART=207
ORFEND=2789
THREADS=4
MINJACCARD=0.2
MINMATCHES=1
MAXDIFF=2
MINQUAL=100
VIRTUALBC=0
QARG=""
BLACKLIST=""


#helper function to print usage information
usage () {
  cat << EOF

pacybara.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

Runs Pacybara
Usage: pacybara.sh [-b|--barcode <BARCODE>] [-s|--orfStart <ORFSTART>] 
   [-e|--orfEnd <ORFEND>] [--minQual <MINQUAL>] 
   [-m|--minMatches <MINMATCHES>] [--maxDiff <MAXDIFF>] 
   [-j|--minJaccard <MINJACCARD>] [-v|--virtualBC]
   [-c|--cpus <NUMBER>] [-q|--queue <QUEUE>] 
   [--blacklist {<NODE>,}]
   <INFASTQ> <FASTA> [<WORKSPACE>]

-b|--barcode   : The barcode degeneracy code sequence, defaults to
                 $BARCODE
-s|--orfStart  : The ORF start position, defaults to $ORFSTART
-e|--orfEnd    : The ORF end position, defaults to $ORFEND
--minQual      : The minimum PHRED quality for variant basecall to be
                 considered real. Defaults to $MINQUAL
-m|--minMatches: The minimum number of variant matches for a merge
                 to occur. Defaults to $MINMATCHES
--maxDiff      : The maxium allowed edit distance between two clusters
                 for a merge to occur. Defaults to $MAXDIFF
-j|--minJaccard: The minimum Jaccard coefficient between to clusters
                 for a merge to occur. Defaults to $MINJACCARD
-v|--virtualBC : Use virtual barcodes (fusion of up- and down-tags) 
                 for clustering. Otherwise only use uptags.
-c|--cpus      : Number of CPUs to use, defaults to $THREADS
-q|--queue     : The queue (slurm partition) to use
--blacklist    : A comma-separated list of compute nodes not to be used
<INFASTQ>      : The input fastq.gz file to process
<FASTA>        : The raw reference fasta file 
<WORKSPACE>    : The workspace directory. Defaults to 'workspace/'

EOF
 exit $1
}


#Parse Arguments
NUMRX='^[0-9]+$'
FLOATRX='^0\.[0-9]+$'
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
    --minQual)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        if ! [[ $2 =~ $NUMRX ]] ; then
           echo "ERROR: minQual must be a positive integer number" >&2
           usage 1
        fi
        MINQUAL=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -m|--minMatches)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        if ! [[ $2 =~ $NUMRX ]] ; then
           echo "ERROR: minMatches must be a positive integer number" >&2
           usage 1
        fi
        MINMATCHES=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    --maxDiff)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        if ! [[ $2 =~ $NUMRX ]] ; then
           echo "ERROR: maxDiff must be a positive integer number" >&2
           usage 1
        fi
        MAXDIFF=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -j|--minJaccard)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        if ! [[ $2 =~ $FLOATRX ]] ; then
           echo "ERROR: minJaccard must be between 0 and 1" >&2
           usage 1
        fi
        MINJACCARD=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -v|--virtualBC)
      VIRTUALBC=1
      shift
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
    -b|--blacklist)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        BLACKLIST=$2
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
INFASTQ=$1
if ! [[ -r "$INFASTQ" ]]; then
  echo "ERROR: Unable to read FASTQ file $INFASTQ">&2
  exit 1
fi
RX='\.fastq\.gz$'
if ! [[ "$INFASTQ" =~ $RX ]]; then
   echo 'ERROR: First parameter must be a *.fastq.gz file' >&2
   exit 1
fi
# INDIR=$1
# if ! [[ -d "$INDIR" ]]; then
#   echo "ERROR: $INDIR is not a directory!">&2
#   exit 1
# fi
# #add trailing "/" to directory name, if missing
# RX='/$'
# if ! [[ "$INDIR" =~ $RX ]]; then
#   INDIR="${INDIR}/"
# fi

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

if [[ -z $BLACKLIST ]]; then
  BLARG=""
else
  BLARG="--blacklist $BLACKLIST"
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
  #count number of barcodes
  NUMBC=$(echo $FBODY | grep -o $BARCODE | wc -l)
  #cut off the first and last barcode match and anything following it
  PREFIX1=${FBODY%%$BARCODE*}
  PREFIX2=${FBODY%$BARCODE*}
  if [[ $NUMBC == 2 ]]; then
    #the lengths of the remaining prefixes (+1) are the barcode positions
    echo "$((${#PREFIX1}+1)),$((${#PREFIX2}+1))"
  elif [[ $NUMBC == 1 ]]; then
    echo "$((${#PREFIX1}+1))"
  else
    echo "Error: Reference must contain either one or two barcodes!">&2
    exit 1
  fi
}

#find barcode position in reference
BCPOS=$(findBarcodPos $REFFASTA $BARCODE)
#remove barcode from reference
REFFASTANOBC=$(removeBarcode $REFFASTA $BARCODE)
#and create an index for the barcodeless reference
samtools faidx "$REFFASTANOBC"
bwa index -a is "$REFFASTANOBC"


#define intermediate files
OUTPREFIX=$(basename $INFASTQ|sed -r "s/.fastq.gz$//")
#extracted barcodes file
CHUNKDIR="${WORKSPACE}/${OUTPREFIX}_chunks/"
mkdir -p $CHUNKDIR
# mkdir -p ${CHUNKDIR}/logs/

#SPLIT FASTQ INTO CHUNKS
echo "Splitting FASTQ file into job packages"
zcat "$INFASTQ"|split -a 3 -l 200000 --additional-suffix .fastq - \
  "${CHUNKDIR}/$OUTPREFIX"

CHUNKS=$(ls ${WORKSPACE}/${OUTPREFIX}_chunks/${OUTPREFIX}*.fastq)

#PROCESS CHUNKS IN PARALLEL
startJobs() {
  CHUNKS=$1
  JOBS=""
  for CHUNK in $CHUNKS; do
    TAG=$(basename $CHUNK|sed -r "s/\\.fastq//g")
    #start barseq.R job and capture the job-id number
    echo "Submitting $TAG"
    RETVAL=$(submitjob.sh -n $TAG -l ${CHUNKDIR}/${TAG}.log \
      -e ${CHUNKDIR}/${TAG}.log -t 24:00:00 \
      -c 4 -m 4G $BLARG \
      pacybara_worker.sh --barcode $BARCODE --barcodePos "$BCPOS" \
      --orfStart $ORFSTART --orfEnd $ORFEND -c 4 \
      "$CHUNK" "$REFFASTANOBC" "$CHUNKDIR")
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

checkForFailedJobs() {
  CHUNKS=$1
  FAILEDCHUNKS=""
  for CHUNK in $CHUNKS; do
    TAG=$(basename $CHUNK|sed -r "s/\\.fastq//g")
    LOG=${CHUNKDIR}/${TAG}.log
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

echo "Running jobs..."
startJobs "$CHUNKS"
echo "Validating jobs..."
FAILEDCHUNKS=$(checkForFailedJobs "$CHUNKS")
#If any jobs failed, try 2 more times
TRIES=1
while ! [[ -z $FAILEDCHUNKS ]]; do
  if [[ $TRIES < 3 ]]; then
    echo "Attempting to re-run failed jobs."
    startJobs "$FAILEDCHUNKS"
    echo "Validating jobs..."
    FAILEDCHUNKS=$(checkForFailedJobs "$CHUNKS")
    ((TRIES++))
  else
    echo "ERROR: Exhausted 3 attempts at re-running failed jobs!">&2
    exit 1
  fi
done

#DIRECTORY FOR STORING EXTRACTION RESULTS
EXTRACTDIR="${WORKSPACE}/${OUTPREFIX}_extract/"
mkdir -p $EXTRACTDIR

#CONSOLIDATE RESULT CHUNKS
echo "Consolidating alignments and genotypes"
for CHUNK in $CHUNKS; do
  #this for loop is just to absolutely make sure the order is preserved
  TAG=$(basename $CHUNK|sed -r "s/\\.fastq//g")
  cat ${CHUNKDIR}${TAG}_extract/bcExtract_1.fastq.gz>>${EXTRACTDIR}/bcExtract_1.fastq.gz
  if [[ -e ${CHUNKDIR}${TAG}_extract/bcExtract_2.fastq.gz ]]; then
    cat ${CHUNKDIR}${TAG}_extract/bcExtract_2.fastq.gz>>${EXTRACTDIR}/bcExtract_2.fastq.gz
  fi
  cat ${CHUNKDIR}${TAG}_extract/bcExtract_combo.fastq.gz>>${EXTRACTDIR}/bcExtract_combo.fastq.gz
  cat ${CHUNKDIR}${TAG}_extract/genoExtract.csv.gz>>${EXTRACTDIR}/genoExtract.csv.gz
done

samtools cat -o "${OUTPREFIX}_align.bam" ${CHUNKDIR}/*bam
tar czf "${OUTPREFIX}_alignLog.tgz" ${CHUNKDIR}/*log

#consolidate exception counts
Rscript -e '
options(stringsAsFactors=FALSE)
lines<-readLines(pipe(paste("cat",paste(commandArgs(TRUE),collapse=" "))))
mat <- do.call(rbind,lapply(strsplit(lines,","),function(elems) {
  s<-do.call(rbind,strsplit(elems,"="))
  setNames(as.integer(trimws(s[,2])),s[,1])
}))
tot <- colSums(mat)
cat(paste(sprintf("%s=%d",names(tot),tot),collapse="\n"))
cat("\n")
# print(colSums(mat))
' ${CHUNKDIR}*/*exceptions.txt>${EXTRACTDIR}/exceptions.txt

#Run quick QC of barcode length distributions
mkdir ${EXTRACTDIR}/qc
zcat "${EXTRACTDIR}/bcExtract_1.fastq.gz"|grep len=|cut -f 3,3 -d'='|\
  sort -n|uniq -c>"${EXTRACTDIR}/qc/bc1len_distr.txt"
zcat "${EXTRACTDIR}/bcExtract_combo.fastq.gz"|grep len=|cut -f 3,3 -d'='|\
  sort -n|uniq -c>"${EXTRACTDIR}/qc/bccombolen_distr.txt"

#CLEAN UP CHUNKS
rm -r $CHUNKDIR

############
#Clustering
############


#DIRECTORY FOR STORING EXTRACTION RESULTS
CLUSTERDIR="${WORKSPACE}/${OUTPREFIX}_clustering/"
mkdir -p $CLUSTERDIR/qc

#pre-clustering (group fully identical barcode reads)
if [[ $VIRTUALBC == 1 ]]; then
  echo "Indexing virtual barcodes..."
  zcat "${EXTRACTDIR}/bcExtract_combo.fastq.gz"|pacybara_precluster.py\
    |gzip -c>"${CLUSTERDIR}/bcPreclust.fastq.gz"
else
  echo "Indexing uptag barcodes..."
  zcat "${EXTRACTDIR}/bcExtract_1.fastq.gz"|pacybara_precluster.py\
    |gzip -c>"${CLUSTERDIR}/bcPreclust.fastq.gz"
fi
#record distribution of pre-cluster sizes
zcat "${CLUSTERDIR}/bcPreclust.fastq.gz"|grep size|cut -f2,2 -d=|\
  sort -n|uniq -c>"${CLUSTERDIR}/qc/bcPreclust_distr.txt"

#build bowtie index
seqret -sequence <(zcat "${CLUSTERDIR}/bcPreclust.fastq.gz")\
  -outseq "${CLUSTERDIR}/bcPreclust.fasta"
mkdir "${CLUSTERDIR}/db"
bowtie2-build "${CLUSTERDIR}/bcPreclust.fasta" \
  "${CLUSTERDIR}/db/bcDB"
rm "${CLUSTERDIR}/bcPreclust.fasta"

#align all vs all
echo "Aligning all vs all barcodes..."
bowtie2 --no-head --norc --very-sensitive --all \
  -x "${CLUSTERDIR}/db/bcDB" \
  -U "${CLUSTERDIR}/bcPreclust.fastq.gz" 2>/dev/null\
  |awk '{if($1!=$3){print}}'|gzip -c> "${CLUSTERDIR}/bcMatches.sam.gz"
#delete bowtie library files; no longer needed
rm -r "${CLUSTERDIR}/db"

#calculate edit distance
echo "Calculating barcode edit distances..."
pacybara_calcEdits.R "${CLUSTERDIR}/bcMatches.sam.gz" \
  "${CLUSTERDIR}/bcPreclust.fastq.gz" --maxErr "$MAXDIFF" \
  --output "${CLUSTERDIR}/editDistance.csv.gz"

zcat "${CLUSTERDIR}/editDistance.csv.gz"|tail -n +2|cut -f5,5 -d,|sort -n\
  |uniq -c>"${CLUSTERDIR}/qc/edDistr.txt"

# #perform actual clustering and form consensus
# # submitjob.sh -n "${OUTPREFIX}_clustering" -c 12 -m 24G -t 14-00:00:00\
# #   -l "${CLUSTERDIR}/qc/clustering.log" -e "${CLUSTERDIR}/qc/clustering_err.log"\
# #   -q guru -- \
# pacybara_runClustering.R "${CLUSTERDIR}/editDistance.csv.gz" \
#   "${EXTRACTDIR}/genoExtract.csv.gz" "${CLUSTERDIR}/bcPreclust.fastq.gz" \
#   --uptagBarcodeFile "${EXTRACTDIR}/bcExtract_1.fastq.gz" \
#   --out "${CLUSTERDIR}/clusters.csv.gz" --minJaccard "$MINJACCARD" \
#   --minMatches "$MINMATCHES" --maxDiff "$MAXDIFF" --minQual "$MINQUAL" 

############################
# PREPARE FOR PRIMARY CLUSTERING
############################

#calculate parallelization strategy
pacybara_parstrat.py "${CLUSTERDIR}/bcPreclust.fastq.gz" \
  "${CLUSTERDIR}/editDistance.csv.gz"|gzip -c >"${CLUSTERDIR}/parstrat.txt.gz"

#re-instate chunkdir
mkdir -p $CHUNKDIR
#we can't use split -n 200 when reading from stdin
#so we need to calculate the number of lines we need per chunk
CHUNKLINES=$(($(zcat "${CLUSTERDIR}/parstrat.txt.gz"|wc -l) / 200 + 1))
zcat "${CLUSTERDIR}/parstrat.txt.gz"|split -a 3 -l $CHUNKLINES - \
  "${CHUNKDIR}/${OUTPREFIX}_parstrat_"

CHUNKS=$(ls "${CHUNKDIR}${OUTPREFIX}_parstrat_"*)


filterEditDistByRID() {
  CHUNKFILE="$1"
  echo '"query","ref","mismatch","indel","dist"'
python3 -c '
import sys
import gzip
ridsfile=str(sys.argv[1])
ridDict = {}
with open(ridsfile,"r") as stream:
  for line in stream:
    line = line.rstrip()
    rids = line.split("|")
    for rid in rids:
      ridDict[rid] = 1
infile = str(sys.argv[2])
with gzip.open(infile,"rb") as stream:
  for line in stream:
    line = line.decode().rstrip()
    fields = line.split(",")
    if len(fields) < 2:
      break
    rids1 = fields[0].strip("\"").split("|")
    rids2 = fields[1].strip("\"").split("|")
    for rid in (rids1+rids2): 
      if rid in ridDict:
        print(line)
        break
' "$CHUNKFILE" "${CLUSTERDIR}/editDistance.csv.gz"
}

filterGenoByRID() {
  CHUNKFILE="$1"
  OUTFILE=$2
python3 -c '
import sys
import gzip
ridsfile=str(sys.argv[1])
ridDict = {}
with open(ridsfile,"r") as stream:
  for line in stream:
    line = line.rstrip()
    rids = line.split("|")
    for rid in rids:
      ridDict[rid] = 1
infile = str(sys.argv[2])
with gzip.open(infile,"rb") as stream:
  for line in stream:
    line = line.decode().rstrip()
    fields = line.split(",")
    if len(fields) < 1:
      break
    rid = fields[0].strip("\"")
    if rid in ridDict:
      print(line)
' "$CHUNKFILE" "${EXTRACTDIR}/genoExtract.csv.gz"|gzip -c>"$OUTFILE"
}

filterFastqByRID() {
  CHUNKFILE=$1
  FASTQ=$2
  OUTFILE=$3
python3 -c '
import sys
import gzip
import re
ridsfile=str(sys.argv[1])
ridDict = {}
with open(ridsfile,"r") as stream:
  for line in stream:
    line = line.rstrip()
    rids = line.split("|")
    for rid in rids:
      ridDict[rid] = 1
infile = str(sys.argv[2])
with gzip.open(infile,"rb") as stream:
  lcount=5
  for line in stream:
    line = line.decode().rstrip()
    lcount+=1
    if lcount <= 4:
      print(line)
    if line.startswith("@"):
      matchObj = re.match(r"@([^ ]+)", line, 0)
      if matchObj:
        rids1=matchObj.group(1).split("|")
        for rid in rids1:
          if rid in ridDict:
            print(line)
            lcount=1
            break
' "$CHUNKFILE" "$FASTQ"|gzip -c>"$OUTFILE"
}


JOBS=""
for CHUNK in $CHUNKS; do
  TAG=${CHUNK##*/}
  # #extract the set of read IDs in this chunk
  # RIDS=$(cat $CHUNK|tr '\n' '|')
  # #remove trailing pipe symbol for valid regex
  # RIDS=${RIDS%'|'}

  #make edit distance subset
  filterEditDistByRID "$CHUNK"|gzip -c>"${CHUNK}_ed.csv.gz"
  #make genoFile subest
  filterGenoByRID "$CHUNK" "${CHUNK}_geno.csv.gz"
  #make preClustfile subset
  filterFastqByRID "$CHUNK" "${CLUSTERDIR}/bcPreclust.fastq.gz" \
    "${CHUNK}_bcPreclust.fastq.gz"
  #make uptagFile subset
  filterFastqByRID "$CHUNK" "${EXTRACTDIR}/bcExtract_1.fastq.gz" \
    "${CHUNK}_upBC.fastq.gz"
  #run chunks
  echo "Submitting job $TAG"
  RETVAL=$(submitjob.sh -n $TAG -l "${CHUNK}.log" \
      -e "${CHUNK}.log" -t 24:00:00 \
      -c 4 -m 4G $BLARG \
      pacybara_runClustering.R "${CHUNK}_ed.csv.gz" \
      "${CHUNK}_geno.csv.gz" "${CHUNK}_bcPreclust.fastq.gz" \
      --uptagBarcodeFile "${CHUNK}_upBC.fastq.gz" \
      --out "${CHUNK}_clusters.csv.gz" --minJaccard "$MINJACCARD" \
      --minMatches "$MINMATCHES" --maxDiff "$MAXDIFF" --minQual "$MINQUAL"
  )
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

#combine results
combineClusterChunks() {
  zcat -1 $(echo "${CHUNKS}"|head -1)_clusters.csv.gz|head -1
  for CHUNK in $CHUNKS; do
    zcat "${CHUNK}_clusters.csv.gz"|tail +2
  done
}
combineClusterChunks|gzip -c>"${CLUSTERDIR}/clusters.csv.gz"

#clean up chunks
rm -r "$CHUNKDIR"

#analyze cluster size distribution
zcat "${CLUSTERDIR}/clusters.csv.gz"|tail -n +2|cut -f 4 -d,|sort -n\
  |uniq -c>"${CLUSTERDIR}/qc/clusters_distr.txt"

#translate variants to amino acid level and separate off-target mutations
pacybara_translate.R "${CLUSTERDIR}/clusters.csv.gz" \
  "$REFFASTA" --orfStart "$ORFSTART" --orfEnd "$ORFEND"

#run QC analysis
pacybara_qc.R "${CLUSTERDIR}/clusters_transl.csv.gz" \
  "${EXTRACTDIR}/" "${CLUSTERDIR}/qc"

echo "Done!"
