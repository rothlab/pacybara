
###########################
# PARALLELIZED PRIMARY CLUSTERING
###########################

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