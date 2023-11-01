
EXTRACTDIR=
CLUSTERDIR=
CHUNKDIR="${WORKSPACE}/${OUTPREFIX}_chunks/"
OUTPREFIX=
MINJACCARD=
MINMATCHES=
MAXDIFF=
MINQUAL=

###########################
# PARALLELIZED PRIMARY CLUSTERING
###########################

if [[ "$CLUSTERMODE" == "downtag" ]]; then
  TAGFILE="${EXTRACTDIR}/bcExtract_2.fastq.gz"
else
  TAGFILE="${EXTRACTDIR}/bcExtract_1.fastq.gz"
fi

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

JOBS=""
for CHUNK in $CHUNKS; do
  TAG=${CHUNK##*/}
  # #extract the set of read IDs in this chunk
  # RIDS=$(cat $CHUNK|tr '\n' '|')
  # #remove trailing pipe symbol for valid regex
  # RIDS=${RIDS%'|'}

  #make edit distance subset
  python3 pacybara_ridFilter.py "$CHUNK" "${CLUSTERDIR}/editDistance.csv.gz" "ed"|gzip -c>"${CHUNK}_ed.csv.gz"
  #make genoFile subest
  python3 pacybara_ridFilter.py "$CHUNK" "${EXTRACTDIR}/genoExtract.csv.gz" "geno"|gzip -c>"${CHUNK}_geno.csv.gz"
  #make preClustfile subset
  python3 pacybara_ridFilter.py "$CHUNK" "${CLUSTERDIR}/bcPreclust.fastq.gz" "fastq"|gzip -c>"${CHUNK}_bcPreclust.fastq.gz"
  #make uptagFile subset
  python3 pacybara_ridFilter.py "$CHUNK" "$TAGFILE" "fastq"|gzip -c>"${CHUNK}_tag.fastq.gz"
  #run chunks
  echo "Submitting job $TAG"

  RETVAL=$(submitjob.sh -n $TAG -l "${CHUNK}.log" \
      -e "${CHUNK}.log" -t 7-00:00:00 \
      -c 2 -m 8G $BLARG \
      pacybara_cluster.py "
      --uptagBarcodeFile "${CHUNK}_tag.fastq.gz" \
      --out "${CHUNK}_clusters.csv.gz" --minJaccard "$MINJACCARD" \
      --minMatches "$MINMATCHES" --maxDiff "$MAXDIFF" --minQual "$MINQUAL" \
      ${CHUNK}_ed.csv.gz" "${CHUNK}_geno.csv.gz" \
      "${CHUNK}_bcPreclust.fastq.gz"
  )
  JOBID=${RETVAL##* }
    if [ -z "$JOBS" ]; then
      JOBS="$JOBID"
    else
      JOBS="${JOBS},${JOBID}"
    fi
done

waitForJobs.sh -v "$JOBS"

#combine results
combineClusterChunks() {
  #print header from first file
  zcat $(echo "${CHUNKS}"|head -1)_clusters.csv.gz|head -1
  #print bodies of all files
  for CHUNK in $CHUNKS; do
    zcat "${CHUNK}_clusters.csv.gz"|tail +2
  done
}
combineClusterChunks|gzip -c>"${CLUSTERDIR}/clusters.csv.gz"

#clean up chunks
rm -r "$CHUNKDIR"