#!/usr/bin/env python3

import sys
import re
import math
import gzip

#pre-clust and edit distance files
preClustFile = str(sys.argv[1])
edFile = str(sys.argv[2])
genoFile = str(sys.argv[3])

rid2bc = {}
rid2geno = {}

#load singleton barcodes into hashmap (indexed by readID)
with gzip.open(preClustFile,"rb") as stream:
  lcount=0
  lastRid=""
  for line in stream:
    line=line.decode().rstrip()
    lcount+=1
    #write entry elements to variables
    if line.startswith('@'):
      lcount=1
      #remove leading @ and metadata (cluster size) 
      ridStr=line[1:].split(" ")[0]
      #split by | delimiter
      rids=ridStr.split("|")
      if (len(rids) == 1): 
        lastRid=rids[0]
      else:
        lastRid=""
    elif lcount==2 and lastRid != "":
      rid2bc[lastRid]=line

#remove entries that are found in edit distance file
with gzip.open(edFile,"rb") as stream:
  #skip first line
  header=stream.readline()
  for line in stream:
    line = line.decode().rstrip()
    fields = line.split(",")
    #pick representative from each cluster for mergers
    rids1 = fields[0].strip('\"').split("|")
    rids2 = fields[1].strip('\"').split("|")
    for rid in (rids1+rids2):
      if rid in rid2bc:
        del rid2bc[rid]

#read genotypes into hashmap
with gzip.open(genoFile,"rb") as stream:
  for line in stream:
    line = line.decode().rstrip()
    fields = line.split(",")
    rid = fields[0].strip('\"')
    geno = fields[1].strip('\"')
    if rid in rid2bc:
      rid2geno[rid]=geno

#write output
#"virtualBarcode","upBarcode","reads","size","geno","collision","upTagCollision"
with sys.stdout as stream:
  for rid in rid2bc:
    bc=rid2bc[rid]
    geno=rid2geno[rid]
    stream.write('"'+bc+'",')
    stream.write('"'+bc+'",')
    stream.write('"'+rid+'","1",')
    stream.write('"'+geno+'","",""\n')

