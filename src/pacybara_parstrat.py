#!/usr/bin/env python3

"""
pacybara_parstrat.py
The purpose of this script is to formulate a parallelization strategy for 
the main clustering step. It calculates the "worst case units of interdependence"

"""

import sys
import re
import math
import gzip

#pre-clust and edit distance files
preClustFile = str(sys.argv[1])
edFile = str(sys.argv[2])

#declare/initialize hashmap
#uoi : "unit of interdependence"; rid : "read ID"
rid2uoi = {}
uoi2rids = {}

lastUOI = 0

#helper function to add cluster linkage between RIDS
def addLink(rid1,rid2): 
  global lastUOI
  if rid1 in rid2uoi:
    if rid2 in rid2uoi:
      #both already have clusters
      #we only need to do anything if they're not already clustered together
      if rid2uoi[rid1] != rid2uoi[rid2]:
        #add entries from rid2 into cluster from rid1
        uoi1 = rid2uoi[rid1]
        uoi2 = rid2uoi[rid2]
        #change uoi associations for all members of second cluster to first
        for ridi in uoi2rids[uoi2]:
          rid2uoi[ridi] = uoi1 
          uoi2rids[uoi1].append(ridi)
        #delete old uoi entry for second cluster
        del uoi2rids[uoi2]
    else:
      #rid2 is new
      uoi = rid2uoi[rid1]
      rid2uoi[rid2] = uoi
      uoi2rids[uoi].append(rid2)
  else:
    if rid2 in rid2uoi:
      #rid1 is new
      uoi = rid2uoi[rid2]
      rid2uoi[rid1] = uoi
      uoi2rids[uoi].append(rid1)
    else:
      #both are new
      lastUOI = lastUOI+1
      uoi = str(lastUOI)
      rid2uoi[rid1]=uoi
      rid2uoi[rid2]=uoi
      uoi2rids[uoi]=[rid1,rid2]


#parse preclust file
with gzip.open(preClustFile,"rb") as stream:
  for line in stream:
    line = line.decode().rstrip()
    #write entry elements to variables
    if line.startswith('@'):
      #remove leading @ and metadata (cluster size) 
      ridStr = line[1:].split(" ")[0]
      #split by | delimiter
      rids = ridStr.split("|")
      if (len(rids) > 1): 
        for rid in rids[1:]:
          addLink(rids[0],rid)

#parse edit distance file
with gzip.open(edFile,"rb") as stream:
  #skip first line
  header=stream.readline()
  for line in stream:
    line = line.decode().rstrip()
    fields = line.split(",")
    #pick representative from each cluster for mergers
    rid1 = fields[0].strip('\"').split("|")[0]
    rid2 = fields[1].strip('\"').split("|")[0]
    addLink(rid1,rid2)

#write results to stdout
with sys.stdout as stream:
  for uoi in uoi2rids.keys():
    line = '|'.join(uoi2rids[uoi])
    stream.write(line+"\n")

