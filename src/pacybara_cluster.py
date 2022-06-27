#!/usr/bin/env python3

"""
pacybara_cluster.py
Performs clustering for pacybara. This is a port of pacybara_runClustering.R
"""

import sys
import re
import math
import gzip
import getopt
import os.path

#define help message
def usage():
  print("""
pacybara_cluster.py v0.0.1

(c) Jochen Weile <jochenweile@gmail.com> 2022

Barcode Clustering for Pacybara

Usage: pacybara_cluster.py [-u|--uptagBarcodeFile <FASTQFILE>]
  [-j|--minJaccard <FLOAT>] [-d|--maxDiff <INT>] [-q|--minQual <INT>]
  [-o|--out <FILE>] [--help] <EDFILE> <GENOFILE> <PRECLUSTFILE>

-u|--uptagBarcodeFile : optional uptag barcode FASTQ file
-j|--minJaccard       : minimum jaccard coefficient to accept cluster
-d|--maxDiff          : maximum edit distance to consider
-q|--minQual          : minimum PHRED quality to accept variant call
-o|--out              : output file. defaults to clusters.csv.gz
--help                : show this help text
EDFILE                : The edit distance file
GENOFILE              : The bcGenos file
PRECLUSTFILE          : The bcPreClust file
""")

#define helper function to throw errors
def usageAndDie(err):
  print(err,file=sys.stderr)
  usage()
  sys.exit(1)

#parse command arguments / options
try: 
  opts, args = getopt.getopt(sys.argv[1:], 
    "u:j:d:q:o:",
    ["uptagBarcodeFile=","minJaccard=","maxDiff=","minQual=","out=","help"]
  )
except getopt.GetoptError as err:
  usageAndDie(err)

#define default options
uptagFile=None
minJaccard=0.3
maxDiff=2
minQual=100
outfile="clusters.csv.gz"

#process the parsed command options
for o, a in opts:
  if o == "--help":
    usage()
    sys.exit()
  elif o in ("-u","--uptagBarcodeFile"):
    uptagFile=a
  elif o in ("-j","--minJaccard"):
    try:
      minJaccard=float(a)
    except:
      usageAndDie("ERROR: minJaccard must be a number")
  elif o in ("-d","--maxDiff"):
    try:
      maxDiff=int(a)
    except:
      usageAndDie("ERROR: maxDiff must be an integer number")
  elif o in ("-q","--minQual"):
    try:
      minQual=int(a)
    except:
      usageAndDie("ERROR: minQual must be an integer number")
  elif o in ("-o","--out"):
    outfile=a

if len(args) < 3:
  usageAndDie("ERROR: Missing required arguments!")
edFile = str(args[0])
genoFile = str(args[1])
preClustFile = str(args[2])

#validate inputs
if not os.path.exists(edFile):
  usageAndDie("ERROR: The file "+edFile+" does not exist")
if not os.path.exists(genoFile):
  usageAndDie("ERROR: The file "+genoFile+" does not exist")
if not os.path.exists(preClustFile):
  usageAndDie("ERROR: The file "+preClustFile+" does not exist")
if not uptagFile is None and not os.path.exists(uptagFile):
  usageAndDie("ERROR: The file "+uptagFile+" does not exist")
if minJaccard > 1 or minJaccard <= 0:
  usageAndDie("ERROR: minJaccard must be between 0 and 1")
if maxDiff < 0 or maxDiff > 3:
  usageAndDie("ERROR: maxDiff must be between 0 and 3")
if minQual < 1:
  usageAndDie("ERROR: minQual must be a positive integer")  


######################
#  FIELDS AND FUNCTION DEFINITIONS
######################

knownRIDPairs = {}

def makePID(rid1,rid2):
  if (rid1 < rid2):
    pid=rid1+"-"+rid2
  else:
    pid=rid2+"-"+rid1
  return pid

def markAsKnown(rid1,rid2):
    knownRIDPairs[makePID(rid1,rid2)]=True

def isKnown(rid1,rid2):
  return (makePID(rid1,rid2) in knownRIDPairs)

#declare/initialize hashmap
#rid : "read ID"
#cid : "cluster ID"
rid2cid = {}
cid2rids = {}

lastCID = 0

#helper function to add cluster linkage between RIDS
def addLink(rid1,rid2): 
  #this means we're accessing the global variable above
  #instead of declaring one in the function scope
  global lastCID
  if rid1 in rid2cid:
    if rid2 in rid2cid:
      #both already have clusters
      #we only need to do anything if they're not already clustered together
      if rid2cid[rid1] != rid2cid[rid2]:
        #add entries from rid2 into cluster from rid1
        cid1 = rid2cid[rid1]
        cid2 = rid2cid[rid2]
        #change cid associations for all members of second cluster to first
        for ridi in cid2rids[cid2]:
          rid2cid[ridi] = cid1 
          cid2rids[cid1].append(ridi)
        #delete old cid entry for second cluster
        del cid2rids[cid2]
    else:
      #rid2 is new
      cid = rid2cid[rid1]
      rid2cid[rid2] = cid
      cid2rids[cid].append(rid2)
  else:
    if rid2 in rid2cid:
      #rid1 is new
      cid = rid2cid[rid2]
      rid2cid[rid1] = cid
      cid2rids[cid].append(rid1)
    else:
      #both are new
      lastCID = lastCID+1
      cid = str(lastCID)
      rid2cid[rid1]=cid
      rid2cid[rid2]=cid
      cid2rids[cid]=[rid1,rid2]

def getClusterForRead(rid):
  return rid2cid.get(rid) #this returns None if it doesn't exist

def getReadsForCluster(cid):
  return cid2rids[cid] #this throws an error if it doesn't exist

#edit distance pairs
#this is a list of dictionaries (for fast searching)
# the list index corresponds to edit distance
# dictionaries link pair IDs to pairs of read IDs
edistPairs = [{} for i in range(maxDiff+1)]
rid2bc = {}
rid2geno = {}
rid2uptag = {}
preClusters = []
preSingles = []


################
#  LOAD DATA
################

print("Loading data")

#process edit distance file
with gzip.open(edFile,"rb") as stream:
  #skip first line
  header=stream.readline()
  for line in stream:
    line = line.decode().rstrip()
    fields = line.split(",")
    rids1=fields[0].strip('\"').split("|")
    rids2 = fields[1].strip('\"').split("|")
    ed = int(fields[4])
    for rid1 in rids1:
      for rid2 in rids2:
        edistPairs[ed][makePID(rid1,rid2)]=(rid1,rid2)

#process preclust file
with gzip.open(preClustFile,"rb") as stream:
  lcount=0
  rids=""
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
      if len(rids)==1:
        preSingles.extend(rids)
      else:
        preClusters.extend([rids])
      #register all pairs of reads in pre-cluster as edit distance 0
      #TODO: This might be a bad idea! could take up tons of memory!
      if len(rids) > 1:
        ridpairs = [(a, b) for i, a in enumerate(rids) for b in rids[i + 1:]]
        for rid1, rid2 in ridpairs:
          edistPairs[0][makePID(rid1,rid2)]=(rid1,rid2)
    #check if the line contains the sequence
    elif lcount==2 and rids!="":
      for rid in rids:
        rid2bc[rid]=line

#process genotypes file
with gzip.open(genoFile,"rb") as stream:
  for line in stream:
    line = line.decode().rstrip()
    fields = line.split(",")
    rid = fields[0].strip('\"')
    geno = fields[1].strip('\"')
    genoDict = {}
    if geno!="=":
      for mut in geno.split(";"):
        m,q=mut.split(":")
        genoDict[m]=int(q)
    rid2geno[rid]=genoDict

#process uptag file
if not uptagFile is None:
  with gzip.open(uptagFile,"rb") as stream:
    lcount=0
    rid=""
    for line in stream:
      line=line.decode().rstrip()
      lcount+=1
      #write entry elements to variables
      if line.startswith('@'):
        lcount=1
        #remove leading @ and metadata (cluster size) 
        rid=line[1:].split(" ")[0]
      #check if the line contains the sequence
      elif lcount==2 and rid!="":
        rid2uptag[rid]=line


###################
# SEED CLUSTERING #
###################

# def varMatrix(rids):
#   joint = {}
#   for rid in rids:
#     for key,val in rid2geno[rid].items():
#       if key in joint:
#         joint[key].extend([val])
#       else:
#         joint[key]=[val]
#   return joint


def varMatrix(rids,penalty=0):
  genos = [rid2geno[rid] for rid in rids]
  muts = []
  for geno in genos:
    muts |= geno.keys()
  varMat=[[geno.get(m,penalty) for m in muts] for geno in genos]
  return list(muts),varMat

# def filterSeqError(varMat):
#   muts=[]
#   for mut,qs in varMat.items():
#     if len(qs) > 1 or sum(qs) > minQual:
#       muts.extend([mut])
#   return { mut : varMat[mut] for mut in muts }

def filterSeqError(muts,varMat):
  idx=[]
  for i in range(len(varMat)):
    observed = len([q for q in varMat[i] if q > 0])
    if observed > 1 or sum(varMat[i]) > minQual:
      idx.append(i)
  return ([muts[i] for i in idx], [varMat[i] for i in idx])


#iterate over pairs at edit distance 0
for rids in preClusters:
  muts,varMat = varMatrix(rids)
  muts,varMat = filterSeqError(muts,varMat)
  

