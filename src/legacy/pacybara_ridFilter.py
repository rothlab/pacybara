#!/usr/bin/env python3
import sys
import gzip
import re

#filter edit distance
def filterEditDistance(ridsfile,infile):
  # ridsfile=str(sys.argv[1])
  ridDict = {}
  with open(ridsfile,"r") as stream:
    for line in stream:
      line = line.rstrip()
      rids = line.split("|")
      for rid in rids:
        ridDict[rid] = 1
  # infile = str(sys.argv[2])
  print('"query","ref","mismatch","indel","dist"')
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

#filter genotypes
def filterGenotypes(ridsfile,infile):
  # ridsfile=str(sys.argv[1])
  ridDict = {}
  with open(ridsfile,"r") as stream:
    for line in stream:
      line = line.rstrip()
      rids = line.split("|")
      for rid in rids:
        ridDict[rid] = 1
  # infile = str(sys.argv[2])
  with gzip.open(infile,"rb") as stream:
    for line in stream:
      line = line.decode().rstrip()
      fields = line.split(",")
      if len(fields) < 1:
        break
      rid = fields[0].strip("\"")
      if rid in ridDict:
        print(line)

#filter fastq
def filterFASTQ(ridsfile,infile):
  # ridsfile=str(sys.argv[1])
  ridDict = {}
  with open(ridsfile,"r") as stream:
    for line in stream:
      line = line.rstrip()
      rids = line.split("|")
      for rid in rids:
        ridDict[rid] = 1
  # infile = str(sys.argv[2])
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


ridsfile=str(sys.argv[1])
infile = str(sys.argv[2])
mode = str(sys.argv[3])

if (mode == "ed"):
  filterEditDistance(ridsfile,infile)
elif (mode == "geno"):
  filterGenotypes(ridsfile,infile)
elif (mode == "fastq"):
  filterFASTQ(ridsfile,infile)
else:
  print('invalid mode',file=sys.stderr)
  sys.exit(1)

