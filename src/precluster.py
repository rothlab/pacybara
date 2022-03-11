#!/usr/bin/env python3

import sys
import re
import math

#declare/initialize hashmaps
quals = {}
rids = {}

#function to write a new record into the maps
def addRecord(rid,seq,qual):
  #register read ID
  if seq in rids:
    rids[seq].append(rid)
  else:
    rids[seq]=[rid]
  #combine qualities
  qualNums = [ord(x)-33 for x in qual]
  if seq in quals:
    quals[seq] = [min(93,qualNums[i] + quals[seq][i]) for i in range(len(qualNums))]
  else:
    quals[seq] = qualNums

#declare/initialize variables
i=0
rid=seq=qual=""

#stream FASTQ from stdin
# print("Indexing...")
with sys.stdin as stream:
  for line in stream:
    line = line.rstrip()
    #write entry elements to variables
    if line.startswith('@'):
      i=0
      matchObj = re.match(r'@([^ ]+)', line, 0)
      if matchObj:
        rid=matchObj.group(1)
      else:
        raise Exception("Unable to extract read id!")
    i+=1
    if i==2:
      seq=line
    elif i==4:
      qual=line
      addRecord(rid,seq,qual)

#write fastq output
# print("Writing results to file...")
# with open("precluster.fastq","w") as fqFile:
with sys.stdout as fqFile:
  for seq in rids.keys():
    title = '|'.join(rids[seq])
    nreads = len(rids[seq])
    qualStr = ''.join([chr(q+33) for q in quals[seq]])
    fqFile.write(f"@{title} size={nreads}\n")
    fqFile.write(seq+"\n")
    fqFile.write("+"+"\n")
    fqFile.write(qualStr+"\n")

# print("Done!")