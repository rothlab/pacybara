#!/usr/bin/env Rscript
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
options(stringsAsFactors=FALSE)
library(yogiseq)
library(yogitools)
library(hash)

# new.fastClust <- function() {
#   cl2elems <- list()
#   elem2cl <- hash()
#   lastClust <- 0
#   linkClusters <- function(a,b) {
#     clA <- elem2cl[[a]]
#     clB <- elem2cl[[b]]
#     if (is.null(clA)) {
#       if (is.null(clB)) {
#         #both unknown
#         lastClust <<- lastClust+1
#         clID <- as.character(lastClust)
#         cl2elems[[clID]] <<- c(a,b)
#         elem2cl[[a]] <<- clID
#         elem2cl[[b]] <<- clID
#       } else {
#         #A unknown, add it to cluster B
#         cl2elems[[clB]] <<- c(cl2elems[[clB]],a)
#         elem2cl[[a]] <<- clB
#       }
#     } else {
#       if (is.null(clB)) {
#         #B unknown, add it to cluster A
#         cl2elems[[clA]] <<- c(cl2elems[[clA]],b)
#         elem2cl[[b]] <<- clA
#       } else {
#         #both known
#         if (clA != clB) {
#           #join clusters
#           cl2elems[[clA]] <<- c(cl2elems[[clA]],cl2elems[[clB]])
#           elem2cl[cl2elems[[clB]]] <<- clA
#           cl2elems[[clB]] <<- NULL
#         } 
#         #otherwise, they're already in the same cluster. we're done
#       }
#     }
#   }
#   list(
#     link=linkClusters,
#     getClusters=function() cl2elems,
#     getElements=function(cl) cl2elems[[cl]],
#     getClusterFor=function(el) {
#       if (has.key(el,elem2cl)) {
#         elem2cl[[el]]
#       } else NA
#     },
#     getClSizeFor=function(el) {
#       if (has.key(el,elem2cl)) {
#         length(cl2elems[[elem2cl[[el]]]])
#       } else 1
#     }
#   )
# }


library(argparser)

p <- arg_parser(
  "calculate edit distance between barcode clusters",
  name="calcEdits.R"
)
p <- add_argument(p, "samFile", help="sam.gz file with cluster alignments")
p <- add_argument(p, "preclustFASTQ", help="fastq.gz file of cluster consensus seqs")
p <- add_argument(p, "--output", help="output file",default="editDistance.csv.gz")
p <- add_argument(p, "--maxErr", help="Maximum allowed number of errors in barcode",default=3L)
p <- add_argument(p, "--chunkSize", help="#reads to process at a time",default=1000L)
args <- parse_args(p)

# samFile <- "bcMatches.sam.gz"
samFile <- args$samFile
# fqFile <- "bcPreclust.fastq.gz"
fqFile <- args$preclustFASTQ
# genoFile <- "genoExtract.csv.gz"
# outFile <- "editDistance.csv.gz"
outFile <- args$output
maxDist <- args$maxErr
chunkSize <- args$chunkSize

#get barcode read lengths from fastq file
fqLines <- readLines(fqFile)
headerLines <- fqLines[seq(1,length(fqLines),4)]
seqLines <- fqLines[seq(2,length(fqLines),4)]
bcReadName <- yogitools::extract.groups(headerLines,"^@(\\S+) ")[,1]
bcReadLen <- nchar(seqLines)
# dups <- which(duplicated(bcReadName))
# if (length(dups) > 0) {
#   headerLines <- headerLines[-dups]
#   bcReadName <- bcReadName[-dups]
# }
# bcReadLen <- as.integer(yogitools::extract.groups(headerLines,"len=(\\d+)")[,1])
# names(seqLines) <- bcReadName
rm(fqLines,headerLines,seqLines)
# bcReadLen <- hash(bcReadName,bcReadLen)
names(bcReadLen) <- bcReadName
# bcReadLen[["*"]] <- NA


con <- gzfile(samFile, open="r")
stream <- new.sam.streamer(con,chunkSize)
outcon <- gzfile(outFile,open="w")
cat('"query","ref","mismatch","indel","dist"\n',file=outcon)
# fastclust <- new.fastClust()
# distMatches <- NULL
progress <- 0

tryCatch({
  while (length(nlines <- stream$nextChunk()) > 0) {

    rmat <- with(stream,data.frame(
      getSamElement(,c("cname","rname","cigar")),
      getFlags()[,c("segmentUnmapped","revComp","failQC")],
      getSamElement(,c("pos","AS","NM"))
    ))
    rmat <- rmat[with(rmat,!segmentUnmapped & !revComp & !failQC),]

    #num matches/mismatches
    # nMM <- sapply(
    #   yogitools::global.extract.groups(rmat$cigar,"(\\d+)M"),
    #   function(ms) sum(as.integer(ms))
    # )
    #num insertions
    nIns <- sapply(
      yogitools::global.extract.groups(rmat$cigar,"(\\d+)I"),
      function(ms) sum(as.integer(ms),na.rm=TRUE)
    )
    #num deletions
    nDel <- sapply(
      yogitools::global.extract.groups(rmat$cigar,"(\\d+)D"),
      function(ms) sum(as.integer(ms),na.rm=TRUE)
    )
    #num mismatches
    nMis <- rmat$NM-(nIns+nDel)
    #first, get query and template lengths 
    queryLen <- bcReadLen[rmat$cname]
    tempLen <- bcReadLen[rmat$rname]
    #overhang lengths (positive or negative)
    nClipped <- (queryLen+nDel)-(tempLen+nDel)
    #overhangs are also indels
    nIndel <- nDel+nIns+abs(nClipped)
    totalEditDist <- nMis+nIndel

    if (any(totalEditDist <= maxDist)) {
      idx <- which(totalEditDist <= maxDist)
      out <- cbind(rmat[idx,c("cname","rname")],
        mismatches=nMis[idx],indels=nIndel[idx],distance=totalEditDist[idx]
      )
      write.table(out,outcon,append=TRUE,quote=TRUE,
        sep=",",row.names=FALSE,col.names=FALSE
      )
    }

    progress <- progress + nlines
    cat(progress," processed\n")

  }
},error=function(e) {
  warning(e)
},finally={
  close(outcon)
  # close(con)
})


