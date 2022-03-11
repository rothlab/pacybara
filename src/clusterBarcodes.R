#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
library(yogiseq)
library(yogitools)
library(hash)

new.fastClust <- function() {
  cl2elems <- list()
  elem2cl <- hash()
  lastClust <- 0
  linkClusters <- function(a,b) {
    clA <- elem2cl[[a]]
    clB <- elem2cl[[b]]
    if (is.null(clA)) {
      if (is.null(clB)) {
        #both unknown
        lastClust <<- lastClust+1
        clID <- as.character(lastClust)
        cl2elems[[clID]] <<- c(a,b)
        elem2cl[[a]] <<- clID
        elem2cl[[b]] <<- clID
      } else {
        #A unknown, add it to cluster B
        cl2elems[[clB]] <<- c(cl2elems[[clB]],a)
        elem2cl[[a]] <<- clB
      }
    } else {
      if (is.null(clB)) {
        #B unknown, add it to cluster A
        cl2elems[[clA]] <<- c(cl2elems[[clA]],b)
        elem2cl[[b]] <<- clA
      } else {
        #both known
        if (clA != clB) {
          #join clusters
          cl2elems[[clA]] <<- c(cl2elems[[clA]],cl2elems[[clB]])
          elem2cl[cl2elems[[clB]]] <<- clA
          cl2elems[[clB]] <<- NULL
        } 
        #otherwise, they're already in the same cluster. we're done
      }
    }
  }
  list(
    link=linkClusters,
    getClusters=function() cl2elems,
    getElements=function(cl) cl2elems[[cl]],
    getClusterFor=function(el) {
      if (has.key(el,elem2cl)) {
        elem2cl[[el]]
      } else NA
    },
    getClSizeFor=function(el) {
      if (has.key(el,elem2cl)) {
        length(cl2elems[[elem2cl[[el]]]])
      } else 1
    }
  )
}


samFile <- "bcMatches.sam.gz"
fqFile <- "bcPreclust.fastq.gz"
genoFile <- "genoExtract.csv.gz"

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
rm(fqLines,headerLines,seqLines)
bcReadLen <- hash(bcReadName,bcReadLen)
bcReadLen[["*"]] <- NA

chunkSize <- 1000L
maxDist <- 3

con <- gzfile(samFile, open="r")
stream <- new.sam.streamer(con,chunkSize)
fastclust <- new.fastClust()
distMatches <- NULL
progress <- 0

while (length(nlines <- stream$nextChunk()) > 0) {

  rmat <- with(stream,data.frame(
    getSamElement(,c("cname","rname","cigar")),
    getFlags()[,c("segmentUnmapped","revComp","failQC")],
    getSamElement(,c("pos","AS","NM"))
  ))
  #num matches/mismatches
  nMM <- sapply(
    yogitools::global.extract.groups(rmat$cigar,"(\\d+)M"),
    function(ms) sum(as.integer(ms))
  )
  nIns <- sapply(
    yogitools::global.extract.groups(rmat$cigar,"(\\d+)I"),
    function(ms) sum(as.integer(ms),na.rm=TRUE)
  )
  nDel <- sapply(
    yogitools::global.extract.groups(rmat$cigar,"(\\d+)D"),
    function(ms) sum(as.integer(ms),na.rm=TRUE)
  )
  #query and template lengths 
  queryLen <- values(bcReadLen,rmat$cname)
  tempLen <- values(bcReadLen,rmat$rname)
  rmat$end <- with(rmat,pos-1+alignLength)
  rmat$diffs <- with(rmat,pos-1+NM+(tempLen-end))

  rmat <- rmat[with(rmat,!segmentUnmapped & !revComp & !failQC),]

  #pre-cluster all identical reads
  invisible(apply(
    rmat[which(rmat$diffs==0),1:2], 1, 
    function(ids) fastclust$link(ids[[1]],ids[[2]])
  ))

  #return more distant matches
  distMatches <- rbind(
    distMatches,
    rmat[with(rmat,diffs > 0 & diffs <= maxDist),c("cname","rname","diffs")]
  )

  progress <- progress + nlines
  cat(progress," processed\n")

}

# close(con)


distMatches$clustC <- sapply(distMatches$cname,function(x)fastclust$getClusterFor(x))
distMatches$clustR <- sapply(distMatches$rname,function(x)fastclust$getClusterFor(x))

singleIdx <- with(distMatches,is.na(clustC) | is.na(clustR))

singletons <- distMatches[singleIdx,]
nonsingles <- distMatches[!singleIdx,]
nonsingles$comboID <- apply(nonsingles[,c("clustC","clustR")],1,function(ids) paste(sort(ids),collapse="-"))
diffsByCluster <- as.df(with(nonsingles,tapply(1:nrow(nonsingles),comboID,function(is){
  i <- is[[1]]
  list(
    clust1=clustC[[i]],
    clust2=clustR[[i]],
    size1=fastclust$getClSizeFor(cname[[i]]),
    size2=fastclust$getClSizeFor(rname[[i]]),
    mdiff=median(diffs[is])
  )
})))

# distMatches$sizeC <- sapply(distMatches$cname,function(x)fastclust$getClSizeFor(x))
# distMatches$sizeR <- sapply(distMatches$rname,function(x)fastclust$getClSizeFor(x))


# fastclust <- new.fastClust()
# fastclust$link("A","B")
# fastclust$link("D","E")
# fastclust$getClusters()

