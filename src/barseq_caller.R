#!/usr/bin/env Rscript
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

options(stringsAsFactors=FALSE)

library(yogiseq)
library(yogitools)
library(hash)
library(argparser)

p <- arg_parser(
  "match barseq data against a barcode library",
  name="barseq_caller.R"
)
# p <- add_argument(p, "fastq", help="input FASTQ file")
p <- add_argument(p, "library", help="barcode library CSV file")
p <- add_argument(p, "--r1", help="R1 read FASTQ file (required)")
p <- add_argument(p, "--r2", help="R2 read FASTQ file (optional)")
p <- add_argument(p, "--bcLen", help="barcode length",default=25L)
p <- add_argument(p, "--output", help="output file. Defaults to <fastq>_hits.csv.gz")
p <- add_argument(p, "--flanking", help="FASTA file containing the barcode flanking sequences",default="flanking.fasta")
p <- add_argument(p, "--rc", help="use reverse complement of reads",flag=TRUE)
p <- add_argument(p, "--maxErr", help="Maximum allowed number of errors in barcode",default=2L)
args <- parse_args(p)

# args <- list(
#   library="../../libraries/LDLR_R01_pacybara.csv",
#   r1="LDLR_test.fastq",
#   r2=NA,bcLen=25L,output=NA,flanking="flanking.fasta",
#   rc=TRUE, maxErr=1L
# )

#intended barcode length
bcLen <- args$bcLen
#load reference barcodes
refFile <- args$library
bcRef <- read.csv(refFile)
cat(nrow(bcRef)," reference barcodes loaded.\n")
# bcRef <- read.csv("pDEST_pool3_subassembly_bc1_varcalls_transl.csv")
#create matcher object
matcher <- new.bc.matcher(bcRef$barcode,errCutoff=args$maxErr)

#load flanking sequences
flFASTA <- args$flanking
flcon <- file(flFASTA,open="r")
flanking <- readFASTA(flcon)
close(flcon)
flseqs <- sapply(flanking,function(ys)ys$toString())

#input file
# r1File <- "barseq_test.fastq.gz"
r1File <- args$r1
if (!yogitools::canRead(r1File)) {
  stop("Unable to find or read FASTQ file ",r1File)
}
filetype <- system2("file",r1File,stdout=TRUE)
if (!grepl("text|gzip",filetype)) {
  stop("R1 must be fastq or fastq.gz formatted.")
}

r2File <- args$r2
if (!is.na(r2File)) {
  if (!yogitools::canRead(r1File)) {
    stop("Unable to find or read FASTQ file ",r1File)
  }
  filetype2 <- system2("file",r1File,stdout=TRUE)
  if (filetype2 != filetype) {
    stop("R1 and R2 must both be either fastq or fastq.gz files!")
  }
}

outfile <- args$output
if (is.na(outfile)) {
  outfile <- sub("\\..*$","_hits.csv.gz",r1File)
}


#open file connection and create FASTQ parser
if (grepl("gzip",filetype)) {
  cat("Detected GZIP FASTQ file.\n")
  con1 <- gzfile(r1File,open="r")
} else if (grepl("text",filetype)) {
  cat("Detected plain FASTQ file.\n")
  con1 <- file(r1File,open="r")
} else {
  stop("Unrecognized file type: ",filetype)
}
fqp1 <- new.fastq.parser(con1)

if (!is.na(r2File)) {
  if (grepl("gzip",filetype)) {
    con2 <- gzfile(r2File,open="r")
  } else {#if this weren't txt it would have failed above
    con2 <- file(r2File,open="r")
  }
  fqp2 <- new.fastq.parser(con2)
}

#open output connection
outcon <- gzfile(outfile,open="w")
cat("hits,diffs,nhits\n",file=outcon)

cat("Processing FASTQ file.\n")
progress <- 0

# extractBCs <- function(rseqs,flseqs,bcLen=25) {
  # starts <- as.vector(regexpr(flseqs[[1]],rseqs)+nchar(flseqs[[1]]))
  # ends <- as.vector(regexpr(flseqs[[2]],rseqs))-1
#   as.vector(mapply(function(str,start,end) {
#     if (start < 0 || end < 0 || end-start+1 != bcLen) NA_character_ else substr(str,start,end)
#   },rseqs,starts,ends))
# }

#extract barcode segments from read sequences
extractBCs <- function(rseqs,flseqs,bcLen=25) {
  #iterate over reads
  sapply(rseqs, function(rseq) {

    #try simple regex match first
    start <- regexpr(flseqs[[1]],rseq)+nchar(flseqs[[1]])
    end <- regexpr(flseqs[[2]],rseq)-1
    if (attr(start,"match.length") > 0 && attr(end,"match.length") > 0 && 
        start >= 0 && end >= 0 && end-start+1 > bcLen-4) {
      return(substr(rseq,start,end))
    }

    #if that didn't work and no Ns exist in the sequence, there's no match
    if (!grepl("N",rseq)) {
      return(NA_character_)
    }

    #othwerwise we'll move to suffix trees with 'N' as wildcards.
    startCandidates <- integer()
    endCandidates <- integer()
    tryCatch({
      #build suffix tree over read sequence
      stree <- yogiseq::suffixTree(rseq)
      #search suffix tree for flanking sequences while allowing N as wildcard
      startCandidates <- yogiseq::searchSuffixTree(stree,flseqs[[1]],wildcard="N")
      endCandidates <- yogiseq::searchSuffixTree(stree,flseqs[[2]],wildcard="N")
    },error=function(e) {
      cat("WARNING: Unable to process read:",rseq,"\n")
    })


    # cat("starts:",startCandidates,"ends:",endCandidates,"\n")

    if (length(startCandidates) < 1 || length(endCandidates) < 1) {
      #no match
      return(NA_character_)
    }

    #find start and end candidates at correct distance
    dists <- do.call(rbind,lapply(startCandidates, function(sc) {
      endCandidates - sc - nchar(flseqs[[1]])
    }))

    usable <- (dists > bcLen-4 & dists < bcLen+4)
    if (sum(usable) != 1) {
      #ambiguous matches
      return(NA_character_)
    }

    idx <- which(usable,arr.ind=TRUE)
    start <- startCandidates[[idx[1,"row"]]]+nchar(flseqs[[1]])
    end <- endCandidates[[idx[1,"col"]]]-1

    return(substr(rseq,start,end))

  })
}

exceptions <- yogitools::new.counter()

#iterate over 100-read chunks
while (length(reads <- fqp1$parse.next(100,ignore.quality=TRUE)) > 0) {

  if (args$rc) {
    #take reverse complement of reads
    reads <- lapply(reads,reverseComplement)
  }
  #and convert to simple strings
  rseqs <- sapply(reads,function(ys)ys$toString())
  #extract barcodes
  bcReads <- extractBCs(rseqs,flseqs,bcLen=bcLen)


  #if there's a R2 read, validate and extract as well
  if (exists("fqp2")) {
    #parse the equivalent 100 R2 reads as well
    reads2 <- fqp2$parse.next(100,ignore.quality=TRUE)
    names1 <- sapply(reads,function(x)x$getID())
    names2 <- sapply(reads2,function(x)x$getID())
    if (!identical(names1,names2)) {
      stop("R1 and R2 files do not contain the same clusters!")
    }

    if (!args$rc) {
      #take reverse complement of reads if appropriate
      reads2 <- lapply(reads2,reverseComplement)
    }
    rseqs2 <- sapply(reads2,function(ys)ys$toString())

    #extract barcodes
    bcReads2 <- extractBCs(rseqs2,flseqs)
    naReads <- is.na(bcReads) & is.na(bcReads2)
  } else {
    bcReads2 <- NA
    naReads <- is.na(bcReads)
  }


  #find matches
  # system.time({
    matches <- matcher$findMatches(bcReads,bcReads2)
  # })

  #record reasons for rejected reads
  exceptions$add("failedExtraction",sum(naReads))
  exceptions$add("noMatch",sum(unlist(matches[,3])==0 & !naReads))
  exceptions$add("ambiguous",sum(unlist(matches[,3])>1))

  #remove ambiguous matches
  if (any(matches$nhits > 1)) {
    matches[which(matches$nhits > 1),"hits"] <- NA
  }
  #write to output stream
  write.table(matches,outcon,col.names=FALSE,quote=FALSE,row.names=FALSE,sep=",")
  flush(outcon)
  progress <- progress+length(rseqs)
  cat("\r",progress," reads processed.     ")
}
cat("\n")

#close file connections
close(con1)
if (exists("con2")) close(con2)
close(outcon)

#print exception counts to log file
cat(exceptions$export(),"\n")

cat("\nDone!\n")
