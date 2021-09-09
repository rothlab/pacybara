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
p <- add_argument(p, "fastq", help="input FASTQ file")
p <- add_argument(p, "library", help="barcode library CSV file")
p <- add_argument(p, "--bcLen", help="barcode length",default=25L)
p <- add_argument(p, "--output", help="output file. Defaults to <fastq>_hits.csv.gz")
p <- add_argument(p, "--flanking", help="FASTA file containing the barcode flanking sequences",default="flanking.fasta")
p <- add_argument(p, "--rc", help="use reverse complement of reads",flag=TRUE)
p <- add_argument(p, "--maxErr", help="Maximum allowed number of errors in barcode",default=2L)
args <- parse_args(p)


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
# fastqFile <- "barseq_test.fastq.gz"
fastqFile <- args$fastq
outfile <- args$output
if (is.na(outfile)) {
	outfile <- sub("\\..*$","_hits.csv.gz",fastqFile)
}
filetype <- system2("file",fastqFile,stdout=TRUE)


#open file connection and create FASTQ parser
if (grepl("gzip",filetype)) {
	cat("Detected GZIP FASTQ file.\n")
	con <- gzfile(fastqFile,open="r")
} else if (grepl("text",filetype)) {
	cat("Detected plain FASTQ file.\n")
	con <- file(fastqFile,open="r")
} else {
	stop("Unrecognized file type: ",filetype)
}
fqp <- new.fastq.parser(con)

#open output connection
outcon <- gzfile(outfile,open="w")
cat("hits,diffs,nhits\n",file=outcon)

cat("Processing FASTQ file.\n")
progress <- 0

while (length(reads <- fqp$parse.next(100,ignore.quality=TRUE)) > 0) {

	if (args$rc) {
		#take reverse complement of reads
		reads <- lapply(reads,reverseComplement)
	}
	#and convert to simple strings
	rseqs <- sapply(reads,function(ys)ys$toString())

	#extract barcodes
	starts <- as.vector(regexpr(flseqs[[1]],rseqs)+nchar(flseqs[[1]]))
	ends <- as.vector(regexpr(flseqs[[2]],rseqs))-1
	bcReads <- as.vector(mapply(function(str,start,end) {
		if (start < 0 || end < 0 || end-start+1 != bcLen) NA else substr(str,start,end)
	},rseqs,starts,ends))

	#find matches
	matches <- matcher$findMatches(bcReads)
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

close(con)
close(outcon)

cat("Done!\n")
