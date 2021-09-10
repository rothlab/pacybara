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

library(yogitools)
library(argparser)

p <- arg_parser(
  "consolidate barseq counts",
  name="barseq_consolidator.R"
)
p <- add_argument(p, "indir", help="input directory")
p <- add_argument(p, "samples", help="sample table file")
p <- add_argument(p, "library", help="barcode library file")
pargs <- parse_args(p)


countFiles <- list.files(pargs$indir,pattern="_counts.txt")

samples <- read.delim(pargs$samples)

#re-order files based on sample sheet
countFiles <- sapply(samples$sample, function(sid) {
  hits <- grep(sid,countFiles)
  if (length(hits)==1) {
    return(paste0(pargs$indir,"/",countFiles[[hits]]))
  } else {
    stop("No unambiguous file match for ",sid)
  }
})

# samples <- as.data.frame(
#   yogitools::extract.groups(countFiles,"([^-]+)-([^-]+)-Rep(\\d+)")
# )
# colnames(samples) <- c("Experiment","Condition","Replicate")

# samplenames <- with(samples,{
#   sprintf("%s.%s.rep%s",Experiment,Condition,Replicate)
# })

# bcFile <- "pDEST_pool3_subassembly_bc1_varcalls_transl.csv"
barcodes <- read.csv(pargs$library)

countData <- do.call(data.frame,lapply(countFiles, function(countsFile) {
  counts <- yogitools::as.df(strsplit(trimws(readLines(countsFile))," +"))
  counts <- counts[counts[,2]!="hits",]
  counts[,1] <- as.integer(counts[,1])
  counts[,2] <- as.integer(counts[,2])

  usefulCounts <- na.omit(counts)
  countVec <- integer(nrow(barcodes))
  countVec[usefulCounts[,2]] <- usefulCounts[,1]
  return(countVec)
}))
colnames(countData) <- samples$sample

noMatch <- yogitools::as.df(lapply(countFiles, function(countsFile) {
  counts <- yogitools::as.df(strsplit(trimws(readLines(countsFile))," +"))
  nom <- as.integer(counts[counts[,2]=="NA",1])
  totReads <- as.integer(sum(as.integer(counts[,1],na.rm=TRUE)))
  list(
    noMatch=nom,
    totalReads=totReads,
    frac=nom/totReads
  )
}))
rownames(noMatch) <- samples$sample
write.csv(noMatch,paste0(pargs$indir,"/noMatch.csv"))

out <- cbind(barcodes,countData)
write.csv(out,paste0(pargs$indir,"/allCounts.csv"),row.names=FALSE)
