#!/usr/bin/env Rscript

#Region combiner

options(stringsAsFactors=FALSE)
library(yogitools)
library(argparser)

p <- arg_parser(
  "draw QC plots for barseq",
  name="barseq_qc.R"
)
p <- add_argument(p, "regionFile", help="regions definition TSV file")
p <- add_argument(p, "--out", help="output file",default="jointMap.csv")
pargs <- parse_args(p)

regionFile <- pargs$regionFile
# regionFile <- commandArgs(TRUE)[[1]]
# regionFile <- "regions.tsv"

regions <- read.delim(regionFile)

jointMap <- do.call(rbind,yogitools::rowApply(regions,
  function(region,startAA,endAA,file) {
    if (!is.na(file)) {
      map <- read.csv(file)
      pos <- as.integer(gsub("\\D+","",map$hgvs_pro))
      inRegion <- which(pos >= startAA & pos <= endAA)
      map[inRegion,]
    } else NULL
  }
))

write.csv(jointMap,pargs$out,row.names=FALSE)
