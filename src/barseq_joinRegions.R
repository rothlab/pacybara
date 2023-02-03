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
