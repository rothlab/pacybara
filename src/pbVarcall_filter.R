#!/usr/bin/Rscript
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

# library(yogitools)
# library(hash)
library(argparser)


p <- arg_parser(
  "filter pacbio reads",
  name="pbVarcall_filter.R"
)
p <- add_argument(p, "infile", help="input file")
p <- add_argument(p, "ccsCounts", help="CCS read count file")
p <- add_argument(p, "--minCCS", help="minimum number of CCS reads to pass filter",default=2)
p <- add_argument(p, "--minIndel", help="minimum number of reads to accept indels",default=2)
pargs <- parse_args(p)

# pargs <- list(
#   infile="LPL-R03_160PM_CELL1-Cell1-CCS_minPass5_subassembly_varcalls.txt",
#   ccsCounts="ccs_count_per_barcode.txt",
#   minCCS=2,
#   minIndel=3
# )

ccsCounts <- read.delim(pargs$ccsCounts,header=FALSE,row.names=1)

varcalls <- read.delim(pargs$infile,header=FALSE)
colnames(varcalls) <- c("barcode","calls")
varcalls$ccs <- ccsCounts[varcalls[,1],1]

varcallsFiltered <- varcalls[which(varcalls$ccs >= pargs$minCCS),]

sepCalls  <- strsplit(varcallsFiltered$calls,";")
idx <- which(varcallsFiltered$ccs < pargs$minIndel)
sepCalls[idx] <- lapply(sepCalls[idx],function(muts) {
  muts <- muts[!grepl("del|ins",muts)]
  if (length(muts)==0) {
    return("=")
  } else {
    return(muts)
  }
})
varcallsFiltered$fixedcalls <- sapply(sepCalls,paste,collapse=";")

outfile <- sub("\\.txt$","_filtered.txt",pargs$infile)
out <- varcallsFiltered[,c("barcode","fixedcalls")]
write.table(out,outfile,row.names=FALSE,sep="\t",quote=FALSE)

cat("Done!\n")

