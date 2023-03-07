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

#"${WORKSPACE}counts/" "$SAMPLES" "$LIBRARY"

library(yogitools)
library(argparser)

p <- arg_parser(
  "consolidate bartender counts",
  name="barseq_consolidator.R"
)
p <- add_argument(p, "indir", help="input directory")
p <- add_argument(p, "samples", help="sample table file")
p <- add_argument(p, "library", help="barcode library file")
pargs <- parse_args(p)

#read library
bcLib <- read.csv(pargs$library)
if (ncol(bcLib) > 9) {
  cat("Trimming excess columns from library file!\n")
  bcLib <- bcLib[,1:9]
}
bcs <- bcLib$barcode

#read sample table
samples <- read.delim(pargs$samples)


missing <- list()
countMat <- do.call(cbind,lapply(samples$sample, function(sname) {
  cat("Processing ",sname,"...\n")
  cfile <- list.files(pargs$indir,pattern=sprintf("%s.*cluster.csv",sname),full.names=TRUE)
  if (length(cfile) != 1) {
    stop("No unambiguous file match for ",sname)
  }
  counts <- read.csv(cfile)
  counts <- setNames(counts$time_point_1,counts$Center)
  #Record non-overlapping clusters
  missBC <- setdiff(names(counts),bcs)
  missing[[sname]] <<- counts[missBC]
  #Record the matches
  overlap <- intersect(names(counts),bcs)
  out <- setNames(numeric(length(bcs)),bcs)
  out[overlap] <- counts[overlap]
  return(out)
}))
colnames(countMat) <- samples$sample
out <- cbind(bcLib,countMat)

missingBCs <- Reduce(union,lapply(missing,names))
missingCounts <- do.call(cbind,lapply(missing,function(mcounts) {
  out <- setNames(numeric(length(missingBCs)),missingBCs)
  out[names(mcounts)] <- mcounts
  return(out)
}))

write.csv(out,paste0(pargs$indir,"/allCounts.csv"),row.names=FALSE)
write.csv(missingCounts,paste0(pargs$indir,"/noMatchCounts.csv"))

