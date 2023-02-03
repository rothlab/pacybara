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

options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE
)

library(argparser)
library(yogitools)
library(hash)

p <- arg_parser(
  "Soft-filter Pacybara libraries",
  name="pacybara_softfilter.R"
)
p <- add_argument(p, "infile", help="clusters csv.gz file")
p <- add_argument(p, "--minRatio","minimum ccs ratio for cluster dominance",default=0.666)
p <- add_argument(p, "--out", help="output file name",default="softfiltered.csv.gz")
pargs <- parse_args(p)

stopifnot({
  pargs$minRatio >= 0.5
  pargs$minRatio <= 1
})

cat("Reading input...\n")
clusters <- read.csv(pargs$infile)

cat("Indexing...\n")
bcIdx <- tapply(clusters$size,clusters$upBarcode,c)
bcIdx <- hash(names(bcIdx),bcIdx)

cat("Calculating clone dominance...\n")
usable <- do.call(c,sapply(1:nrow(clusters), function(row) {
  bc <- clusters[row,"upBarcode"]
  size <- clusters[row,"size"]
  totSize <- sum(bcIdx[[bc]],na.rm=TRUE)
  if (size/totSize > pargs$minRatio) {
    return(row)
  } else {
    return(NULL)
  }
}))

cat("Filtering...\n")
out <- clusters[usable,]
out <- out[out$size > 1,]

cat("Writing output to file ",pargs$out,"\n")
con <- gzfile(pargs$out,open="w")
write.csv(out,file=con,row.names=FALSE)
close(con)


cat("Plotting result...\n")
rest <- clusters[-usable,]
rest <- rest[rest$size > 1,]

fates <- c(
  clean=sum(out$collision==""),
  dominant=sum(out$collision=="collision"),
  `non-dominant`=nrow(rest)
)
outpdf <- gsub("\\.csv\\.gz$",".pdf",pargs$out)
pdf(outpdf,5,5)
barplot(fates,
  border=NA,ylab="barcodes with >2CCS",
  names.arg=c("clean","has dominant\nclone","no dominant\nclone"),
  col=c("darkolivegreen3","gold","firebrick3")
)
dev.off()|>invisible()

cat("Done!\n")

