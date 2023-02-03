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

library(argparser)


p <- arg_parser(
  "Prep Pacybara libraries for use in barseq",
  name="pacybara_preplib.R"
)
p <- add_argument(p, "clusters", help="translated clusters csv.gz file")
p <- add_argument(p, "outfile", help="output file name")
p <- add_argument(p, "--mode", help="translation mode. Valid arguments: up2up, down2down, virt2up, virt2down.
'up2up' directly writes uptag-based clusters to an uptag library.
'down2down' directly writes downtag-based clusters to a downtag library.
'virt2up' converts virtual-barcode-based clusters into an uptag library.
'virt2down' converts virtual-barcode-based clusters into a downtag library.",default="up2up")
pargs <- parse_args(p)

pargs$mode <- match.arg(pargs$mode,c("up2up","down2down","virt2up","virt2down"))

cat("Reading input...\n")
clusters <- read.csv(pargs$clusters)

cat("Converting...\n")
barcodes <- switch(pargs$mode,
  up2up={
    clusters$virtualBarcode
  },
  down2down={
    clusters$virtualBarcode
  },
  virt2up={
    clusters$upBarcode
  },
  virt2down={
    with(clusters,substr(virtualBarcode, nchar(upBarcode)+1,nchar(virtualBarcode)))
  }
)

cat("Writing output...\n")
out <- with(clusters,data.frame(
  barcode=barcodes,
  hgvsc=hgvsc,
  hgvsp=hgvsp,
  codonChanges=codonChanges,
  codonHGVS=codonHGVS,
  aaChanges=aaChanges,
  aaChangeHGVS=aaChangeHGVS,
  offTarget=offTarget,
  size=size
))
write.csv(out,pargs$out,row.names=FALSE)

cat("Done!\n")


