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
options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE
)

library(argparser)
library(yogitools)
# library(tileseqMave)
library(hgvsParseR)
library(yogiseq)
library(pbmcapply)

p <- arg_parser(
  "translate library",
  name="pbVarcall_translate.R"
)
p <- add_argument(p, "varcalls", help="variant call file")
# p <- add_argument(p, "parameters", help="parameter json file")
p <- add_argument(p, "cds", help="fasta file containing CDS sequence")
pargs <- parse_args(p)

outfile <- sub("\\.txt$","_transl.csv",pargs$varcalls)

fcon <- file(pargs$cds,open="r")
cdsSeq <- yogiseq::readFASTA(fcon)[[1]]$toString()
close(fcon)

builder <- hgvsParseR::new.hgvs.builder.p(3)
cbuilder <- hgvsParseR::new.hgvs.builder.c()

vc <- read.delim(pargs$varcalls,header=FALSE)
colnames(vc) <- c("barcode","snvs")
out <- yogitools::as.df(pbmclapply(
# out <- yogitools::as.df(lapply(
  vc[,2],
  function(mut) {
    tryCatch({
      if (mut=="=") {
        list(hgvsc="c.=",hgvsp="p.=",codonChanges="WT",
          codonHGVS="c.=",aaChanges="WT",aaChangeHGVS="p.=")
      } else if (grepl(";",mut)) {
        hgvsParseR::translateHGVS(paste0("c.[",mut,"]"),cdsSeq,builder,cbuilder)
      } else {
        hgvsParseR::translateHGVS(paste0("c.",mut),cdsSeq,builder,cbuilder)
      }
    }, error=function(e) {
      list(hgvsc=conditionMessage(e),hgvsp=NA_character_,codonChanges=NA_character_,
          codonHGVS=NA_character_,aaChanges=NA_character_,aaChangeHGVS=NA_character_)
    })
  }
  ,mc.cores=8
))

out <- cbind(vc,out)

write.csv(out,outfile,row.names=FALSE)
