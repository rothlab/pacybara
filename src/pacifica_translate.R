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
library(hash)

p <- arg_parser(
  "translate library",
  name="pacifica_translate.R"
)
p <- add_argument(p, "clusters", help="clusters csv.gz file")
# p <- add_argument(p, "parameters", help="parameter json file")
p <- add_argument(p, "reference", help="fasta file containing reference sequence")
p <- add_argument(p, "--orfStart", help="ORF start position",default=207L)
p <- add_argument(p, "--orfEnd", help="ORF end position",default=2789L)
pargs <- parse_args(p)

outfile <- sub("\\.csv.gz$","_transl.csv.gz",pargs$clusters)

fcon <- file(pargs$reference,open="r")
refSeq <- yogiseq::readFASTA(fcon)[[1]]$toString()
close(fcon)

orfSeq <- substr(refSeq,pargs$orfStart,pargs$orfEnd)
orfLen <- nchar(orfSeq)

builder <- hgvsParseR::new.hgvs.builder.p(3)
cbuilder <- hgvsParseR::new.hgvs.builder.c()

#separate out mutations outside of ORF
cleanGenos <- strsplit(clusters$geno,";") |> pbmclapply(function(muts) {
  if (all(muts=="=")) {
    return(c(main=muts,outside=NA))
  }

  #figure out which variants are outside of the ORF
  pos <- as.integer(gsub("\\D+","",muts))
  outsideIdx <- which(grepl("^-",muts) | pos > orfLen)

  #fix insertions
  isIns <- grepl("ins",muts)
  if (any(isIns)) {
    muts[isIns] <- sapply(muts[isIns], function(mut) {
      pos <- as.integer(sub("ins.+","",mut))
      rest <- sub(".+ins","",mut)
      sprintf("%d_%dins%s",pos-1,pos,rest)
    })
  }

  if (length(outsideIdx) == 0) {
    return(c(main=paste(muts,collapse=";"),outside=NA))
  } else if (length(outsideIdx) == length(muts)) {
    return(c(main="=",outside=paste(muts,collapse=";")))
  } else {
    return(c(
      main=paste(muts[-outsideIdx],collapse=";"),
      outside=paste(muts[outsideIdx],collapse=";")
    ))
  }
},mc.cores=6) |> yogitools::as.df()
# })


transl <- pbmclapply(
  cleanGenos$main,
  function(mut) {
    #remove leading or trailing semicolons
    # mut <- gsub("^;|;$","",mut)
    if (mut=="=") {
      return(list(hgvsc="c.=",hgvsp="p.=",codonChanges="WT",
        codonHGVS="c.=",aaChanges="WT",aaChangeHGVS="p.="))
    }
    tryCatch({
      if (grepl(";",mut)) {
        hgvsParseR::translateHGVS(paste0("c.[",mut,"]"),orfSeq,builder,cbuilder)
      } else {
        hgvsParseR::translateHGVS(paste0("c.",mut),orfSeq,builder,cbuilder)
      }
    }, error=function(e) {
      list(hgvsc=conditionMessage(e),hgvsp=NA_character_,codonChanges=NA_character_,
          codonHGVS=NA_character_,aaChanges=NA_character_,aaChangeHGVS=NA_character_)
    })
  }, mc.cores=8
) |> yogitools::as.df()

out <- cbind(clusters[,c(1:4,6:7,5)],transl,offTarget=cleanGenos$outside)

outcon <- gzfile(outfile,open="w")
write.csv(out,outcon,row.names=FALSE)
close(outcon)

