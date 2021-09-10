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

library(argparser)

p <- arg_parser(
  "calculate barseq enrichments",
  name="barseq_enrichment.R"
)
p <- add_argument(p, "counts", help="counts file")
p <- add_argument(p, "samples", help="sample table file")
p <- add_argument(p, "outdir", help="output directory")
pargs <- parse_args(p)

counts <- read.csv(pargs$counts)
samples <- read.delim(pargs$samples)

#substitute dashes in sample names (to account for read.csv)
samples$sample <- gsub("-",".",samples$sample)

#' join multiple datapoints weighted by stdev
#' @param ms data means
#' @param sds data stdevs
#' @param dfs degrees of freedom for each data point
#' @param ws datapoint weights. By default these get auto-calculated from sds
#' @return a vector containing the joint mean and joint stdev
join.datapoints <- function(ms,sds,dfs,ws=(1/sds)/sum(1/sds)) {
  #weighted mean
  mj <- sum(ws*ms)
  #weighted joint variance
  vj <- sum(ws*(sds^2+ms^2)) -mj^2
  #return output
  return(c(mj=mj,sj=sqrt(vj),dfj=sum(dfs)))
}

counts$snvs <- NULL
freqs <- apply(counts[,-(1:7)],2,function(col) col/sum(col))

conds <- with(samples,tapply(sample,with(samples,paste(assay,condition,sep=".")),c))
msd <- do.call(cbind,lapply(conds,function(smpl) {
  data.frame(
    mean=apply(freqs[,smpl],1,mean),
    sd=apply(freqs[,smpl],1,sd)
  )
}))

minFreq <- unique(sort(unlist(msd[,grep("mean$",colnames(msd))])))[[2]]
pseudoCount <- 10^(floor(log10(minFreq)))

assays <- unique(samples$assay)

lrs <- do.call(cbind,setNames(lapply(assays, function(assay) {
  sConds <- setdiff(unique(samples[samples$assay==assay,"condition"]),"All")
  do.call(cbind,setNames(lapply(sConds, function(sCond) {
    smeancol <- sprintf("%s.%s.mean",assay,sCond)
    nsmeancol <- sprintf("%s.All.mean",assay)
    ssdcol <- sprintf("%s.%s.sd",assay,sCond)
    nssdcol <- sprintf("%s.All.sd",assay)
    lr <- log10(msd[,smeancol]+pseudoCount)-log10(msd[,nsmeancol])
    logs.sd <- abs(msd[,ssdcol]/(log(10)*(msd[,smeancol]+pseudoCount)))
    logns.sd <- abs(msd[,nssdcol]/(log(10)*msd[,nsmeancol]))
    lr.sd <- sqrt(logs.sd^2 + logns.sd^2)
    data.frame(lr=lr,sd=lr.sd,allfreq=msd[,nsmeancol])
  }),sConds))
}),assays))

lrs <- cbind(counts[,1:7],lrs)

write.csv(lrs,paste0(pargs$outdir,"allLRs.csv"),row.names=FALSE)

isSingleMut <- grepl("^p\\.[A-Za-z]{3}\\d+[A-Za-z]{3}$|^p\\.[A-Za-z]{3}\\d+=$",lrs$aaChangeHGVS)

nreps <- length(unique(samples$replicate))

#Export to MaveDB format
simple.aa.lr <- do.call(c,lapply(assays, function(assay) {
  sConds <- setdiff(unique(samples[samples$assay==assay,"condition"]),"All")
  setNames(lapply(sConds, function(sCond) {
    lrcol <- sprintf("%s.%s.lr",assay,sCond)
    sdcol <- sprintf("%s.%s.sd",assay,sCond)
    fcol <- sprintf("%s.%s.allfreq",assay,sCond)
    hqLRs <- lrs[
      which(isSingleMut & lrs[,sdcol] < 0.6 & lrs[,fcol] > 1e-6),
    ]
    out <- yogitools::as.df(tapply(1:nrow(hqLRs), hqLRs$aaChangeHGVS, function(idxs) {
      join.datapoints(
        hqLRs[idxs,lrcol],
        hqLRs[idxs,sdcol],
        rep(nreps,length(idxs))
      )
    }))
    out <- cbind(data.frame(hgvs_pro=rownames(out)),out)
    out
  }),sprintf("%s.%s",assay,sConds))
}))

scores <- lapply(simple.aa.lr, function(map) {
  nonmed <- median(map$mj[grep("Ter$",map$hgvs_pro)])
  score <- (map$mj-nonmed)/(-nonmed)
  sd <- map$sj/abs(nonmed)
  se <- sd/sqrt(map$dfj)
  data.frame(hgvs_pro=map$hgvs_pro,score=score,sd=sd,se=se,df=map$dfj)
})

invisible(lapply(names(scores),function(mapname) {
  write.csv(scores[[mapname]],paste0(pargs$outdir,"scores_",mapname,".csv"),row.names=FALSE)
}))

