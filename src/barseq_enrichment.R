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
library(pbmcapply)
library(hash)
library(yogitools)

p <- arg_parser(
  "calculate barseq enrichments",
  name="barseq_enrichment.R"
)
p <- add_argument(p, "counts", help="counts file")
p <- add_argument(p, "samples", help="sample table file")
p <- add_argument(p, "outdir", help="output directory")
p <- add_argument(p, "--wtMed", help="manual override for wildtype median")
p <- add_argument(p, "--nsMed", help="manual override for nonsense median")
p <- add_argument(p, "--freqFilter", help="frequency filter cutoff",default=5e-7)
p <- add_argument(p, "--bnFilter", help="bottleneck filter count cutoff",default=-Inf)
pargs <- parse_args(p)

if (!is.na(pargs$wtMed)) {
  pargs$wtMed <- as.numeric(pargs$wtMed)
  if (is.na(pargs$wtMed)) {
    stop("--wtMed must be numerical")
  }
}
if (!is.na(pargs$nsMed)) {
  pargs$nsMed <- as.numeric(pargs$nsMed)
  if (is.na(pargs$nsMed)) {
    stop("--nsMed must be numerical")
  }
}
if (!is.numeric(pargs$freqFilter) || pargs$freqFilter >= 1 || pargs$freqFilter < 0) {
  stop("--freqFilter must be numeric and must be between 0 and 1!")
}
if (!is.numeric(pargs$bnFilter)) {
  stop("--bnFilter must be numeric!")
}

# pargs <- list(counts="counts/allCounts.csv",samples="~/tmp/upsamples.txt",outdir="~/tmp/",wtMed=NA,nsMed=NA,freqFilter=5e-7,bnFilter=-Inf)
# pargs <- list(counts="counts/allCounts.csv",samples="samples.tsv",outdir="scores/")

counts <- read.csv(pargs$counts)
samples <- read.delim(pargs$samples)

#substitute dashes in sample names (to account for read.csv)
samples$sample <- gsub("-",".",samples$sample)
#adjust to R-convention for numerical column names
if (!all(is.na(as.numeric(samples$sample)))) {
 toFix <- which(!is.na(as.numeric(samples$sample))) 
 samples$sample[toFix] <- paste0("X",samples$sample[toFix])
}
#remove snvs column
# counts$snvs <- NULL

#quick lookup for number of replicates
nreps <- length(unique(samples$replicate))


#' join multiple datapoints weighted by stdev
#' @param ms data means
#' @param sds data stdevs
#' @param dfs degrees of freedom for each data point
#' @param ws datapoint weights. By default these get auto-calculated from sds
#' @return a vector containing the joint mean and joint stdev
join.datapoints <- function(ms,sds,dfs,ws=(1/sds)/sum(1/sds)) {
  l <- length(ms)
  if (l==0) {
    return(c(mj=numeric(0),sj=numeric(0),dfj=numeric(0)))
  } else if (l==1) {
    return(c(mj=ms,sj=sds,dfj=dfs))
  }
  #weighted mean
  mj <- sum(ws*ms)
  #weighted joint variance
  vj <- sum(ws*(sds^2+ms^2)) -mj^2
  #return output
  return(c(mj=mj,sj=sqrt(vj),dfj=sum(dfs)))
}

#calculate frequencies
countsOnly <- counts[,-(1:9)]
freqs <- apply(countsOnly,2,function(col) col/sum(col))

#collapse replicates
conds <- with(samples,tapply(sample,with(samples,paste(assay,condition,sep=".")),c))
msd <- do.call(cbind,pbmclapply(conds,function(smpl) {
  data.frame(
    mean=apply(freqs[,smpl],1,mean),
    sd=apply(freqs[,smpl],1,sd)
  )
},mc.cores=6))
countMeans <- do.call(cbind,pbmclapply(conds,function(smpl){
  data.frame(
    mean=apply(countsOnly[,smpl],1,mean),
    sd=apply(freqs[,smpl],1,sd)
  )
},mc.cores=6))

#determine size of pseudocount
minFreq <- unique(sort(unlist(msd[,grep("mean$",colnames(msd))])))[[2]]
pseudoCount <- 10^(floor(log10(minFreq)))

#calculate LRs
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
lrs <- cbind(counts[,1:9],lrs)

#save to file
write.csv(lrs,paste0(pargs$outdir,"allLRs.csv"),row.names=FALSE)

#normalize to lrs to scores
scores <- do.call(cbind,setNames(lapply(assays, function(assay) {
  sConds <- setdiff(unique(samples[samples$assay==assay,"condition"]),"All")
  wellmeasured <- msd[,paste0(assay,".All.mean")] > pargs$freqFilter
  noBN <- apply(countMeans[,grep(paste0(assay,".*mean"),colnames(countMeans))],1,function(x) all(x >= pargs$bnFilter))
  wtclones <- which(lrs$hgvsc=="c.=" & wellmeasured & noBN)
  stopclones <- which(grepl("Ter$",lrs$hgvsp) & wellmeasured & noBN)
  do.call(cbind,setNames(lapply(sConds, function(sCond) {
    lrcol <- sprintf("%s.%s.lr",assay,sCond)
    sdcol <- sprintf("%s.%s.sd",assay,sCond)
    fcol <- sprintf("%s.%s.allfreq",assay,sCond)
    wtmed <- median(yogitools::fin(lrs[wtclones,lrcol]))
    if (!is.na(pargs$wtMed)) {
      wtmed <- pargs$wtMed
    }
    stopmed <- median(yogitools::fin(lrs[stopclones,lrcol]))
    if (!is.na(pargs$nsMed)) {
      stopmed <- pargs$nsMed
    }
    score <- (lrs[,lrcol]-stopmed)/(wtmed-stopmed)
    score.sd <- lrs[,sdcol]/abs(wtmed-stopmed)
    data.frame(score=score,sd=score.sd,allfreq=lrs[,fcol])
  }),sConds))
}),assays))
scores <- cbind(counts[,1:9],scores)

write.csv(scores,paste0(pargs$outdir,"allScores.csv"),row.names=FALSE)


#single-mutant filter
isSingleMut <- grepl("^p\\.[A-Za-z]{3}\\d+[A-Za-z]{3}$|^p\\.[A-Za-z]{3}\\d+=$",scores$aaChangeHGVS)

#build list of mutation components for multi-mutants
splitMuts <- lapply(
  strsplit(gsub("^p\\.|\\[|\\]","",scores$aaChangeHGVS),";"),
  function(m) paste0("p.",m)
)
#Build index of variant occurrence in clones (~marginal index)
margIdx <- hash::hash()
for (i in 1:length(splitMuts)) {
  if (!any(grepl("Ter$|fs$",splitMuts[[i]]))) {
    for (mut in splitMuts[[i]]) {
      margIdx[[mut]] <- c(margIdx[[mut]],i)
    }
  }
}
#and index of single-mutant clones
singleIdx <- hash::hash()
for (i in which(isSingleMut)) {
  mut <- scores$aaChangeHGVS[[i]]
  singleIdx[[mut]] <- c(singleIdx[[mut]],i)
}

#Calculate variant-centric scores
# allmuts <- keys(margIdx)
allmuts <- unique(do.call(c,splitMuts))
allmuts <- allmuts[grep("^p\\.[A-Za-z]{3}\\d+[A-Za-z]{3}$|^p\\.[A-Za-z]{3}\\d+=$|^p\\.=$",allmuts)]

varscores <- do.call(c,lapply(assays, function(assay) {
  sConds <- setdiff(unique(samples[samples$assay==assay,"condition"]),"All")
  setNames(lapply(sConds, function(sCond) {
    cat("Processing",assay,"-",sCond,"\n")
    scol <- sprintf("%s.%s.score",assay,sCond)
    sdcol <- sprintf("%s.%s.sd",assay,sCond)
    fcol <- sprintf("%s.%s.allfreq",assay,sCond)
    ccol <- sprintf("%s.%s.mean",assay,sCond)

    #for each aa change mutation:
    out <- do.call(rbind,pbmclapply(allmuts, function(mut) {
      
      #first check for single mutants
      if (mut=="p.=") {
        idxs <- which(scores$hgvsp=="p.=")
      } else {
        idxs <- singleIdx[[mut]]
      }
      ffPass <- scores[idxs,fcol] >= pargs$freqFilter
      bnPass <- countMeans[idxs,ccol] >= pargs$bnFilter
      idxs <- idxs[ffPass & bnPass]
      if (length(idxs) > 0) {
        singleOut <- setNames(join.datapoints(
          scores[idxs,scol],
          scores[idxs,sdcol],
          rep(nreps,length(idxs))
        ),c("score.single","score.single.sd","score.single.df"))
      }else {
        singleOut <- c(score.single=NA,score.single.sd=NA,score.single.df=NA)
      }

      #next check for multi-mutant to average over
      midxs <- margIdx[[mut]]
      ffPass <- scores[midxs,fcol] >= pargs$freqFilter
      bnPass <- countMeans[midxs,ccol] >= pargs$bnFilter
      midxs <- midxs[which(ffPass & bnPass)]
      if (length(midxs) > 0) {
        multiOut <- setNames(join.datapoints(
          scores[midxs,scol],
          scores[midxs,sdcol],
          rep(nreps,length(midxs))
        ),c("score.multi","score.multi.sd","score.multi.df"))
      } else {
        multiOut <- c(score.multi=NA,score.multi.sd=NA,score.multi.df=NA)
      }

      #next check if multis can be used to infer singles
      mutLists <- splitMuts[midxs]
      #co-mutation products (cpds)
      cpds <- do.call(rbind,lapply(mutLists, function(ms) {
        if (length(ms) < 2) {
          return(c(NA,NA))
        }
        if (!all(has.key(setdiff(ms,mut),singleIdx))) {
          return(c(NA,NA))
        }
        #find the co-occurring mutations and their respective single-mutant scores
        comuts <- do.call(rbind,lapply(hash::values(singleIdx,setdiff(ms,mut)), function(pidxs) {
          pidxs <- pidxs[which(!is.na(scores[pidxs,sdcol]))]
          join.datapoints(
            scores[pidxs,scol],
            scores[pidxs,sdcol],
            rep(nreps,length(pidxs))
          )
        }))
        if (length(comuts)==0 || nrow(comuts) == 0) {
          return(c(NA,NA))
        }
        #form the product
        coprod <- prod(comuts[,1],na.rm=TRUE)
        coprod.sd <- abs(coprod)*sqrt(sum((comuts[,2]/comuts[,1])^2,na.rm=TRUE))
        if (coprod < 0.1) {
          return(c(NA,NA))
        }
        c(coprod=coprod,coprod.sd=coprod.sd)
      }))
      #only clones with non-NA products are usable
      usable <- which(!is.na(cpds[,1]))
      if (length(usable) > 0) {
        #the inferred score is the multi-mutant score divided by the co-mutation products
        inferred.scores <- scores[midxs[usable],scol]/cpds[usable,1]
        inferred.sds <- abs(inferred.scores)*sapply(usable, function(i) {
          sqrt((scores[midxs[[i]],sdcol]/scores[midxs[[i]],scol])^2 + (cpds[i,2]/cpds[i,1])^2)
        })
        #average over all the inferences
        inferOut <- setNames(join.datapoints(
          inferred.scores,
          inferred.sds,
          rep(nreps,length(inferred.scores))
        ),c("score.infer","score.infer.sd","score.infer.df"))
      } else {
        inferOut <- c(score.infer=NA,score.infer.sd=NA,score.infer.df=NA)
      }

      return(c(singleOut,multiOut,inferOut))
    },mc.cores=6))
    out <- data.frame(hgvs=allmuts,out)
    return(out)
  }),sprintf("%s.%s",assay,sConds))
}))

#write to file
invisible(lapply(names(varscores),function(mapname) {
  write.csv(varscores[[mapname]],paste0(pargs$outdir,"all_aa_scores_",mapname,".csv"),row.names=FALSE)
}))

#generate MaveDB output
jointscores <- setNames(lapply(names(varscores),function(mapname) {
  map <- varscores[[mapname]]
  out <- yogitools::as.df(lapply(1:nrow(map), function(i) with(map[i,],{
    if (!is.na(score.single)) {
      list(hgvs_pro=hgvs,score=score.single,sd=score.single.sd,
        se=score.single.sd/sqrt(score.single.df),
        df=score.single.df
      )
    } else if (!is.na(score.multi)) {
      list(hgvs_pro=hgvs,score=score.multi,sd=score.multi.sd,
        se=score.multi.sd/sqrt(score.multi.df),
        df=score.multi.df
      )
    } else {
      list(hgvs_pro=hgvs,score=NA,sd=NA,se=NA,df=NA)
    }
  })))
  wtrow <- which(out$hgvs_pro=="p.=")
  #remove wt row
  out[-wtrow,]
}),names(varscores))

invisible(lapply(names(jointscores),function(mapname) {
  write.csv(jointscores[[mapname]],paste0(pargs$outdir,"joint_scores_",mapname,".csv"),row.names=FALSE)
}))

