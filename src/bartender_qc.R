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
library(yogitools)

p <- arg_parser(
  "draw QC plots for pacybartender",
  name="bartender_qc.R"
)
p <- add_argument(p, "lrs", help="lrs file")
p <- add_argument(p, "counts", help="counts file")
p <- add_argument(p, "samples", help="sample table txt file")
p <- add_argument(p, "outdir", help="output directory")
p <- add_argument(p, "--logfolder", help="folder containing bartender log files")
p <- add_argument(p, "--deCumufy",help="counteract cumulative counts",flag=TRUE)
p <- add_argument(p, "--freqFilter", help="frequency filter cutoff",default=5e-7)
p <- add_argument(p, "--bnFilter", help="bottleneck filter count cutoff",default=-Inf)
pargs <- parse_args(p)

# pargs <- list(lrs="scores/allLRs.csv",counts="counts/allCounts.csv",samples="samples.txt",outdir="qc/",freqFilter=5e-7,bnFilter=-Inf)

dir.create(pargs$outdir,recursive=TRUE)

lrs <- read.csv(pargs$lrs)
counts <- read.csv(pargs$counts)
samples <- read.delim(pargs$samples)
#substitute dashes in sample names (to account for read.csv)
samples$sample <- gsub("-",".",samples$sample)
#adjust to R-convention for numerical column names
if (!all(is.na(as.numeric(samples$sample)))) {
 toFix <- which(!is.na(as.numeric(samples$sample))) 
 samples$sample[toFix] <- paste0("X",samples$sample[toFix])
}


cmat <- counts[,samples$sample]
# cmat <- counts[,grep("Rep",colnames(counts))]
snames <- with(samples,tapply(sample,paste(assay,condition,sep="."),c))
# snames <- tapply(colnames(cmat),gsub("\\.Rep.*","",colnames(cmat)),c)
cmeans <- do.call(cbind,lapply(snames, function(cols) {
  apply(cmat[,cols],1,mean)
}))

# Plot count distribution ---------------------

fcols <- grep("\\.allfreq$",colnames(lrs))
assays <- sub("\\..+$","",colnames(lrs)[fcols])
fcols <- fcols[!duplicated(assays)]

minFreq <- unique(sort(unlist(lrs[,fcols])))[[2]]
pseudoCount <- 10^(floor(log10(minFreq)))
maxFreq <- max(lrs[,fcols])
breaks <- seq(log10(pseudoCount),ceiling(log10(maxFreq)),0.05)

pdf(paste0(pargs$outdir,"bcFreqDistr.pdf"),8.5,11)
op <- par(mfrow=c(length(fcols),1),oma=c(2,2,2,2))
for (fcol in fcols) {
  hist(log10(lrs[,fcol]),
    breaks=breaks,freq=FALSE,ylim=c(0,2),
    col="steelblue3",border=NA,
    xlab=expression(log[10](barcode~frequency)),
    main=sub("\\..+$","",colnames(lrs)[[fcol]])
  )
  grid(NULL,NULL)
}
par(op)
invisible(dev.off())



synCountIdx <- sapply(strsplit(gsub("p\\.|\\[|\\]","",counts$hgvsp),";"), function(muts) {
  all(grepl("=$",muts))
})
wtCountIdx <- counts$codonChanges=="WT"
stopCountsIdx <- grepl("Ter",counts$hgvsp)
fsCountIdx <- grepl("fs",counts$aaChanges)


conditionCombos <- do.call(c,lapply(unique(samples$assay),function(assay) {
  conds <- unique(samples[samples$assay==assay,"condition"])
  selConds <- setdiff(conds,"All")
  allName <- sprintf("%s.%s",assay,"All")
  lapply(selConds,function(selCond) {
    c(allName,sprintf("%s.%s",assay,selCond))
  })
}))

pdf(paste0(pargs$outdir,"countChange.pdf"),8.5,11)
op <- par(mfrow=c(2,2))
lapply(conditionCombos, function(combo) {
  plot(cmeans[which(synCountIdx | wtCountIdx),combo]+.1,col=yogitools::colAlpha("chartreuse3",.2),pch=20,log="xy")
  points(cmeans[which(stopCountsIdx | fsCountIdx),combo]+.1,col=yogitools::colAlpha("firebrick3",.2),pch=20)
})|>invisible()
par(op)
invisible(dev.off())


# panelfun <- function(x,y,...) {
#   points(x,y,...)
#   grid(NULL,NULL)
#   abline(v=0:1,h=0:1,col=c("firebrick3","chartreuse3"),lty="dashed")
#   corr <- cor(yogitools::fin(cbind(x,y)))[1,2]
#   text(0.5,0.5,sprintf("R=%.02f",corr),col="blue",cex=1.2)
# }
# with(varscores[["Uptake.F5"]],{
#   scm <- score.multi
#   scm[scm==score.single] <- NA
#   pairs(
#     cbind(score.single,scm,score.infer),
#     labels=c("single-mutants","average over\nmulti-mutants","inverse\nmultiplicative\nmodel"),
#     pch=".",col=yogitools::colAlpha(1,0.3),
#     xlim=c(-.5,1.5),ylim=c(-.5,1.5),lower.panel=panelfun,upper.panel=panelfun
#   )
# })

# bioreps <- with(scores[isSingleMut & scores$Uptake.F4.allfreq > 5e-7,],{
#   do.call(rbind,tapply(1:length(hgvsp),hgvsp,function(is) {
#     if (length(is)>1) {
#       t(combn(Uptake.F5.score[is],2))
#     } else NULL
#   }))
# })
# plot(bioreps,
#   xlab="replicate clone A score",ylab="replicate clone B score",
#   pch=".",col=yogitools::colAlpha(1,0.3)
# )
# grid(NULL,NULL)
# abline(v=0:1,h=0:1,col=c("firebrick3","chartreuse3"),lty="dashed")
# corr <- cor(yogitools::fin(bioreps))[1,2]
# text(0.5,0.5,sprintf("R=%.02f",corr),col="blue",cex=1.2)


#plot AllFreq vs CV
sampleNames <- unique(sub("\\.allfreq","",colnames(lrs)[grep("allfreq$",colnames(lrs))]))
pdf(paste0(pargs$outdir,"freqCV.pdf"),8.5,11)
op <- par(mfrow=c(3,2),oma=c(2,2,2,2))
for (sample in sampleNames) {
  fcol <- paste0(sample,".allfreq")
  sdcol <- paste0(sample,".sd")
  lrcol <- paste0(sample,".lr")
  plot(
    lrs[,fcol]+pseudoCount, lrs[,sdcol]/lrs[,lrcol],
    xlab="read frequency 'all'",ylab="CV(LR)",
    log="xy",pch=".",col=yogitools::colAlpha(1,0.2),
    main=sample
  )
}
par(op)
invisible(dev.off())


#find all conditions
conds <- sub("\\.allfreq$","",colnames(lrs)[grep("allfreq",colnames(lrs))])
# relCond <- conds[order(conds)][[length(conds)]]
#iterate over conditions
for (relCond in conds) {

  afCol <- paste0(relCond,".allfreq")
  lrCol <- paste0(relCond,".lr")

  #apply frequency filter
  # hqLRs <- lrs[which(lrs$Uptake.F5.allfreq > 1e-6),]
  ffPass <- lrs[,afCol] > pargs$freqFilter
  # bnPass <- apply(cmeans,1,function(x) all(x >= pargs$bnFilter))
  bnPass <- cmeans[,relCond] >= pargs$bnFilter
  hqLRs <- lrs[which(ffPass & bnPass),]

  #get subgroups of variants (WT, stop, frameshift, synonymous, missense)
  isSyn <- sapply(strsplit(gsub("p\\.|\\[|\\]","",hqLRs$hgvsp),";"), function(muts) {
    all(grepl("=$",muts))
  })
  wtLRs <- hqLRs[hqLRs$codonChanges=="WT",]
  stopLRs <- hqLRs[grepl("Ter",hqLRs$hgvsp),]
  fsLRs <- hqLRs[grepl("fs",hqLRs$aaChanges),]
  synLRs <- hqLRs[isSyn,]
  misLRs <- hqLRs[!(isSyn | hqLRs$codonChanges=="WT" | grepl("Ter",hqLRs$hgvsp) | grepl("fs",hqLRs$aaChanges)),]


  pdf(paste0(pargs$outdir,"barseqDistributions_",relCond,".pdf"),8.5,11)
  op <- par(mfcol=c(3,2),oma=c(10,2,2,2),mar=c(5,4,1,1))
  breaks <- seq(-5,5,0.05)
  hist(
    wtLRs[,lrCol],breaks=breaks,
    col="darkolivegreen3",border=NA,main="WT clones",
    xlab=sprintf("log10(%s/All)",relCond)
  )
  abline(v=mean(yogitools::fin(wtLRs[,lrCol])),lty="dotted")
  hist(
    stopLRs[,lrCol],breaks=breaks,
    col="firebrick3",border=NA,main="Nonsense clones",
    xlab=sprintf("log10(%s/All)",relCond)
  )
  abline(v=mean(yogitools::fin(stopLRs[,lrCol])),lty="dotted")
  hist(
    synLRs[,lrCol],breaks=breaks,
    col="darkolivegreen2",border=NA,main="Synonymous clones",
    xlab=sprintf("log10(%s/All)",relCond)
  )
  abline(v=mean(yogitools::fin(synLRs[,lrCol])),lty="dotted")
  hist(
    fsLRs[,lrCol],breaks=breaks,
    col="firebrick4",border=NA,main="Frameshift clones",
    xlab=sprintf("log10(%s/All)",relCond)
  )
  abline(v=mean(yogitools::fin(fsLRs[,lrCol])),lty="dotted")
  hist(
    misLRs[,lrCol],breaks=breaks,
    col="gray",border=NA,main="Missense clones",
    xlab=sprintf("log10(%s/All)",relCond)
  )
  abline(v=mean(yogitools::fin(misLRs[,lrCol])),lty="dotted")
  par(op)
  invisible(dev.off())

}

### Gather data for read fate plot

logsBySample <- sapply(samples$sample, function(sname) {
  list.files(pargs$logfolder,pattern=sname,full.name=TRUE)
})
extractStats <- yogitools::as.df(lapply(logsBySample, function(lfile) {
  lines <- readLines(lfile, 3L)
  yogitools::extract.groups(lines,"there are (\\d+) ")[,1] |> 
    as.integer()|>setNames(c("total","extract","passfilter"))
}))
extractStats$failedExtraction <- with(extractStats,total-extract)
extractStats$failedFilter <- with(extractStats,extract-passfilter)

nmFile <- sub("allCounts.csv$","noMatchCounts.csv",pargs$counts)
nmTable <- read.csv(nmFile,row.names=1)
nmStats <- t(apply(nmTable,2,function(cs) c(singletons=sum(cs==1),clustered=sum(cs[cs>1]))))
nmStats <- nmStats[rownames(extractStats),]

fateStats <- cbind(extractStats,noMatch=nmStats)
fateStats$success <- with(fateStats, passfilter-(noMatch.singletons+noMatch.clustered))

plotData <- t(fateStats[,-(1:3)])

pdf(paste0(pargs$outdir,"readFates.pdf"),8.5,11)
op <- par(las=3,mar=c(12,4,1,1),oma=c(10,2,2,2))
plotcols=c("firebrick3","firebrick2","gold","orange","chartreuse3")
barplot(plotData,col=plotcols,border=NA,ylab="reads")
grid(NA,NULL)
legend("right",rownames(plotData),fill=plotcols,bg="white")
par(op)
dev.off()

