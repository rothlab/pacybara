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
  "draw QC plots for barseq",
  name="barseq_qc.R"
)
p <- add_argument(p, "lrs", help="lrs file")
p <- add_argument(p, "counts", help="counts file")
p <- add_argument(p, "outdir", help="output directory")
pargs <- parse_args(p)

dir.create(pargs$outdir,recursive=TRUE)

lrs <- read.csv(pargs$lrs)
counts <- read.csv(pargs$counts)


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
samples <- unique(sub("\\.allfreq","",colnames(lrs)[grep("allfreq$",colnames(lrs))]))
pdf(paste0(pargs$outdir,"freqCV.pdf"),8.5,11)
op <- par(mfrow=c(2,2),oma=c(2,2,2,2))
for (sample in samples) {
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





hqLRs <- lrs[which(lrs$Uptake.F5.allfreq > 1e-6),]

isSyn <- sapply(strsplit(gsub("p\\.|\\[|\\]","",hqLRs$hgvsp),";"), function(muts) {
  all(grepl("=$",muts))
})

wtLRs <- hqLRs[hqLRs$codonChanges=="WT",]
stopLRs <- hqLRs[grepl("Ter",hqLRs$hgvsp),]
fsLRs <- hqLRs[grepl("fs",hqLRs$aaChanges),]
synLRs <- hqLRs[isSyn,]
misLRs <- hqLRs[!(isSyn | hqLRs$codonChanges=="WT" | grepl("Ter",hqLRs$hgvsp) | grepl("fs",hqLRs$aaChanges)),]


pdf(paste0(pargs$outdir,"barseqDistributions_UptF5.pdf"),7,6)
op <- par(mfcol=c(3,2))
breaks <- seq(-5,5,0.05)
hist(
  wtLRs$Uptake.F5.lr,breaks=breaks,
  col="darkolivegreen3",border=NA,main="WT clones",
  xlab="log10(F4/All)"
)
abline(v=mean(yogitools::fin(wtLRs$Uptake.F5.lr)),lty="dotted")
hist(
  stopLRs$Uptake.F5.lr,breaks=breaks,
  col="firebrick3",border=NA,main="Nonsense clones",
  xlab="log10(F4/All)"
)
abline(v=mean(yogitools::fin(stopLRs$Uptake.F5.lr)),lty="dotted")
hist(
  synLRs$Uptake.F5.lr,breaks=breaks,
  col="darkolivegreen2",border=NA,main="Synonymous clones",
  xlab="log10(F4/All)"
)
abline(v=mean(yogitools::fin(synLRs$Uptake.F5.lr)),lty="dotted")
hist(
  fsLRs$Uptake.F5.lr,breaks=breaks,
  col="firebrick4",border=NA,main="Frameshift clones",
  xlab="log10(F4/All)"
)
abline(v=mean(yogitools::fin(fsLRs$Uptake.F5.lr)),lty="dotted")
hist(
  misLRs$Uptake.F5.lr,breaks=breaks,
  col="gray",border=NA,main="Missense clones",
  xlab="log10(F4/All)"
)
abline(v=mean(yogitools::fin(misLRs$Uptake.F5.lr)),lty="dotted")
par(op)
invisible(dev.off())


stops <- unlist(sapply(strsplit(stopLRs[which(stopLRs$Surface.F4.lr > -0.5),"aaChanges"],"\\|"),function(muts){
  muts <- muts[grep("\\*$",muts)]
  pos <- as.integer(gsub("\\D+","",muts))
  muts[which.min(pos)]
}))
sort(table(stops),decreasing=TRUE)
