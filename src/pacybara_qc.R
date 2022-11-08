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
library(pbmcapply)

p <- arg_parser(
  "pacybara QC",
  name="pacybara_qc.R"
)
p <- add_argument(p, "clusters", help="translated clusters csv.gz file")
p <- add_argument(p, "extractDir", help="pacybara extract directory")
p <- add_argument(p, "outdir", help="output directory")
p <- add_argument(p, "--softFilter", help="indicates that the clusters file is soft-filtered.", flag=TRUE)
pargs <- parse_args(p)

# pargs <- list(clusters="clusters_transl.csv.gz",outdir="qc/")

if (!grepl("/$",pargs$outdir)) {
  pargs$outdir <- paste0(pargs$outdir,"/")
}
dir.create(pargs$outdir,showWarnings=FALSE,recursive=TRUE)

clusters <- read.csv(pargs$clusters)


#######################################
#analyze distribution of barcode sizes
uptagDistro <- clusters$upBarcode|>table()|>table()
sizeDistro <- clusters$size|>table()

if (!pargs$softFilter) {
  sizes <- 1:(c(uptagDistro|>names(),sizeDistro|>names())|>as.integer()|>max())
  distros <- sizes|>as.character()|>sapply(\(x)c(up=uptagDistro[x],clust=sizeDistro[x]))
  distros[is.na(distros)] <- 0

  if (max(sizes) > 20) {
    distros <- cbind(
      distros[,1:19],
      `>=20`=rowSums(distros[,20:ncol(distros)])
    )
  }

  pdf(paste0(pargs$outdir,"clusterSizes.pdf"),7,5)
  opar <- par(las=3)
  distros|>barplot(beside=TRUE,border=NA,
    xlab="size",ylab="count",
    col=c("steelblue3","chartreuse3")
  )
  par(opar)
  legend("topright",c("uptag counts","final cluster sizes"),fill=c("steelblue3","chartreuse3"))
  dev.off()
  
} else {
  sizes <- 1:(sizeDistro|>names()|>as.integer()|>max())
  sizeDistro <- setNames(sizeDistro[as.character(1:max(sizes))],1:max(sizes))
  sizeDistro[is.na(sizeDistro)] <- 0
  if (max(sizes) > 20) {
    sizeDistro <- c(
      sizeDistro[1:19],
      `>=20`=sum(sizeDistro[20:length(sizeDistro)])
    )
  }
  pdf(paste0(pargs$outdir,"clusterSizes_softFilter.pdf"),7,5)
  opar <- par(las=3)
  sizeDistro|>barplot(border=NA,
    xlab="size",ylab="count",
    col="chartreuse3"
  )
  par(opar)
  # legend("topright",c("uptag counts","final cluster sizes"),fill=c("steelblue3","chartreuse3"))
  dev.off()
}

if (!pargs$softFilter) {
  cleanBarcodes <- with(clusters, upBarcode[upTagCollision=="" & size > 1] |> unique() )
  compromisedBarcodes <- with(clusters, upBarcode[upTagCollision=="collision"] |> unique() )

  pdf(paste0(pargs$outdir,"compromisedBC.pdf"),5,5)
  barplot(
    c(clean=length(cleanBarcodes),compromised=length(compromisedBarcodes)),
    border=NA,col=c("steelblue3","firebrick3"),ylab="#uptag barcodes"
  )
  dev.off()
}

#########################################
# Filter for usable barcodes and CCS > 1
if (pargs$softFilter) {
  clusters <- clusters[clusters$size > 1,]
} else {
  clusters <- clusters[with(clusters,upTagCollision=="" & size > 1),]
}

##########################################
# Draw coverage plot
cat("Drawing coverage plot...\n")
aaChanges <- strsplit(clusters$aaChanges[!grepl("-|fs|WT|silent|[A-Z*]{2,}",clusters$aaChanges)],"\\|")
aaChangeTally <- table(do.call(c,aaChanges))
aacs <- names(aaChangeTally)

toAA <- substr(aacs,nchar(aacs),nchar(aacs))
pos <- as.integer(substr(aacs,2,nchar(aacs)-1))
aaScale <- yogitools::toChars("AVLIMFYWRHKDESTNQGCP*")
y <- sapply(toAA,function(a)22-which(aaScale==a))
colvals <- colorRampPalette(c("orange","firebrick3"))(max(aaChangeTally))[aaChangeTally]

pdf(paste0(pargs$outdir,"pb_coverage.pdf"),60,5)
plot(NA,type="n",xlim=c(0,max(pos)+1),ylim=c(0,22),
  axes=FALSE,xlab="position",ylab="AA",xaxs="i"
)
axis(1)
axis(2,21:1,aaScale)
rect(.5,.5,max(pos)+.5,21.5,col="gray",border=NA)
rect(pos-.5,y-.5,pos+.5,y+.5,col=colvals,border=NA)
invisible(dev.off())


#draw Census plot
cat("Drawing census plots...\n")
mutsPerClone <- c(
  fs=sum(grepl("fs",clusters$aaChanges)),
  if.indel=sum(grepl("-",clusters$aaChanges))+sum(grepl("\\d+[A-Z*]{2,}",clusters$aaChanges)),
  WT=sum(grepl("WT|silent",clusters$aaChanges)),
  table(sapply(aaChanges,length))
)
plotCols=c("firebrick3","orange","gray",rep("chartreuse3",length(mutsPerClone)-3))

pdf(paste0(pargs$outdir,"/pb_census.pdf"),7,5)
op <- par(las=3)
barplot(
  mutsPerClone,
  beside=TRUE,col=plotCols,border=NA,
  xlab="#variants in a given clone",
  ylab="#clones",main=sprintf("%d useful clones out of %d passing filter",length(aaChanges),nrow(clusters))
)
par(op)
invisible(dev.off())

#extraction QC
cat("Performing extraction QC...\n")
successCon <- pipe(sprintf("zcat %s/genoExtract.csv.gz|wc -l",pargs$extractDir))
success <- readLines(successCon)
close(successCon)
success <- as.integer(success)
exceptions <- readLines(sprintf("%s/exceptions.txt",pargs$extractDir))
exceptions <- strsplit(exceptions,"=")
extractQC <- sapply(exceptions,`[[`,2)|>as.integer()|>setNames(sapply(exceptions,`[[`,1))
extractQC <- c(success=success,extractQC)

pdf(sprintf("%s/extractionQC.pdf",pargs$outdir),5,7)
# png(sprintf("%s/extractionQC.png",pargs$outdir),200*5,200*7,res=200)
op <- par(mar=c(10,4,4,1),las=3)
barplot(extractQC,border=NA,col=c(3,2,2,2),ylab="CCS reads",main="barcode and genotype extraction")
par(op)
dev.off()


#### JACKPOT PLOT
cat("Performing jackpot QC...\n")
cvars <- gsub("c\\.|\\[|\\]$","",clusters$hgvsc)|>strsplit(";")
cvarTally <- unlist(cvars)|>table()
cvarTally <- cvarTally[-which(names(cvarTally)=="=")]|>sort(decreasing=TRUE)

pdf(sprintf("%s/jackpots.pdf",pargs$outdir),10,5)
cmap <- yogitools::colmap(
  c(0,2,4,nrow(clusters)/1e4,nrow(clusters)/5e3,nrow(clusters)/1e3),
  c("firebrick3","gold","chartreuse3","chartreuse3","gold","firebrick3")
)
plotcol <- cmap(cvarTally)
plot(
  seq(0,1,length.out=length(cvarTally)),as.vector(cvarTally),
  xlab="fraction of variants",ylab="#clones",col=plotcol,pch=20
)
topVars <- cvarTally[cvarTally > nrow(clusters)/5e3]|>head(20)
labelx <- 1:length(topVars)%%5/6
segments((1:20)/length(cvarTally),topVars,labelx,col="gray",lty="dotted")
text(labelx,topVars,names(topVars),pos=4,cex=.7)
dev.off()


#### Nucleotide bias plot
cat("Analyzing nucleotide change biases...\n")
ccTable <- clusters$codonChanges[!grepl("WT|indel|silent",clusters$codonChanges)]|>
  strsplit("\\|")|>unlist()|>na.omit()|>
  yogitools::extract.groups("([ACGT]{3})(\\d+)([ACGT]+)")|>
  as.data.frame()
colnames(ccTable) <- c("from","pos","to")
ccTable$pos <- as.integer(ccTable$pos)
changes <- ccTable[,c(1,3)] |> apply(1,function(row) {
  ks <- sapply(1:3,function(k){
    substr(row[[1]],k,k)!=substr(row[[2]],k,k)
  }) |> which()
  lapply(ks, function(k) {
    list(k,substr(row[[1]],k,k),substr(row[[2]],k,k))
  }) 
},simplify=FALSE)
nchanges <- sapply(changes,length)
singleChanges <- do.call(c,changes[nchanges==1])
multiChanges <- do.call(c,changes[nchanges>1])

countChanges <- function(changeList) {
  mat <- array(0, dim=list(4,4,3),
    dimnames=list(c("A","C","G","T"),c("A","C","G","T"),1:3)
  )
  lapply(changeList, function(ch) {
    mat[ch[[2]],ch[[3]],ch[[1]]] <<- mat[ch[[2]],ch[[3]],ch[[1]]]+1
  }) |> invisible()
  apply(mat,c(1,3),function(x)x/sum(x))
}

singleMat <- countChanges(singleChanges)
multiMat <- countChanges(multiChanges)
cmap <- yogitools::colmap(c(0,1),c("white","steelblue3"))

drawMat <- function(mat,main="") {
  plot(NA,type="n",xlim=c(0,4),ylim=c(0,4),
    xlab="from",ylab="to",axes=FALSE,main=main
  )
  axis(1,at=(1:4)-.5,c("A","C","G","T"))
  axis(2,at=(4:1)-.5,c("A","C","G","T"))
  xs <- rep(1:4,each=4)
  ys <- rep(4:1,4)
  rect(xs-1,ys-1,xs,ys,col=cmap(mat),border=NA)
  text(xs-.5,ys-.5,sprintf("%.01f%%",mat*100))
}

pdf(sprintf("%s/nuclBias.pdf",pargs$outdir),9,6)
layout(cbind(1:2,3:4,5:6))
for (i in 1:3) {
  drawMat(singleMat[,,i],paste("SNP pos.",i))
  drawMat(multiMat[,,i],paste("POP pos.",i))
}
dev.off()


### Barcode bias plots
parseFASTQ <- function(incon) {
  seqs <- character(0)
  while(length(lines <- readLines(incon,4000L)) > 0) {
    seqs <- c(seqs,lines[seq(2,length(lines),4)])
  }
  return(seqs)
}
drawBigMat <- function(mat,main="",cmap=yogitools::colmap(c(0,0.5,1),c("white","steelblue3","steelblue4"))) {
  plot(NA,type="n",xlim=c(0,ncol(mat)),ylim=c(0,nrow(mat)),
    xlab="barcode nucleotide position",ylab="base",axes=FALSE,main=main
  )
  axis(1,at=(1:ncol(mat))-.5,1:ncol(mat))
  axis(2,at=(5:1)-.5,c("A","C","G","T","-"))
  xs <- rep(1:ncol(mat),each=5)
  ys <- rep(5:1,ncol(mat))
  rect(xs-1,ys-1,xs,ys,col=cmap(mat),border=NA)
  text(xs-.5,ys-.5,sprintf("%.01f%%",mat*100),cex=0.7)
}

cat("Analyzing barcode base biases...\n")

fqfiles <- list.files(pargs$extractDir,pattern="bcExtract_\\d+\\.fastq\\.gz",full.names=TRUE)

for (fqfile in fqfiles) {
  bcname <- gsub("bcExtract_|.fastq.gz$","",basename(fqfile))

  cat("--> Processing barcode ",bcname,"\n")

  incon <- gzfile(fqfile,open="r")
  bcSeqs <- parseFASTQ(incon)
  close(incon)
  maxLen <- quantile(nchar(bcSeqs),.99)
  freqs <- do.call(cbind,pbmclapply(1:maxLen,function(i) {
    table(factor(sapply(bcSeqs,substr,i,i),levels=c("A","C","G","T","")))
  },mc.cores=8))

  pdf(sprintf("%s/barcodeBias_%s.pdf",pargs$outdir,bcname),14,4)
  drawBigMat(freqs/length(bcSeqs),main=sprintf("Barcode %s nucleotide bias",bcname))
  invisible(dev.off())

  # plot(NA,type="n",xlim=c(0,ncol(freqs)),ylim=c(0,1))
  # grid(NA,NULL)
  # for (i in 1:nrow(freqs)) {
  #   lines(1:ncol(freqs),freqs[i,]/length(bcSeqs),col=i)
  # }
}

cat("Done!\n")
