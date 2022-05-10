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

p <- arg_parser(
  "pacybara QC",
  name="pacybara_qc.R"
)
p <- add_argument(p, "clusters", help="translated clusters csv.gz file")
p <- add_argument(p, "extractDir", help="pacybara extract directory")
p <- add_argument(p, "outdir", help="output directory")
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

sizes <- 1:(uptagDistro|>names()|>as.integer()|>max())
distros <- sizes|>as.character()|>sapply(\(x)c(up=uptagDistro[x],clust=sizeDistro[x]))
distros[is.na(distros)] <- 0

pdf(paste0(pargs$outdir,"clusterSizes.pdf"),5,5)
distros|>barplot(beside=TRUE,border=NA,
  xlab="size",ylab="count",
  col=c("steelblue3","chartreuse3")
)
legend("topright",c("uptag counts","final cluster sizes"),fill=c("steelblue3","chartreuse3"))
dev.off()

#########################################
# Filter for usable barcodes and CCS > 1
cleanBarcodes <- with(clusters, upBarcode[upTagCollision=="" & size > 1] |> unique() )
compromisedBarcodes <- with(clusters, upBarcode[upTagCollision=="collision"] |> unique() )

pdf(paste0(pargs$outdir,"compromisedBC.pdf"),5,5)
barplot(
  c(clean=length(cleanBarcodes),compromised=length(compromisedBarcodes)),
  border=NA,col=c("steelblue3","firebrick3"),ylab="#uptag barcodes"
)
dev.off()

clusters <- clusters[with(clusters,upTagCollision=="" & size > 1),]

##########################################
# Draw coverage plot
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
successCon <- pipe(sprintf("zcat %s/genoExtract.csv.gz|wc -l",pargs$extractDir))
success <- readLines(successCon)
close(successCon)
success <- as.integer(success)
exceptions <- readLines(sprintf("%s/exceptions.txt",pargs$extractDir))
exceptions <- strsplit(exceptions,"=")
extractQC <- sapply(exceptions,`[[`,2)|>as.integer()|>setNames(sapply(exceptions,`[[`,1))
extractQC <- c(success=success,extractQC)

pdf(sprintf("%s/extractionQC.png",pargs$outdir),5,7)
# png(sprintf("%s/extractionQC.png",pargs$outdir),200*5,200*7,res=200)
op <- par(mar=c(10,4,4,1),las=3)
barplot(extractQC,border=NA,col=c(3,2,2,2),ylab="CCS reads",main="barcode and genotype extraction")
par(op)
dev.off()

# tallyPos <- gsub("\\D+","",names(aaChangeTally))|>as.integer()
# inReg <- tallyPos >=706 & tallyPos <= 879
# aaChangeTally[inReg]|>table()

cat("Done!\n")
