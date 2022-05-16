#!/usr/bin/env Rscript

#pacbio QC

options(stringsAsFactors=FALSE)


library(argparser)
# library(hgvsParseR)

p <- arg_parser(
  "Pacbio library QC",
  name="barseq_qc.R"
)
p <- add_argument(p, "library", help="lrs file")
p <- add_argument(p, "outdir", help="output directory")
pargs <- parse_args(p)

# pargs <- list(
#   library="libraries/LDLR-R03_85PM_CELL1-Cell1-CCS_minPass5_subassembly_varcalls_filtered_transl.csv",
#   outdir="workspace/qc/LDLR-R03/"
# )

#Read file
barcodes <- read.csv(pargs$library)

dir.create(pargs$outdir,showWarnings=FALSE,recursive=TRUE)


#Draw coverage plot
aaChanges <- strsplit(barcodes$aaChanges[!grepl("-|fs|WT|silent|[A-Z*]{2,}",barcodes$aaChanges)],"\\|")
aaChangeTally <- table(do.call(c,aaChanges))
aacs <- names(aaChangeTally)

toAA <- substr(aacs,nchar(aacs),nchar(aacs))
pos <- as.integer(substr(aacs,2,nchar(aacs)-1))
aaScale <- yogitools::toChars("AVLIMFYWRHKDESTNQGCP*")
y <- sapply(toAA,function(a)22-which(aaScale==a))
colvals <- colorRampPalette(c("orange","firebrick3"))(max(aaChangeTally))[aaChangeTally]

pdf(paste0(pargs$outdir,"/pb_coverage.pdf"),60,5)
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
  fs=sum(grepl("fs",barcodes$aaChanges)),
  if.indel=sum(grepl("-",barcodes$aaChanges))+sum(grepl("\\d+[A-Z*]{2,}",barcodes$aaChanges)),
  WT=sum(grepl("WT|silent",barcodes$aaChanges)),
  table(sapply(aaChanges,length))
)
plotCols=c("firebrick3","orange","gray",rep("chartreuse3",length(mutsPerClone)-3))

pdf(paste0(pargs$outdir,"/pb_census.pdf"),7,5)
op <- par(las=3)
barplot(
  mutsPerClone,
  beside=TRUE,col=plotCols,border=NA,
  xlab="#variants in a given clone",
  ylab="#clones",main=sprintf("%d useful clones out of %d passing filter",length(aaChanges),nrow(barcodes))
)
par(op)
invisible(dev.off())


cat("Done!\n")
