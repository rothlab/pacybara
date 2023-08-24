#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE,ignore.interactive=TRUE)

if (as.integer(version$major) < 4) {
  stop("This script requires R version 4 or higher!")
}

library(yogitools)
library(argparser)
library(pbmcapply)

p <- arg_parser(
  "Convert cluster file to PacRat output",
  name="bara2rat.R"
)
p <- add_argument(p, "parameters", help="input directory")
p <- add_argument(p, "clustersFile", help="sample table file (*.csv.gz)")
p <- add_argument(p, "--outFile", help="output file. Defaults to <clustersFile>_consensusSeqs.txt")
pargs <- parse_args(p)
# pargs <- list(parameters="/home/jweile/projects/pacybara/templates/20230725_pacybaraparam_CYP2C9_codonOptimized.txt",clustersFile="clusters_transl_softfilter.csv.gz")

if (is.na(pargs$outFile)) {
  outfile <- paste0(sub(".csv.gz$","",pargs$clustersFile),"_consensusSeqs.txt")
} else {
  outfile <- pargs$outFile
}

cat("Parsing parameters...\n")

#read parameter sheet and extract ORF sequence from it
paramLines <- readLines(pargs$parameters)
orfStart <- as.integer(sub("^.+=","",paramLines[grep("ORFSTART=",paramLines)]))
orfEnd <- as.integer(sub("^.+=","",paramLines[grep("ORFEND=",paramLines)]))
ampRange <- grep("^#(BEGIN|END) AMPLICON SEQUENCE",paramLines)
fasta <- paramLines[(ampRange[[1]]+1):(ampRange[[2]]-1)]
ampSeq <- paste(fasta[-1],collapse="")
orfSeq <- substr(ampSeq,orfStart,orfEnd)

cat("Reading clusters file...\n")
#read cluster table
clusters <- read.csv(pargs$clustersFile)

cat("Processing...\n")
genos <- strsplit(clusters$geno,";")

mutSeqs <- do.call(c,pbmclapply(genos, function(muts) {
    #check for WT case
    if (all(muts == "=")) {
      return(orfSeq)
    }
    #extract change information
    mutdata <- data.frame(
      pos=as.integer(gsub("\\D+","",muts)),
      ref=extract.groups(muts,"\\d+([ACTG]{1})>")[,1],
      alt=extract.groups(muts,"([ACTG]+)$")[,1]
    )
    #sort in reverse order to preserve positions after deletions
    mutdata <- mutdata[order(mutdata$pos,decreasing=TRUE),]
    #make a copy of the reference sequence
    mutSeq <- orfSeq
    #implement mutations
    for (i in 1:nrow(mutdata)) {
      pos <- mutdata[i,"pos"]
      if (pos > 0 && pos <= nchar(orfSeq)) {
        if (is.na(mutdata[i,"alt"])) {#deletion
          mutSeq <- paste0(substr(mutSeq,1,pos-1),substr(mutSeq,pos+1,nchar(mutSeq)))
        } else if (is.na(mutdata[i,"ref"])) {#insertion
          mutSeq <- paste0(substr(mutSeq,1,pos-1),mutdata[i,"alt"],substr(mutSeq,pos,nchar(mutSeq)))
        } else {#substitution
          substr(mutSeq,pos,pos) <- mutdata[i,"alt"]
        }
      } else {
        warning("Ignoring off-target mutation at position",mutdata[i,"pos"])
      }
    }
    return(mutSeq)
},mc.cores=6))

out <- data.frame(barcode=clusters$upBarcode,consensus=mutSeqs)

cat("Writing output to",outfile,"\n")

write.table(out,outfile,sep="\t",row.names=FALSE,col.names=FALSE)

cat("Success!\n")

