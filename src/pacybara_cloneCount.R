#!/usr/bin/env Rscript
options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE
)

library(argparser)
library(yogitools)

p <- arg_parser(
  "pacybara QC",
  name="pacybara_cloneCount.R"
)
p <- add_argument(p, "clusters", help="translated clusters csv.gz file")
p <- add_argument(p, "--regionStart",help="AA start position of mutagenesis region.",default=706L)
p <- add_argument(p, "--regionEnd",help="AA start position of mutagenesis region.",default=879L)
pargs <- parse_args(p)

regionStart <- as.integer(pargs$regionStart)
regionEnd <- as.integer(pargs$regionEnd)

clusters <- read.csv(pargs$clusters)
#filter for usable clones
clusters <- clusters[with(clusters,upTagCollision=="" & size > 1),]

aaChanges <- strsplit(clusters$aaChanges[!grepl("-|fs|WT|silent|[A-Z*]{2,}",clusters$aaChanges)],"\\|")
aaChangeTally <- table(do.call(c,aaChanges))
aacs <- names(aaChangeTally)

tallyPos <- as.integer(gsub("\\D+","",names(aaChangeTally)))
inReg <- tallyPos >= regionStart & tallyPos <= regionEnd

cat(sprintf("On-target variants: %d\n",length(aaChangeTally[inReg])))
table(aaChangeTally[inReg])

cat(sprintf("Off-target variants: %d\n",length(aaChangeTally[!inReg])))
table(aaChangeTally[!inReg])

