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

