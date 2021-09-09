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

countfile <- commandArgs(TRUE)[[1]]

counts <- read.csv(countfile)

#re-base variant calls to correct reference sequence
rebase <- function(counts) {
	rebaseHGVS <- function(values,prefix="c",purge="1773C>T") {
		sapply(
			strsplit(gsub(paste0(prefix,"\\.|\\[|\\]"),"",values),";"), 
			function(strs) {
				out <- setdiff(strs,purge)
				if (length(out)==0) {
					paste0(prefix,".=")
				} else if (length(out)==1) {
					paste0(prefix,".",out)
				} else {
					paste0(prefix,".[",paste(out,collapse=";"),"]")
				}
			}
		)
	}
	rebaseTxt <- function(values,purge="AAC591AAT") {
		sapply(
			strsplit(values,"\\|"),
			function(strs) {
				out <- setdiff(strs,purge)
				if (length(out)==0) {
					"WT"
				} else {
					paste(out,collapse="|")
				}
			}
		)
	}
	counts$snvs <- NULL
	counts$hgvsc <- rebaseHGVS(counts$hgvsc)
	counts$hgvsp <- rebaseHGVS(counts$hgvsp,"p","Asn591=")
	counts$codonChanges <- rebaseTxt(counts$codonChanges)
	counts$codonHGVS <- rebaseHGVS(counts$codonHGVS)
	counts$aaChanges <- rebaseTxt(counts$aaChanges,"N591N")
	counts$aaChangeHGVS <- rebaseHGVS(counts$aaChangeHGVS,"p","Asn591=")
	return(counts)
}

counts <- rebase(counts)

write.csv(counts,sub("\\.csv$","_rebased.csv",countfile),row.names=FALSE)
