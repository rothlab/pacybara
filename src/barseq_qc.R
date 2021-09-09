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

#plot AllFreq vs Stdev
with(lrs,plot(Uptake.F5.allfreq,Uptake.F5.sd,log="x",pch=20,col=yogitools::colAlpha(1,0.2)))


hqLRs <- lrs[which(lrs$Uptake.F5.sd < 0.6 & lrs$Uptake.F5.allfreq > 1e-6),]

isSyn <- sapply(strsplit(gsub("p\\.|\\[|\\]","",hqLRs$hgvsp),";"), function(muts) {
	all(grepl("=$",muts))
})

wtLRs <- hqLRs[hqLRs$codonChanges=="WT",]
stopLRs <- hqLRs[grepl("Ter",hqLRs$hgvsp),]
fsLRs <- hqLRs[grepl("fs",hqLRs$aaChanges),]
synLRs <- hqLRs[isSyn,]
misLRs <- hqLRs[!(isSyn | hqLRs$codonChanges=="WT" | grepl("Ter",hqLRs$hgvsp) | grepl("fs",hqLRs$aaChanges)),]


pdf("barseqDistributions_UptF5.pdf",7,6)
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
dev.off()


stops <- unlist(sapply(strsplit(stopLRs[which(stopLRs$Surface.F4.lr > -0.5),"aaChanges"],"\\|"),function(muts){
	muts <- muts[grep("\\*$",muts)]
	pos <- as.integer(gsub("\\D+","",muts))
	muts[which.min(pos)]
}))
sort(table(stops),decreasing=TRUE)
