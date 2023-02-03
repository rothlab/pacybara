#/bin/env Rscript
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

# infiles <- list.files(pattern="Region1.*\\.csv")

# bins <- yogitools::extract.groups(infiles,"(Q\\d+)")[,1]
# maps <- lapply(infiles, function(infile) {
#   map <- read.csv(infile,comment.char="#")
#   rownames(map) <- map$hgvs_pro
#   map
# })
# names(maps) <- bins

# allMut <- Reduce(union,lapply(maps,`[`,,1))

# scores <- do.call(data.frame,lapply(maps, function(map) map[allMut,"score"]))
# colnames(scores) <- bins
# scores <- cbind(hgvs=allMut,scores)

bins <- paste0("Q",1:4)


lrs <- read.csv("allLRs.csv")
isSingleMut <- grepl("^p\\.[A-Za-z]{3}\\d+[A-Za-z]{3}$|^p\\.[A-Za-z]{3}\\d+=$",lrs$aaChangeHGVS)
lrs <- lrs[isSingleMut,]
lrvals <- lrs[,sprintf("Region1.%s.lr",bins)]
lrsds <- lrs[,sprintf("Region1.%s.sd",bins)]


nonempty <- apply(lrvals,1,function(row)all(is.finite(row)))
wellmeasured <- lrs[,"Region1.Q1.allfreq"] > 5e-7
lrs <- lrs[which(nonempty & wellmeasured),]
lrvals <- lrvals[nonempty & wellmeasured,]
lrsds <- lrsds[nonempty & wellmeasured,]

nsIdx <- grepl("Ter$",lrs$aaChangeHGVS)
synIdx <- grepl("=$",lrs$aaChangeHGVS)
misIdx <- !nsIdx & !synIdx

drawBinCurve <- function(idx,...) {
  plot(NA,type="n",
    xlim=c(1,4),ylim=range(lrvals),
    xlab="bin",ylab="LR",
    axes=FALSE,...
  )
  axis(2)
  axis(1,at=1:4,labels=bins)
  invisible(lapply(idx, function(i){
    lines(1:4,lrvals[i,],...)
  }))
  abline(h=0,lty="dashed")
}

pdf("binCurves.pdf",5,10)
op <- par(mfrow=c(3,1))
drawBinCurve(which(synIdx),
  col=yogitools::colAlpha("chartreuse3",0.2),
  main="Synonymous"
)
drawBinCurve(which(nsIdx),
  col=yogitools::colAlpha("firebrick3",0.2),
  main="Nonsense"
)
drawBinCurve(sample(which(misIdx),1000),
  col=yogitools::colAlpha("royalblue3",0.2),
  main="Missense"
)
par(op)
dev.off()

dotcol <- misIdx + 2*nsIdx + 3*synIdx
res <- 200
png("binDots.png",width=res*5,height=res*5,res=res)
pairs(lrvals,col=dotcol,pch=c(".","x","o")[dotcol],
  labels=bins
)
dev.off()


