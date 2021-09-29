#/bin/env Rscript

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

lrs <- read.csv("allLRs.csv")
isSingleMut <- grepl("^p\\.[A-Za-z]{3}\\d+[A-Za-z]{3}$|^p\\.[A-Za-z]{3}\\d+=$",lrs$aaChangeHGVS)
lrs <- lrs[isSingleMut,]
lrvals <- lrs[,sprintf("Region1.Q%d.lr",1:4)]
lrsds <- lrs[,sprintf("Region1.Q%d.sd",1:4)]

nonempty <- apply(lrvals,1,function(row)all(is.finite(row)))
wellmeasured <- lrs[,"Region1.Q1.allfreq"] > 5e-7
lrs <- lrs[which(nonempty & wellmeasured),]
lrvals <- lrvals[nonempty,]
lrsds <- lrsds[nonempty,]

nsIdx <- grepl("Ter$",lrs$aaChangeHGVS)
synIdx <- grepl("=$",lrs$aaChangeHGVS)
misIdx <- !nsIdx & !synIdx

drawBinCurve <- function(idx,...) {
  plot(NA,type="n",
    xlim=c(1,4),ylim=range(lrvals),
    xlab="bin",ylab="LR",
    axes=FALSE
  )
  axis(2)
  axis(1,at=1:4,labels=bins)
  invisible(lapply(idx, function(i){
    lines(1:4,lrvals[i,],...)
  }))
}

op <- par(mfrow=c(3,1))
drawBinCurve(which(synIdx),col=yogitools::colAlpha("chartreuse3",0.2))
drawBinCurve(which(nsIdx),col=yogitools::colAlpha("firebrick3",0.2))
drawBinCurve(sample(which(misIdx),1000),col=yogitools::colAlpha(1,0.2))
par(op)


dotcol <- misIdx + 2*nsIdx + 3*synIdx
pairs(lrvals,col=dotcol,pch=".")
