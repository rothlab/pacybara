assay <- "Uptake"
sCond <- "F5"

scol <- sprintf("%s.%s.score",assay,sCond)
sdcol <- sprintf("%s.%s.sd",assay,sCond)
fcol <- sprintf("%s.%s.allfreq",assay,sCond)


combos <- lapply(allmuts, function(mut) {

  cat(mut,"\n")
  singles <- singleIdx[[mut]]
  multis <- setdiff(margIdx[[mut]],singles)

  singles <- singles[which(scores[singles,fcol] > 5e-7)]
  multis <- multis[which(scores[multis,fcol] > 5e-7)]

  cat("  ",length(singles)," ",length(multis),"\n")

  if (length(singles) > 0 && length(multis) > 0) {

    singleMean <- join.datapoints(
      scores[singles,scol],
      scores[singles,sdcol],
      rep(nreps,length(singles))
    )

    lapply(1:length(multis), function(n) {
      # cat(n,"\n")
      if (choose(length(multis),n)*n > 1e6) {
        matrix(sample(multis,10*n,replace=TRUE),nrow=n,ncol=10)
      } else {
        pairings <- combn(multis,n)
        if (ncol(pairings) > 10) {
          pairings <- pairings[,sample(ncol(pairings),10),drop=FALSE]
        }
      }
      
      t(apply(pairings,2,function(midxs){
        j <- join.datapoints(
          scores[midxs,scol],
          scores[midxs,sdcol],
          rep(nreps,n)
        )
        c(single=singleMean[["mj"]],multi=j[["mj"]])
      }))
    })

  } else NULL

})

# ldist <- sapply(combos,length)

multicompare <- lapply(1:20,function(n) {
  do.call(rbind,lapply(combos, function(cb) {
    if (length(cb)>= n) cb[[n]] else NULL
  }))
})

op <- par(mfrow=c(4,5)) 
for (i in 1:20) {
  pcc <- cor(yogitools::fin(multicompare[[i]]))[1,2]
  plot(multicompare[[i]],pch=".",
    xlim=c(-1,2),ylim=c(-1,2),
    main=paste("n =",i),col=yogitools::colAlpha(1,0.2)
  )
  abline(h=0:1,v=0:1,col=2:3,lty="dashed")
  text(.5,1.8,sprintf("R=%.02f",pcc),col="blue")
}
par(op)


# for (i in 1:20) {
#   png(sprintf("gif/n%02d.png",i),4*200,4*200,res=200)
#   pcc <- cor(yogitools::fin(multicompare[[i]]))[1,2]
#   plot(multicompare[[i]],pch=".",
#     xlim=c(-1,2),ylim=c(-1,2),
#     main=paste("n =",i),col=yogitools::colAlpha(1,0.2)
#   )
#   abline(h=0:1,v=0:1,col=2:3,lty="dashed")
#   text(.5,1.8,sprintf("R=%.02f",pcc),col="blue")
#   dev.off()
# }


# ffmpeg -i "n%02d.png"  -framerate 24 -filter_complex  "fps=24,split=2[palette_in][gif];[palette_in]palettegen[palette_out];[gif]fifo[gif_fifo]; [gif_fifo][palette_out]paletteuse" -y image.gif



# DOWNSAMPLING ---------------------------------

mpos <- as.integer(gsub("\\D+","",counts$hgvsp))
r1single <- which(isSingleMut & (mpos <= 187))

relfreqs <- msd$Uptake.All.mean[r1single]
relmuts <- counts$hgvsp[r1single]

fcuts <- 10^seq(-8,-3,0.2)
possible <- 21*187

downsample <- do.call(rbind,lapply(fcuts, function(fcut) {
  survivors <- relmuts[which(relfreqs > fcut)]
  out <- setNames(numeric(41),0:40)
  observ <- table(table(survivors))
  out[names(observ)] <- observ
  out[["0"]] <- possible-sum(observ)
  out
}))

cumuObs <- t(apply(downsample/possible,1,cumsum))
plotcols <- c("white",colorRampPalette(c("yellow","red","red"))(40))

plot(NA,type="n",
  xlim=range(1/fcuts),ylim=c(0,1),log="x",
  xlab="depth",ylab="coverage"
)
for (i in ncol(cumuObs):1) {
  polygon(
    c(1/fcuts[[length(fcuts)]],rev(1/fcuts),1/fcuts[[1]]),
    c(0,rev(cumuObs[,i]),0),
    border="gray",col=plotcols[[i]],lty="dotted"
  )
}
text(1/fcuts[[2]],cumuObs[2,],sprintf("%sx",colnames(cumuObs)),pos=1)


