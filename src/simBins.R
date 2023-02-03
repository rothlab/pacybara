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

logistic <- function(x) exp(x)/(1+exp(x))

n <- 1000

trueFitness <- logistic(rnorm(n,1,3))
noisyFitness <- trueFitness + rnorm(n,0,.1)

fluorBins <- seq(-2,3,0.1)
distr <- sapply(fluorBins, function(bin) {
  sum(dnorm(bin,noisyFitness,0.5))
})
cumu <- cumsum(distr/sum(distr))

nbins <- 10
qborders <- sapply(seq(0,1,length.out=nbins+1)[c(-1,-(nbins+1))], function(q) fluorBins[[min(which(cumu > q))]])

dMat <- do.call(cbind,lapply(1:nbins, function(bin) {
  if (bin==1) {
    pnorm(qborders[[1]],noisyFitness,0.5)
  } else if (bin==nbins) {
    pnorm(qborders[[nbins-1]],noisyFitness,0.5,lower.tail=FALSE)
  } else {
    pnorm(qborders[[bin]],noisyFitness,0.5)-pnorm(qborders[[bin-1]],noisyFitness,0.5)
  }
  dnorm(bin,noisyFitness,0.5)
}))
# colnames(dMat) <- fluorBins

# dTot <- colSums(dMat)

fMat <- apply(dMat,2,function(x)x/sum(x))
lrMat <- log10(fMat)-log10(1/n)

op <- par(mfrow=c(2,1))
plot(NA,xlim=c(1,nbins),ylim=c(range(lrMat)),xlab="bin",ylab="simul LR")
for (i in which(trueFitness > .9)) {
  lines(1:nbins,lrMat[i,])
}
plot(NA,xlim=c(1,nbins),ylim=c(range(lrMat)),xlab="bin",ylab="simul LR")
for (i in which(trueFitness <.1)) {
  lines(1:nbins,lrMat[i,])
}
par(op)


