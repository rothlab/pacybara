#!/bin/bash

FQ=$1
if [[ -z $FQ || ! $FQ =~ .fastq.gz$ ]]; then
  echo "First argument is not a FASTQ file!"
  exit 1
fi 

IDX=$2
RE='.fa$|.fasta$'
if [[ -z $IDX || ! $IDX =~ $RE ]]; then
  echo "Second argument is not a FASTA file!"
  exit 1
fi

OUT=${3:-$(basename ${FQ%.fastq.gz}_aln.png)}

#index template (just in case that wasn't done yet)
bwa index -a is "$IDX"

#run quick alignment of 1000 reads against template (filtering out non-aligned reads)
# ALN=$(basename ${FQ%.fastq.gz}_aln.sam)
ALN=$(mktemp)
bwa mem -t 8 -L 80 "$IDX" <(zcat "$FQ"|head -4000)|samtools view -F 0x4 -o $ALN -

#visualize quality scores
Rscript -e '
infile <- commandArgs(TRUE)[[1]]
outfile <- commandArgs(TRUE)[[2]]
#read sam file and filter out invalid entries
sam <- read.delim(infile,header=FALSE)
sam <- sam[!is.na(sam[,4]),]
#extract start positions and calculate plotting offset
poss <- sam[,4]
offset <- poss-min(poss,na.rm=TRUE)
#extract quality scores
qualStr <- sam[,11]
#convert to integer Q scores
quals <- lapply(qualStr,function(str) charToRaw(str)|>as.integer()-33)
#assemble into positional matrix
quals <- mapply(function(x,o)c(rep(NA,o),x),quals,offset)
maxlen <- max(sapply(quals,length))
qmat <- do.call(rbind,lapply(quals,function(x) c(x,rep(NA,maxlen-length(x)))))

#plot the data
png(outfile,4000,2000,res=200)
layout(rbind(1,2),heights=c(2,8))
op <- par(mar=c(0,4,1,1),xaxs="i")
#the top plot is a line plot of average quality scores
plot(1:maxlen,colMeans(qmat,na.rm=TRUE),type="l",axes=FALSE,ylab="mean Q",xlab="",ylim=c(70,100))
axis(2)
par(op)
op <- par(mar=c(5,4,0,1))
#the bottom is a heat map of the matrix
image(x=1:maxlen,y=1:nrow(qmat),z=t(apply(qmat,2,rev)),useRaster=TRUE,xlab="Read position",ylab="Reads")
par(op)
dev.off()|>invisible()
' "$ALN" "$OUT"

rm "$ALN"

echo "Done!"
