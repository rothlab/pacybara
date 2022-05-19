#!/usr/bin/env Rscript

options(
  stringsAsFactors=FALSE,
  ignore.interactive=TRUE
)

library(argparser)
library(yogitools)
library(hash)

p <- arg_parser(
  "Merge Pacybara libraries",
  name="pacybara_merge.R"
)
p <- add_argument(p, "lib1", help="first translated clusters csv.gz file")
p <- add_argument(p, "lib2", help="second translated clusters csv.gz file")
p <- add_argument(p, "--out", help="output file name",default="merged.csv.gz")
pargs <- parse_args(p)

# pargs <- list(
# 	lib1="LDLR_R03/LDLR-R03_RQ998_clustering/clusters_transl.csv.gz",
# 	lib2="LDLR_R03_2204_pool1/R03_bc1001_RQ998_clustering/clusters_transl.csv.gz",
# 	out="merged.csv.gz"
# )

cat("Reading inputs...\n")
lib1 <- read.csv(pargs$lib1)
lib2 <- read.csv(pargs$lib2)

cat("Indexing...\n")
#index the first library
bcIdx1 <- with(lib1,tapply(1:nrow(lib1),lib1$virtualBarcode,c))
bcIdx1 <- hash(names(bcIdx1),bcIdx1)

cat("Comparing...\n")
toMerge <- list()

toAdd <- do.call(c,lapply(1:nrow(lib2), function(row) {
	vbc <- lib2[row,"virtualBarcode"]
	hits <- bcIdx1[[vbc]]
	if (is.null(hits)) {
		return(row)
	} else {
		if (lib2[row,"hgvsc"] %in% lib1[hits,"hgvsc"]) {
			#then we're already covered
			mergeHit <- which(lib1[hits,"hgvsc"] == lib2[row,"hgvsc"])
			toMerge[[length(toMerge)+1]] <<- c(row,hits[mergeHit])
			return(NULL)
		} else {
			#it's a collision
			return(row)
		}
	} 
}))

cat("Merging...\n")
#make a copy of library one to modify for output
out <- lib1

#merge identical hits
for (mergePair in toMerge) {
	row <- mergePair[[1]]
	hit <- mergePair[-1]
	if (length(hit) > 1) {
		hit <- hit[which.max(out[hit,"size"])][[1]]
	}
	out[hit,"size"] <- out[hit,"size"]+lib2[row,"size"]
	out[hit,"reads"] <- paste(out[hit,"reads"],lib2[row,"reads"],sep="|")
}

#add unique ones
out <- rbind(out,lib2[toAdd,])

cat("Tagging collisions...\n")
#re-evaluate collisions
tagCollisions <- function(bcs) {
	bcIdx <- tapply(1:length(bcs),bcs,c)
	c("","collision")[(sapply(bcIdx[bcs],length) > 1)+1]
}

out$collision <- tagCollisions(out$virtualBarcode)
out$upTagCollision <- tagCollisions(out$upBarcode)

cat("Writing output to file ",pargs$out,"\n")

con <- gzfile(pargs$out,open="w")
write.csv(out,file=con,row.names=FALSE)
close(con)

cat("Done!\n")
