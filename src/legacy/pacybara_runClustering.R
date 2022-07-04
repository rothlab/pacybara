#!/usr/bin/env Rscript
options(
  stringsAsFactors=FALSE
)

library(pbmcapply) 
library(hash)
library(yogitools)
library(yogiseq)
library(argparser)

p <- arg_parser(
  "run clustering of barcoded pacbio data",
  name="runClustering.R"
)
p <- add_argument(p, "edFile", help="edit distance csv.gz file")
p <- add_argument(p, "genoFile", help="genotypes csv.gz file")
p <- add_argument(p, "preClustFile", help="pre-clustered virtual barcode file")
p <- add_argument(p, "--uptagBarcodeFile", help="up-tag barcode fastq.gz file")
p <- add_argument(p, "--minJaccard",help="minimum jaccard cofficient to merge clusters",default=0.2)
p <- add_argument(p, "--minMatches",help="minimum number of variant matches to merge clusters",default=1L)
p <- add_argument(p, "--maxDiff",help="maxiumum allowed edit distance between barcodes",default=2L)
p <- add_argument(p, "--minQual", help="minimum variant Q-score to be accepted as real variant.",default=100L)
p <- add_argument(p, "--verbose", help="output details about each cluster merge decision.",flag=TRUE)
p <- add_argument(p, "--out", help="output file (csv.gz)",default="clusters.csv.gz")
args <- parse_args(p)

edFile <- args$edFile
genoFile <- args$genoFile
preClustFile <- args$preClustFile
uptagFile <- args$uptagBarcodeFile
maxDiff <- args$maxDiff
minJaccard <- args$minJaccard
minMatches <- args$minMatches
minQual <- args$minQual
outfile <- args$out

options(pacybaraVerbose=args$verbose)
if (args$verbose){
  options(
    ignore.interactive=TRUE
  )
}

#function composition operator
`%o%` <- function(f,g) {
  function(...) g(f(...))
}

new.fastClust <- function(minJaccard=0.3) {

  knownPairs <- hash()
  cl2members <- hash()
  member2cl <- hash()
  sentinel <- "S"

  clCounter <- 0L
  newClusterID <- function() {
    clCounter <<- clCounter+1L
    sprintf("%012d",clCounter)
  }
  pairID <- function(id1,id2) {
    paste(sort(c(id1,id2)),collapse="-")
  }
  isKnown <- function(id1,id2) {
    has.key(pairID(id1,id2),knownPairs)
  }
  clusterOf <- function(id) {
    if (hash::has.key(id,member2cl)) {
      member2cl[[id]]
    } else {
      sentinel
    }
  }
  membersOf <- function(clID) {
    if (hash::has.key(clID,cl2members)) {
      cl2members[[clID]]
    } else {
      character(0)
    }
  }
  suggestMerge <- function(id1,id2,bcDist,genoDist,jaccard) {
    clID1 <- clusterOf(id1)
    clID2 <- clusterOf(id2)
    size1 <- max(length(membersOf(clID1)),1)
    size2 <- max(length(membersOf(clID2)),1)
    totalDist <- bcDist+genoDist

    #size disproportion. 0 = equal, 1 = double, 2 = quadruple
    sizeDispro <- abs(log2(size1/size2))
    #total edit distance = 0: accept
    #total ED 1: Min. dispro 0.25?
    #total ED 2: Min. dispro 0.5
    #total ED 3: Min. dispro 1?
    #total ED 4: reject
    #FIXME: Conservatively, both being of small size should be no excuse.
    # This just invites chaining!!
    # sizeCriterion <- sizeDispro >= 2^(bcDist-3)
    sizeCriterion <- max(size1,size2) <= 4 || sizeDispro >= 2^(bcDist-3)

    if (totalDist == 0 || (sizeCriterion && jaccard >= minJaccard)) {
      acceptMerge(id1,id2,bcDist,genoDist,jaccard)
    } else {
      rejectMerge(id1,id2,bcDist,genoDist,jaccard)
    }
  }
  acceptMerge <- function(id1,id2,bcDist,genoDist,jaccard) {
    #mark pairing as known
    knownPairs[[pairID(id1,id2)]] <<- TRUE
    #retrieve cluster IDs
    clID1 <- clusterOf(id1)
    clID2 <- clusterOf(id2)
    if (clID1==sentinel) {
      if (clID2==sentinel) {
        #both are new
        jointID <- newClusterID()
        cl2members[[jointID]] <<- c(id1,id2)
        member2cl[[id1]] <<- jointID
        member2cl[[id2]] <<- jointID
        if (getOption("pacybaraVerbose")) {
          cat(sprintf(
            "New cluster %s. Members: %s, %s. bcDist:%d; genoDist:%d; jaccard:%.02f\n",
            jointID,id1,id2,bcDist,genoDist,jaccard
          ))
        }
      } else {
        #read1 is new
        cl2members[[clID2]] <<- c(cl2members[[clID2]],id1)
        member2cl[[id1]] <<- clID2
        if (getOption("pacybaraVerbose")) {
          cat(sprintf(
            "Added %s to %s. Size:%d; bcDist:%d; genoDist:%d; jaccard:%.02f\n",
            id1,clID2,length(cl2members[[clID2]]),bcDist,genoDist,jaccard
          ))
        }
      }
    } else {
      if (clID2==sentinel) {
        #read2 is new
        cl2members[[clID1]] <<- c(cl2members[[clID1]],id2)
        member2cl[[id2]] <<- clID1
        if (getOption("pacybaraVerbose")) {
          cat(sprintf(
            "Added %s to %s. Size:%d; bcDist:%d; genoDist:%d; jaccard:%.02f\n",
            id2,clID1,length(cl2members[[clID1]]),bcDist,genoDist,jaccard
          ))
        }
      } else {
        #both are already with existing clusters
        if (clID1==clID2) {
          #already done
          if (getOption("pacybaraVerbose")) {
            cat(sprintf(
              "Rediscovered %s. Size:%d; bcDist:%d; genoDist:%d; jaccard:%.02f\n",
              clID1,length(cl2members[[clID1]]),bcDist,genoDist,jaccard
            ))
          }
        } else {
          #merge clusters under first cluster ID, delete second
          old2Members <- cl2members[[clID2]]
          cl2members[[clID1]] <<- c(cl2members[[clID1]],old2Members)
          values(member2cl,keys=old2Members) <<- clID1
          cl2members[[clID2]] <<- NULL
          if (getOption("pacybaraVerbose")) {
            cat(sprintf(
              "Merged %s and %s. Size:%d; bcDist:%d; genoDist:%d; jaccard:%.02f\n",
              clID1,clID2,length(cl2members[[clID1]]),bcDist,genoDist,jaccard
            ))
          }
        }
      }
    }
  }
  rejectMerge <- function(id1,id2,bcDist,genoDist,jaccard) {
    knownPairs[[pairID(id1,id2)]] <<- TRUE
    if (getOption("pacybaraVerbose")) {
      cat(sprintf(
        "Rejected pairing %s and %s. bcDist:%d; genoDist:%d; jaccard:%.02f\n",
        id1,id2,bcDist,genoDist,jaccard
      ))
    }
  }
  list(
    isKnown=isKnown,
    suggestMerge=suggestMerge,
    acceptMerge=acceptMerge,
    rejectMerge=rejectMerge,
    getClusters=function() cl2members,
    getMembersOf=membersOf,
    getClusterOf=function(id) {
        if (has.key(id,member2cl)) {
        member2cl[[id]]
      } else {
        NA
      }
    }
  )
}

# Helper function to build variant matrix for list of read IDs
toVarMatrix <- function(ids) {
  gl <- genoList[ids]
  if (length(ids) > 1) {
    suppressWarnings({
      Reduce(function(x,y) merge(x,y,by="var",all=TRUE), gl)
    })
  } else {
    gl[[1]]
  }
}

#Helper function to lookup any existin clusters for a read and build the
# genotype matrix
getClusterVarMatrix <- function(rid) {
  with(fastclust, {
    clid <- getClusterOf(rid)
    if (!is.na(clid)) {
      getMembersOf(clid) |> toVarMatrix()
    } else {
      rid |> toVarMatrix()
    }
  })
}

cat("Loading barcode files\n")

#LOAD BARCODES AND EXTRACT PRE-CLUSTERS
loadBarcodes <- function(preClustFile) {
  fqLines <- readLines(preClustFile)#yes, this automatically deals with gzip
  headerLines <- fqLines[seq(1,length(fqLines),4)]
  seqLines <- fqLines[seq(2,length(fqLines),4)]
  bcReadName <- yogitools::extract.groups(headerLines,"^@(\\S+) ")[,1]
  setNames(seqLines,bcReadName)
}
barcodes <- loadBarcodes(preClustFile)
#List of lists of all of read IDs with identical barcodes
preClustsAll <- names(barcodes)|>strsplit("\\|")
preClusts <- preClustsAll[sapply(preClustsAll,length) > 1]

if (!is.na(uptagFile)){
  uptags <- loadBarcodes(uptagFile)
} else {
  uptags <- NA
}


#LOAD GENOTYPES AND EDIT DISTANCES
csvgz <- function(filename,header=TRUE) {
  con <- gzfile(filename,open="r")
  on.exit({close(con)})
  return(read.csv(con,header=header))
}
cat("Loading edit distance file\n")
edist <- csvgz(edFile)
cat("Loading genotypes\n")
genos <- csvgz(genoFile,header=FALSE)
genos <- setNames(genos[,2],genos[,1])
#parse genos into list of dataframes
genoList <- strsplit(genos,";")|>pbmclapply(function(gs) {
  if (all(gs=="=")) {
    return(data.frame(var=character(0),qual=numeric(0)))
  }
  l <- strsplit(gs,":")
  data.frame(var=sapply(l,`[[`,1),qual=as.numeric(sapply(l,`[[`,2)))
},mc.cores=8)

cat("Clustering: Finding seed clusters...\n")

#create cluster object
fastclust <- new.fastClust()

#FIND SEED CLUSTERS
preClusts |> lapply(function(ids) {
  varMatrix <- toVarMatrix(ids)
  #remove putative sequencing errors
  # - Observed only once, with qual < 100
  observations <- apply(varMatrix[,-1], 1, (is.na %o% `!` %o% sum))
  totalQual <- apply(varMatrix[,-1], 1, sum, na.rm=TRUE)
  varMatrix <- varMatrix[which(observations > 1 | totalQual > minQual),]
  #pick subclusters
  for (i in 2:length(ids)) {
    for (j in 1:(i-1)) {
      vars1 <- !is.na(varMatrix[,i+1])
      vars2 <- !is.na(varMatrix[,j+1])
      totVar <- sum(vars1 | vars2)
      if (totVar == 0) {#both WT
        fastclust$acceptMerge(ids[[i]],ids[[j]],
          bcDist=0,genoDist=0,jaccard=1
        )
        next
      }
      matches <- sum(vars1 & vars2)
      mismatches <- sum(xor(vars1,vars2))
      jaccard <- if (totVar > 0) matches/totVar else 1
      if (matches >= minMatches && jaccard >= minJaccard) {
        fastclust$acceptMerge(ids[[i]],ids[[j]],
          bcDist=0,genoDist=mismatches,jaccard=jaccard
        )
      } else {
        fastclust$rejectMerge(ids[[i]],ids[[j]],
          bcDist=0,genoDist=mismatches,jaccard=jaccard
        )
      }
    }
  }
}) |> invisible()

judgeConnection <- function(id1,id2,bcDist) {
  #if the pairing was already previously examined, we can skip it
  if (fastclust$isKnown(id1,id2)) {
    return()
  }
  #check if they are are already part of clusters
  varMatrix1 <- getClusterVarMatrix(id1)
  varMatrix2 <- getClusterVarMatrix(id2)

  #TODO: check if they are already part of the same cluster and if so, skip

  #cluster sizes
  clSize1 <- ncol(varMatrix1)-1
  clSize2 <- ncol(varMatrix2)-1

  #test effects of possible merger
  combinedMatrix <- suppressWarnings({
    merge(varMatrix1,varMatrix2,by="var",all=TRUE)
  })
  #remove putative seq errors
  observations <- apply(combinedMatrix[,-1], 1, (is.na %o% `!` %o% sum))
  totalQual <- apply(combinedMatrix[,-1], 1, sum, na.rm=TRUE)
  combinedMatrix <- combinedMatrix[which(observations > 1 | totalQual > minQual),]
  
  # print(combinedMatrix)

  # #extract surviving variants from both sides
  vars1 <- combinedMatrix[,(1:clSize1)+1,drop=FALSE] |> apply(1, is.na %o% `!` %o% any)
  vars2 <- !is.na(combinedMatrix[,clSize1+2])
  totVar <- sum(vars1 | vars2)
      
  if (totVar == 0) {#both WT 
    totDist <- bcDist
    fastclust$suggestMerge(id1,id2,
      bcDist=bcDist,genoDist=0,jaccard=1
    )
  } else {
    matches <- sum(vars1 & vars2)
    mismatches <- sum(xor(vars1,vars2))
    jaccard <- if (totVar > 0) matches/totVar else 1
    # cat(sprintf("%d matches; %d mismatches; jaccard=%.03f\n\n",matches,mismatches,jaccard))
    fastclust$suggestMerge(id1,id2,
      bcDist=bcDist,genoDist=mismatches,jaccard=jaccard
    )
  }
}

#Iteratively step through connections with greater barcode differences
for (bcDist in 1:maxDiff) {
  cat("Clustering: Examining merges at barcode distance",bcDist,"...\n")
  edist[edist$dist==bcDist,1:2] |> apply(1, function(pair) {
    #the pair of pre-clusters still needs to be resolved into a list of reads
    ids <- strsplit(pair,"\\|")
    for (id1 in ids[[1]]) {
      for (id2 in ids[[2]]) {
        judgeConnection(id1,id2,bcDist)
      }
    }
  }) |> invisible()
}


#CALCULATE CONSENSUS GENOTYPES
cat("Calculating consensus genotypes\n")

clusterGenos <- fastclust$getClusters()|>as.list()|>pbmclapply(function(members) {
  varMatrix <- toVarMatrix(members)
  #remove "obvious" sequencing error
  observations <- apply(varMatrix[,-1], 1, (is.na %o% `!` %o% sum))
  totalQual <- apply(varMatrix[,-1], 1, sum, na.rm=TRUE)
  varMatrix <- varMatrix[which(observations > 1 | totalQual > minQual),]
  #if no variants survive, we're done
  if (nrow(varMatrix)==0) {
    return("=")
  } 
  #assign WT voting weights
  vars <- varMatrix$var
  weights <- as.matrix(varMatrix[,-1])
  weights[is.na(weights)] <- -minQual
  #filter based on votes
  passFilter <- rowSums(weights) > 0
  # if (any(!passFilter)) {
  #   print(varMatrix)
  # }
  #if no variants survive, we're done
  if (sum(passFilter)==0) {
    return("=")
  }
  #sort the rest and output
  passedVars <- vars[passFilter]
  poss <- as.integer(yogitools::extract.groups(passedVars,"^(\\d+)"))
  # poss <- as.integer(sub("\\D+","",varMatrix$var))
  return(paste(passedVars[order(poss)],collapse=";"))
},mc.cores=8)|>unlist()


#CALCULATE CONSENSUS BARCODES

#helper function to calculate a consensus sequence via MUSCLE
alignmentConsensus <- function(bcs) {
  fastaFile <- tempfile()
  alnFile <- tempfile()
  con <- file(fastaFile,open="w")
  yogiseq::writeFASTA(con,bcs)
  close(con)
  retVal <- system2("muscle",
    args=c(
      "-in",fastaFile,
      "-out",alnFile
    ),
    stdout=FALSE,
    stderr=FALSE
  )
  if (retVal == 0) {
    con <- file(alnFile,open="r")
    alnLines <- yogiseq::readFASTA(con)|>sapply(function(x)x$toString())
    close(con)
    out <- paste(sapply(1:nchar(alnLines[[1]]), function(i) {
      contab <- table(substr(alnLines,i,i))
      outchar <- names(which.max(contab))
      if (outchar == "-") "" else outchar
    }),collapse="")
  } else {
    out <- "alignment_error"
  }
  file.remove(fastaFile,alnFile)
  return(out)
}

calcConsensus <- function(rid2bc) {
  fastclust$getClusters()|>as.list()|>sapply(function(members) {
    bcs <- values(rid2bc,members)
    contable <- table(bcs)
    if (sum(contable==max(contable))==1) {
      return(names(which.max(contable)))
    } else {
      return(alignmentConsensus(bcs))
    }
  })
}

cat("Calculating consensus barcodes\n")

#Calc consensus for virtual barcodes
rid2vbc <- hash(
  preClustsAll|>unlist(),
  lapply(1:length(barcodes),function(i) rep(barcodes[[i]],length(preClustsAll[[i]]))) |>unlist()
)
bcConsensus <- calcConsensus(rid2vbc)

#Calc consensus for uptags if available)
if (!is.na(uptags[[1]])) {
  rid2upbc <- hash(names(uptags),uptags)
  upConsensus <- calcConsensus(rid2upbc)
} else {
  upConsensus <- NA_character_
}


cat("Compiling results\n")

#identify singleton reads
clusteredReadIdx <- hash(fastclust$getClusters()|>as.list()|>unlist(),TRUE)
allReads <- unlist(preClustsAll)
singletonReads <- allReads[which(sapply(
  allReads,
  function(rid) is.null(clusteredReadIdx[[rid]])
))]
singletonGenos <- genoList[singletonReads]|>sapply(function(tbl){
  if (nrow(tbl)==0) return("=")
  outvar <- with(tbl,var[qual >= minQual])
  if (length(outvar)==0) return("=")
  return(paste(outvar,collapse=";"))
})

#build output table
out1 <- lapply(names(clusterGenos),function(cname) {
  rids <- fastclust$getClusters()[[cname]]
  list(
    virtualBarcode=bcConsensus[[cname]],
    upBarcode=if(is.na(upConsensus[[1]])) NA_character_ else upConsensus[[cname]],
    reads=paste(rids,collapse="|"),
    size=length(rids),
    geno=clusterGenos[cname]
  )
})|>yogitools::as.df()
out2 <- data.frame(
  # virtualBarcode=barcodes[singletonReads],
  virtualBarcode=sapply(singletonReads,function(rid) {
    bc <- rid2vbc[[rid]]
    if (is.null(bc)) NA_character_ else bc
  }),
  upBarcode=if(is.na(uptags[[1]])) NA_character_ else uptags[singletonReads],
  reads=singletonReads,
  size=1,
  geno=singletonGenos
)
out <- rbind(out1,out2)

cat("Tagging barcode collisions\n")

#check for barcode collisions
tagDups <- function(bcs) {
  contab <- table(bcs)
  dupIdx <- hash(names(which(contab > 1)),TRUE)
  bcs|>sapply(function(bc) if (is.null(dupIdx[[bc]])) "" else "collision")
}
out$collision <- tagDups(out$virtualBarcode)
if (!is.na(uptagFile)) {
  out$upTagCollision <- tagDups(out$upBarcode)
}

cat("Writing output to file\n")

#write output to file
outcon <- gzfile(outfile,open="w")
write.csv(out,outcon,row.names=FALSE)
close(outcon)

# inspectConnection <- function(id1,id2,bcDist) {
#   #if the pairing was already previously examined, we can skip it
#   if (fastclust$isKnown(id1,id2)) {
#     return()
#   }
#   #check if they are are already part of clusters
#   varMatrix1 <- getClusterVarMatrix(id1)
#   varMatrix2 <- getClusterVarMatrix(id2)

#   #cluster sizes
#   clSize1 <- ncol(varMatrix1)-1
#   clSize2 <- ncol(varMatrix2)-1

#   #test effects of possible merger
#   combinedMatrix <- suppressWarnings({
#     merge(varMatrix1,varMatrix2,by="var",all=TRUE)
#   })
#   #remove putative seq errors
#   observations <- apply(combinedMatrix[,-1], 1, (is.na %o% `!` %o% sum))
#   totalQual <- apply(combinedMatrix[,-1], 1, sum, na.rm=TRUE)
#   combinedMatrix <- combinedMatrix[which(observations > 1 | totalQual > minQual),]
  
#   print(combinedMatrix)

#   # #extract surviving variants from both sides
#   vars1 <- combinedMatrix[,(1:clSize1)+1,drop=FALSE] |> apply(1, is.na %o% `!` %o% any)
#   vars2 <- !is.na(combinedMatrix[,clSize1+2])
#   totVar <- sum(vars1 | vars2)
      
#   if (totVar > 0) {
#     matches <- sum(vars1 & vars2)
#     mismatches <- sum(xor(vars1,vars2))
#     jaccard <- if (totVar > 0) matches/totVar else 1
#     cat(sprintf(
#       "%d matches; %d mismatches; jaccard=%.03f\n\n",
#       matches,mismatches,jaccard
#     ))
#   }
# }

# bcDist<-2
# edist[edist$dist==bcDist,1:2] |> head() |> apply(1, function(pair) {
#   #the pair of pre-clusters still needs to be resolved into a list of reads
#   ids <- strsplit(pair,"\\|")
#   for (id1 in ids[[1]]) {
#     for (id2 in ids[[2]]) {
#       inspectConnection(id1,id2,bcDist)
#     }
#   }
# }) |> invisible()

