#!/usr/bin/Rscript

library(pbmcapply) 
library(hash)

#function composition operator
`%o%` <- function(f,g) {
  function(...) g(f(...))
}

new.fastClust <- function(maxRatio=2,minJaccard=0.8,maxDiff=3) {

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
    if (has.key(id,member2cl)) {
      member2cl[[id]]
    } else {
      sentinel
    }
  }
  membersOf <- function(clID) {
    if (has.key(clID,cl2members)) {
      cl2members[[clID]]
    } else {
      character(0)
    }
  }
  suggestMerge <- function(id1,id2,bcDist,genoDist,jaccard) {
    clID1 <- clusterOf(id1)
    clID2 <- clusterOf(id2)
    size1 <- length(membersOf(clID1))
    size2 <- length(membersOf(clID2))
    totalDist <- bcDist+genoDist
    sizeRatio <- max(size1,size2)/min(size1,size2)

    if ((genoDist <= maxDiff && jaccard >= minJaccard) && 
         (totalDist == 0 || max(size1,size2) < 4 || sizeRatio <= maxRatio)) {
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
        cat(sprintf(
          "New cluster %s. Members: %s, %s. bcDist:%d; genoDist:%d; jaccard:%.02f\n",
          jointID,id1,id2,bcDist,genoDist,jaccard
        ))
      } else {
        #read1 is new
        cl2members[[clID2]] <<- c(cl2members[[clID2]],id1)
        member2cl[[id1]] <<- clID2
        cat(sprintf(
          "Added %s to %s. Size:%d; bcDist:%d; genoDist:%d; jaccard:%.02f\n",
          id1,clID2,length(cl2members[[clID2]]),bcDist,genoDist,jaccard
        ))
      }
    } else {
      if (clID2==sentinel) {
        #read2 is new
        cl2members[[clID1]] <<- c(cl2members[[clID1]],id2)
        member2cl[[id2]] <<- clID1
        cat(sprintf(
          "Added %s to %s. Size:%d; bcDist:%d; genoDist:%d; jaccard:%.02f\n",
          id2,clID1,length(cl2members[[clID1]]),bcDist,genoDist,jaccard
        ))
      } else {
        #both are already with existing clusters
        if (clID1==clID2) {
          #already done
          cat(sprintf(
            "Rediscovered %s. Size:%d; bcDist:%d; genoDist:%d; jaccard:%.02f\n",
            clID1,length(cl2members[[clID1]]),bcDist,genoDist,jaccard
          ))
        } else {
          #merge clusterss under first cluster ID, delete second
          old2Members <- cl2members[[clID2]]
          cl2members[[clID1]] <<- c(cl2members[[clID1]],old2Members)
          values(member2cl,keys=old2Members) <<- clID1
          cl2members[[clID2]] <<- NULL
          cat(sprintf(
            "Merged %s and %s. Size:%d; bcDist:%d; genoDist:%d; jaccard:%.02f\n",
            clID1,clID2,length(cl2members[[clID1]]),bcDist,genoDist,jaccard
          ))
        }
      }
    }
  }
  rejectMerge <- function(id1,id2,bcDist,genoDist,jaccard) {
    knownPairs[[pairID(id1,id2)]] <<- TRUE
    cat(sprintf(
      "Rejected pairing %s and %s. bcDist:%d; genoDist:%d; jaccard:%.02f\n",
      id1,id2,bcDist,genoDist,jaccard
    ))
  }
  list(
    isKnown=isKnown,
    suggestMerge=suggestMerge,
    acceptMerge=acceptMerge,
    rejectMerge=rejectMerge,
    getClusters=function() cl2members
  )
}

# fc <- new.fastClust()
# fc$acceptMerge("foo","bar",0,0,1)
# fc$acceptMerge("foo","baz",0,0,1)
# fc$acceptMerge("aaa","bbb",0,0,1)
# fc$acceptMerge("foo","aaa",0,0,1)

edFile <- "editDistance.csv.gz"
genoFile <- "genoExtract.csv.gz"

csvgz <- function(filename,header=TRUE) {
  con <- gzfile(filename,open="r")
  on.exit({close(con)})
  return(read.csv(con,header=header))
}

genos <- csvgz(genoFile,header=FALSE)
genos <- setNames(genos[,2],genos[,1])
#parse genos into list of dataframes
genoList <- pbmclapply(strsplit(genos,";"),function(gs) {
  if (all(gs=="=")) {
    return(data.frame(var=character(0),qual=numeric(0)))
  }
  l <- strsplit(gs,":")
  data.frame(var=sapply(l,`[[`,1),qual=as.numeric(sapply(l,`[[`,2)))
},mc.cores=8)

edist <- csvgz(edFile)

#get list of (distance-0) pre-clusters
preIDs <- unique(c(edist$ref,edist$edist$query))
preClusts <- strsplit(preIDs,"\\|")
preClusts <- preClusts[sapply(preClusts,length) > 1]

#create cluster object
fastclust <- new.fastClust()

# invisible(lapply(preClusts,function(ids) {
#   apply(combn(ids,2),2,function(idpair) {
#     judgeConnection(idpair[[1]],idpair[[2]],0)
#   })
# }))

minMatches <- 1
minJaccard <- 0.2

#examine seed cluster candidates
preClusts |> lapply(function(ids) {
  gl <- genoList[ids]
  varMatrix <- Reduce(function(x,y) merge(x,y,by="var",all=TRUE), gl)
  #remove putative sequencing errors
  # - Observed only once, with qual < 100
  observations <- apply(varMatrix[,-1], 1, (is.na %o% `!` %o% sum))
  totalQual <- apply(varMatrix[,-1], 1, sum, na.rm=TRUE)
  varMatrix <- varMatrix[which(observations > 1 | totalQual > 100),]
  #pick subclusters
  for (i in 2:length(ids)) {
    for (j in 1:(i-1)) {
      vars1 <- !is.na(varMatrix[,i+1])
      vars2 <- !is.na(varMatrix[,j+1])
      matches <- sum(vars1 & vars2)
      mismatches <- sum(xor(vars1,vars2))
      jaccard <- matches/sum(vars1 | vars2)
      if (matches >= minMatches && jaccard >= minJaccard) {
        fastclust$acceptMerge(ids[[i]],ids[[j]],
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
  varMatrix1 <- with(fastclust, getClusterOf(id1) |> getMembersOf() |> toVarMatrix())
  varMatrix2 <- with(fastclust, getClusterOf(id1) |> getMembersOf() |> toVarMatrix())
  #extract the genotypes
  genoPair <- genoList[c(id1,id2)]
}

#Iteratively step through connections with greater barcode differences
for (bcDist in 1:maxDiff) {
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



# judgeConnection <- function(id1,id2,bcDist) {
#   #if the pairing was already previously examined, we can skip it
#   if (fastclust$isKnown(id1,id2)) {
#     return()
#   }
#   #extract the genotypes
#   genoPair <- genoList[c(id1,id2)]
#   allVars <- union(genoPair[[1]][,1],genoPair[[2]][,1])
#   if (length(allVars) == 0) {
#     #both are WT
#     fastclust$suggestMerge(
#       id1,id2,
#       bcDist=bcDist,genoDist=0,jaccard=1
#     )
#     return()
#   }
#   #get list of matching variants
#   matches <- intersect(genoPair[[1]][,1],genoPair[[2]][,1])
#   mismatches <- setdiff(allVars,matches)
#   jaccard <- length(matches)/length(allVars)
#   # if (length(mismatches) <= maxDiff && jaccard >= minJaccard) {
#     fastclust$suggestMerge(
#       id1,id2,
#       bcDist=bcDist,genoDist=length(mismatches),jaccard=jaccard
#     )
#   # } else {
#   #   fastclust$rejectMerge(
#   #     id1,id2,
#   #     bcDist=bcDist,genoDist=length(mismatches),jaccard=jaccard
#   #   )
#   # }
#   return()
# }

# for (bcDist in 1:maxDiff) {
#   invisible(apply(edist[edist$dist==bcDist,1:2],1,function(pair) {
#     #the pair of pre-clusters still needs to be resolved into a list of reads
#     ids <- strsplit(pair,"\\|")
#     for (id1 in ids[[1]]) {
#       for (id2 in ids[[2]]) {
#         judgeConnection(id1,id2,bcDist)
#       }
#     }
#   }))
# }



# foo <- lapply(head(preClusts,50),function(ids) {
#   gl <- genoList[ids]
#   out <- Reduce(function(x,y) merge(x,y,by="var",all=TRUE), gl)
#   colnames(out) <- c("var",paste0("qual",1:(ncol(out)-1)))
#   out
# })

# dir.create("export")
# con <- file("export/export.csv",open="w")
# for (i in 1:length(foo)) {
#   write.table(foo[[i]],con,na="",row.names=FALSE)
#   cat("\n\n",file=con)
# }
# close(con)
