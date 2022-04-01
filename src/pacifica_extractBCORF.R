#!/usr/bin/env Rscript
library(bitops)
library(yogitools)
library(yogiseq)
library(hgvsParseR)

options(stringsAsFactors=FALSE)

#wrapper function to extract variant base calls en-masse
extractBasecalls <- function(samStream,passfilter) {
  lapply(which(passfilter),function(i) yogiseq::varcallCandidates(
    samStream$getFlagRow(i),
    samStream$getMDZ(i),
    samStream$getCigar(i),
    as.integer(samStream$getSamElement(i,"pos")),
    samStream$getSamElement(i,"seq"),
    samStream$getSamElement(i,"qual")
  ))
}

#' separate barcode insertions from actual variants
#' 
#' also repairs mis-identified barcode starts
#' 
#' @param bcs basecall table from extractBasecalls()
#' @param bcLen expected barcode length 
#' @param bcPoss vector of expected barcode start positions
#' @param refSeq reference sequence (to repair shifted barcodes)
#' @return list containing two tables: barcodes and variants
separateBCs <- function(bcs,bcLen,bcPoss,refSeq) {

  if (nrow(bcs)==0) {
    return(NULL)
  }
  insIdx <- which(bcs$op=="ins")
  insCalls <- bcs[insIdx,]
  if (nrow(insCalls) == 0) {
    return(NULL)
  }

  #expected barcode positions after being shifted due to 
  #exclusion from reference sequence
  pShift <- bcPoss-(seq_along(bcPoss)-1)*bcLen

  #lengths
  ls <- nchar(insCalls$readbase)
  #positions
  ps <- insCalls$refpos+1
  #find barcode ID candidates
  bcIDs <- sapply(ps, function(p) {
    dists <- abs(p-pShift)
    out <- which(dists < 3)
    if (length(out) != 1) NA else out
  })
  #remove candidates that are too short
  if (any(ls < (bcLen-3))) {
    bcIDs[which(ls < (bcLen-3))] <- NA
  }
  bcIdx <- insIdx[which(!is.na(bcIDs))]

  #separate basecall lines into barcodes and genotype details
  barcodeLines <- cbind(bcs[bcIdx,],bcId=na.omit(bcIDs),bcLen=ls[!is.na(bcIDs)])
  genoLines <- bcs[-bcIdx,]

  #repair barcodes that were shifted
  if (nrow(barcodeLines) > 0) {
    for (j in 1:nrow(barcodeLines)) {
      # bcLen <- barcodeLines$bcLen[[j]]
      expected <- pShift[barcodeLines$bcId[[j]]]
      found <- barcodeLines$refpos[[j]]
      shift <- found-expected
      if (shift < 0) {
        bcSeq <- barcodeLines$readbase[[j]]
        prefix <- substr(bcSeq,1,-shift)
        refEquiv <- substr(refSeq,expected+shift,expected-1)
        if (prefix == refEquiv) {
          #it would be better to have the full read sequence and qual string as an argument 
          # and re-extract as appropriate
          barcodeLines$readbase[[j]] <- paste0(substr(bcSeq,1-shift,nchar(bcSeq)),prefix)
          barcodeLines$refpos[[j]] <- expected
          qstr <- barcodeLines$qual[[j]]
          #especially this is technically not quite right.
          barcodeLines$qual[[j]] <- paste0(substr(qstr,1-shift,nchar(qstr)),substr(qstr,1,-shift))
        }
      } else if (shift > 0) {
        bcSeq <- barcodeLines$readbase[[j]]
        suffix <- substr(bcSeq,nchar(bcSeq)-shift+1,nchar(bcSeq))
        refEquiv <- substr(refSeq,expected+1,expected+shift)
        if (suffix == refEquiv) {
          barcodeLines$readbase[[j]] <- paste0(suffix,substr(bcSeq,1,nchar(bcSeq)-shift))
          barcodeLines$refpos[[j]] <- expected
          qstr <- barcodeLines$qual[[j]]
          qsuff <- substr(qstr,nchar(qstr)-shift+1,nchar(qstr))
          barcodeLines$qual[[j]] <- paste0(qsuff,substr(qstr,1,nchar(qstr)-shift))
        }
      }
    }
  }

  colnames(barcodeLines) <- c("op","ref","pos","seq","qual","bcID","bcLen")
  barcodeLines <- barcodeLines[order(barcodeLines$bcID),]

  list(barcodes=barcodeLines[,-(1:2)],genotype=genoLines)

}

fetchBC <- function(bcGenos,bcID) {
  do.call(rbind,lapply(bcGenos, function(bg) {
    if (bcID %in% bg$barcodes$bcID) {
      idx <- which(bg$barcodes$bcID==bcID)
      bg$barcodes[idx[[1]],-4]
    } else {
      data.frame(pos=NA_integer_,seq="",qual="",bcLen=0L)
    }
  }))
}

combineBCs <- function(bcGenos) {
  as.df(lapply(bcGenos, function(bg) {
    if (!is.null(bg$barcodes) && length(bg$barcodes) > 0 && nrow(bg$barcodes) > 0) {
      with(bg$barcodes,list(
        pos=NA_integer_,
        seq=paste(seq,collapse=""),
        qual=paste(qual,collapse=""),
        bcLen=sum(bcLen)
      ))
    } else {
      list(pos=NA_integer_,seq="",qual="",bcLen=0L)
    }
  }))
}

toFASTQ <- function(bcs) {
  headers <- with(bcs,sprintf("@%s pos=%d len=%d",read,pos,bcLen))
  with(bcs,paste(paste(headers,seq,"+",qual,sep="\n"),collapse="\n"))
}

#wrapper function to extract variant base calls en-masse
extractMatchRanges <- function(samStream) {
  do.call(rbind,lapply(1:length(samStream$getSamElement(,1)),function(i) {
    start <- as.integer(samStream$getSamElement(i,"pos"))
    # cigar <- samStream$getCigar(i)
    if (is.na(start)) {
      return(c(NA,NA))
    }
    mdz <- samStream$getMDZ(i)
    if (any(is.na(mdz))) {
      return(c(NA,NA))
    }
    bpCovered <- sum(sapply(mdz, function(elem) {
      if (is.numeric(elem)) return(elem)
      else if (substr(elem,1,1)=="^") return(nchar(elem)-1)
      else return(nchar(elem))
    }))
    return(c(start,start+bpCovered-1))
  }))
}


processSAMs <- function(sam.file,refFasta,outdir,chunkSize=100,bcLen=25,
      bcPoss=c(153,2812),orfStart=207,orfEnd=2789) {

  #adjust ORF positions relative to barcode cutouts
  if (any(bcPoss < orfStart)) {
    orfStart <- orfStart - bcLen*sum(bcPoss < orfStart)
    orfEnd <- orfEnd - bcLen*sum(bcPoss < orfStart)
  }
  #calculate the minimum range of positions that needs to be covered
  bcAdj <- bcPoss-(seq_along(bcPoss)-1)*bcLen
  minRange <- range(c(orfStart,orfEnd,bcAdj,bcAdj+bcLen))

  refCon <- file(refFasta,open="r")
  refSeq <- yogiseq::readFASTA(refCon)[[1]]$toString()
  close(refCon)

  #open stream
  stream <- yogiseq::new.sam.streamer(sam.file,chunkSize)

  #open output streams
  bcOuts <- lapply(1:length(bcPoss),function(i) {
    gzfile(sprintf("%s/bcExtract_%d.fastq.gz",outdir,i),open="w")
  })
  bcComboOut <- gzfile(sprintf("%s/bcExtract_combo.fastq.gz",outdir),open="w")
  genoOut <- gzfile(sprintf("%s/genoExtract.csv.gz",outdir),open="w")

  #error counter
  exceptions <- yogitools::new.counter()

  #tracker for the number of lines processed
  linesDone <- 0

  cat("Processing SAM file...\n\n")

  #main loop for processing stream
  done <- FALSE
  while(!done) {

    #load next chunk in the stream
    nlines <- stream$nextChunk()
    if (nlines==0) {
      break
    }
    #check if we've reached the end of the stream
    if (nlines < chunkSize) {
      done <- TRUE
    }

    #vector of line numbers
    # lineNums <- linesDone + (1:nlines)
    flags <- stream$getFlags()
    passfilter <- with(flags, !segmentUnmapped & !revComp & !failQC & !secondary)
    exceptions$add("failQC",sum(!passfilter))

    #check whether reads cover the entire range of barcodes and ORF
    matchRanges <- extractMatchRanges(stream)
    inRange <- matchRanges[,1] <= minRange[[1]] & matchRanges[,2] >= minRange[[2]]
    inRange[is.na(inRange)] <- FALSE
    exceptions$add("insufficientCoverage",sum(!inRange & passfilter))
    passfilter <- passfilter & inRange

    #if nothing passes filter, skip to next chunk
    if (!any(passfilter)) {
      linesDone <- linesDone + nlines
      cat("\r",linesDone,"lines processed")
      next
    }

    rIDs <-  stream$getSamElement(,"cname")
    rIDs <- rIDs[passfilter]

    #extract variant call candidates
    #barcodes will appear as long insertions
    baseCalls <- extractBasecalls(stream,passfilter)

    #list of corrected barcode-genotype associations
    bcGenos <- lapply(1:length(baseCalls), function(i) {
      separateBCs(baseCalls[[i]],bcLen,bcPoss,refSeq)
    })

    #check if barcodes were found
    numBCFound <- sapply(bcGenos,function(x) {
      if (is.null(x) || !("barcodes" %in% names(x))) return(0) 
      else nrow(x$barcodes)
    })
    missingBC <- numBCFound < length(bcPoss)
    #remove entries with missing barcodes
    if (any(missingBC)) {
      exceptions$add("missingBarcode",sum(missingBC))
      rIDs <- rIDs[which(!missingBC)]
      bcGenos <- bcGenos[which(!missingBC)]
    }
    #if no barcodes were found, skip to next chunk
    if (all(missingBC)) {
      linesDone <- linesDone + nlines
      cat("\r",linesDone,"lines processed")
      next
    }

    #write barcodes to output streams
    for (bcID in 1:length(bcPoss)) {
      bcs <- cbind(read=rIDs,fetchBC(bcGenos,bcID))
      cat(toFASTQ(bcs),"\n",file=bcOuts[[bcID]],sep="")
    }
    comboBCs <- cbind(read=rIDs,combineBCs(bcGenos))
    cat(toFASTQ(comboBCs),"\n",file=bcComboOut,sep="")

    #construct genotype descriptor strings (based on HGVS syntax)
    genoStrs <- sapply(bcGenos,function(bg) {
      if (!is.null(bg$genotype) && nrow(bg$genotype) > 0) {
        qualNum <- sapply(bg$genotype$qual,function(q) 
          as.integer(round(mean(as.integer(charToRaw(q)))))
        )
        muts <- sapply(1:nrow(bg$genotype), function(i) with(bg$genotype[i,],{
          #convert reference position to relative position (wrt ORF)
          relpos <- refpos-orfStart+1
          switch(op,
            sub=sprintf("%d%s>%s:%d",relpos,refbase,readbase,qualNum[[i]]),
            ins=sprintf("%dins%s:%d",relpos,readbase,qualNum[[i]]),
            del=sprintf("%ddel:%d",relpos,qualNum[[i]])
          )
        }))
        paste(muts,collapse=";")
      } else "="
    })
    #write genotypes to output file
    write.table(
      data.frame(read=rIDs,geno=genoStrs), genoOut, 
      sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE
    )
    
    linesDone <- linesDone + nlines

    cat("\r",linesDone,"lines processed")

  }

  cat("\n")

  cat(exceptions$export(),"\n")

  lapply(bcOuts,close)
  close(bcComboOut)
  close(genoOut)

}

library(argparser)

# pargs <- list(
#   sam="LPLM13_head.sam",
#   ref="LPL_M13.fa",
#   outdir="./",
#   chunkSize=100L,
#   bcLen=25L,
#   bcPos="153,1762",
#   orfStart=207L,
#   orfEnd=1631L
# )
p <- arg_parser(
  "Extract barcodes and genotypes from alignments",
  name="extractBCORF.R"
)
p <- add_argument(p, "sam", help="sam file (no header)")
p <- add_argument(p, "ref", help="reference FASTA file (barcodes excised)")
p <- add_argument(p, "outdir", help="output directory")
p <- add_argument(p, "--chunkSize", help="number of lines to process per iteration.",default=100L)
p <- add_argument(p, "--bcLen", help="barcode length",default=25L)
p <- add_argument(p, "--bcPos", help="comma-separated list of barcode positions",default="153,2812")
p <- add_argument(p, "--orfStart", help="ORF start position",default=207L)
p <- add_argument(p, "--orfEnd", help="ORF end position",default=2789L)
pargs <- parse_args(p)


bcPoss <- as.integer(strsplit(pargs$bcPos,",")[[1]])

stopifnot(exprs={
  yogitools::canRead(pargs$sam)
  yogitools::canRead(pargs$ref)
  dir.exists(pargs$outdir)
  is.integer(pargs$chunkSize)
  !any(is.na(bcPoss))
  is.integer(pargs$bcLen)
  is.integer(pargs$orfStart)
  is.integer(pargs$orfEnd)
  with(pargs,(orfEnd-orfStart+1)%%3==0)
})

processSAMs(pargs$sam,pargs$ref,pargs$outdir,
  chunkSize=pargs$chunkSize,bcLen=pargs$bcLen,bcPoss=bcPoss,
  orfStart=pargs$orfStart,orfEnd=pargs$orfEnd
)

cat("\nDone!\n")
