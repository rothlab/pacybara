options(stringsAsFactors=FALSE)

library(argparser)
library(yogitools)
library(tileseqMave)
library(pbmcapply)

p <- arg_parser(
	"translate library",
	name="pbVarcall_translate.R"
)
p <- add_argument(p, "varcalls", help="variant call file")
p <- add_argument(p, "parameters", help="parameter json file")
pargs <- parse_args(p)

outfile <- sub("\\.txt$","_transl.csv",pargs$varcalls)
params <- tileseqMave::parseParameters(pargs$parameters)

builder <- hgvsParseR::new.hgvs.builder.p(3)
cbuilder <- hgvsParseR::new.hgvs.builder.c()

vc <- read.delim(pargs$varcalls,header=FALSE)
colnames(vc) <- c("barcode","snvs")
out <- yogitools::as.df(pbmclapply(
	vc[,2],
	function(mut) {
		if (mut=="=") {
			list(hgvsc="c.=",hgvsp="p.=",codonChanges="WT",
				codonHGVS="c.=",aaChanges="WT",aaChangeHGVS="p.=")
		} else if (grepl(";",mut)) {
			if (grepl("ins",mut)) {
				ms <- strsplit(mut,";")[[1]]
				pos <- as.integer(gsub("\\D+","",ms))
				rest <- gsub("\\d+","",ms)
				idx <- grepl("ins",ms)
				ms[idx] <- paste0(pos[idx]-1,"_",pos[idx],rest[idx])
				mut <- paste(ms,collapse=";")
			}
			tileseqMave::translateHGVS(paste0("c.[",mut,"]"),params,builder,cbuilder)
		} else {
			if (grepl("ins",mut)) {
				pos <- as.integer(gsub("\\D+","",mut))
				rest <- gsub("\\d+","",mut)
				mut <- paste0(pos-1,"_",pos,rest)
			}
			tileseqMave::translateHGVS(paste0("c.",mut),params,builder,cbuilder)
		}
	}
	,mc.cores=8
))

out <- cbind(vc,out)

write.csv(out,outfile,row.names=FALSE)
