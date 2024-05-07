suppressMessages(library(data.table))
suppressMessages(library(cn.mops))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")


################################################################################
## MERGE - works
# Add extra columns to GRanges Metadata
# https://support.bioconductor.org/p/43654/ WOW!!!!

setMethod("$", "GRanges", function(x, name) { # {{{
  elementMetadata(x)[, name]
}) # }}}

setMethod("$<-", "GRanges", function(x, name, value) { # {{{
  elementMetadata(x)[[ name ]] <- value
  return(x)
}) # }}}

################################################################################
## RUN
run_all <- function(args){
  #read all params as variables
  BED <- args[1]
  Rdata <- args[2]
  BAMS <- args[3:length(args)]
  
  ################################################################################
  ## load exons
  exons <- fread(BED,sep="\t")
  colnames(exons) <- c("chromosome","start","end","gene")
  exons <- exons[exons$end>exons$start,] # end > start
  exons$width <- exons$end-exons$start
  exons <- exons[exons$width>10,] #or get segmentation error https://support.bioconductor.org/p/9145581/
  exons <- subset(exons, select = -c(width,gene))
  exons <- exons[exons$chromosome!="X",]
  exons <- exons[exons$chromosome!="Y",]
  exons <- unique(exons,by=c("chromosome", "start", "end"))
  exons <- as.data.frame(exons)
  # exons <- GRanges(exons[,1],IRanges(exons[,2],exons[,3]))
  
  # Results can also be improved if you extend your target regions by a small amount of bases to # the left and to the right (in the following case it is 30bp)
  exons <- GRanges(exons[,1],IRanges(exons[,2]-30,exons[,3]+30))
  exons <- reduce(exons)
  
  ################################################################################
  ## load COHORT BAMS
  BAMcohort <- BAMS
  
  X_cohort <- getSegmentReadCountsFromBAM(BAMcohort, GR=exons, parallel=8)
  
  save(X_cohort,exons,file=Rdata)

  
  
  
  
}



# develop and test WES113
# args <- character(3)
# args[1] <- "/home/rj/4TB/CNV_DATA/WES113/AT-PRO-07krev.bam"
# args[2] <- "HyperExome_GRCh38"
# args[3] <- "AT-PRO-07krev"
# args[4] <- "/home/rj/4TB/CNV_DATA/WES113/AT-PRO-07krev.tsv"
#run as Rscript

script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)