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
  BAM_file <- args[1]
  lib_ROI <- args[2]
  sampleNameWildCard <- args[3]
  output_file <- args[4]
  
  if(lib_ROI=="HyperExome_GRCh38"){
    Rdata_file <- paste0(script_dir,"/COHORT_exonsHyperExome_cn.MOPS.Rdata")
    load(Rdata_file)
    exons <- exons.GRCh38
    }
  if(lib_ROI=="TruSeq_Exome"){
    Rdata_file <- paste0(script_dir,"/COHORT_exonsTruSeqExome_cn.MOPS.Rdata")
    load(Rdata_file)
  }
  
  if(lib_ROI=="HyperExome_w_CNVbb_GRCh38"){
    Rdata_file <- paste0(script_dir,"/COHORT_exonsHyperExome_w_CNVbb_cn.MOPS.Rdata")
    load(Rdata_file)
    exons <- exons.GRCh38
  }
  # if(lib_ROI=="XXX"){
  #   load("COHORT_exonsXXXX_cn.MOPS.Rdata")
  #   exons <- exons.GRCh38
  # }
  
  
  
  ################################################################################
  ## load SAMPLE BAMS
  X_sample <- getSegmentReadCountsFromBAM(BAM_file,sampleNames="SAMPLENAME", GR=exons)
  X_sampleDF <- data.frame(X_sample)
  
  X_cohort$SAMPLENAME <- X_sampleDF$SAMPLENAME
  
  ################################################################################
  ## CN.MOPS
  resCNMOPS <- exomecn.mops(X_cohort)
  resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
  
  ################################################################################
  # results
  
  # segm <- as.data.frame(segmentation(resCNMOPS))
  # CNVRegions <- as.data.frame(cnvr(resCNMOPS))
  # segm <- segm[segm$sampleName=="SAMPLENAME",]
  # segm$sampleName <- sampleNameWildCard
  
  CNVs <- as.data.frame(cnvs(resCNMOPS))
  
  CNVs <- CNVs[CNVs$sampleName=="SAMPLENAME",]
  
  if(nrow(CNVs)>0){
    CNVs$sampleName <- sampleNameWildCard
    setnames(CNVs,"seqnames","chromosome")
    dir.create(file.path(dirname(output_file)), showWarnings = FALSE)
    fwrite(CNVs,file = output_file,sep = "\t")
    
  }else{
    CNVs <- data.frame(seqnames=character(0),start=character(0),end=character(0),width=character(0),strand=character(0),sampleName=character(0),
                       median=character(0),mean=character(0),CN=character(0))
    dir.create(file.path(dirname(output_file)), showWarnings = FALSE)
    fwrite(CNVs,file = output_file,sep = "\t")
    
  }
  
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