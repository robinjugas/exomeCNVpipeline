suppressMessages(library(data.table))
suppressMessages(library(panelcn.mops))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")


################################################################################
## RUN
run_all <- function(args){
  #read all params as variables
  BAM_file <- args[1]
  lib_ROI <- args[2]
  sampleNameWildCard <- args[3]
  output_file <- args[4]
  
  if(lib_ROI=="HyperExome_GRCh38"){
    Rdata_file <- paste0(script_dir,"/COHORT_exonsHyperExome_panelcn.mops.Rdata")
    load(Rdata_file)
    exons <- exons.GRCh38
  }
  if(lib_ROI=="TruSeq_Exome"){
    Rdata_file <- paste0(script_dir,"/COHORT_exonsTruSeqExome_panelcn.mops.Rdata")
    load(Rdata_file)
  }
  # if(lib_ROI=="XXX"){
  #   load("COHORT_exonsXXXX_cn.MOPS.Rdata")
  #   exons <- exons.GRCh38
  # }
  
  
  
  ################################################################################
  ## load SAMPLE BAMS
  X_sample <- countBamListInGRanges(countWindows = exons,
                                    bam.files = BAM_file, read.width = FALSE)
  
  ################################################################################
  ## PANEL CN.MOPS
  XandCB <- X_sample
  elementMetadata(XandCB) <- cbind(elementMetadata(XandCB),
                                   elementMetadata(X_cohort))
  
  resultlist <- runPanelcnMops(XandCB, testiv = c(1), countWindows = exons)
  
  ################################################################################
  ## Results
  sampleNames <- colnames(elementMetadata(X_sample))
  
  RESULTS <- createResultTable(resultlist = resultlist, XandCB = XandCB,
                               countWindows = exons,
                               sampleNames = sampleNames)
  
  RESULTS <- as.data.frame(RESULTS)
  RESULTS <- RESULTS[RESULTS$CN!="CN2",] #filter out nonCNVs
  setnames(RESULTS,c("Chr","Start","End"),c("chromosome","start","end"))
  RESULTS <- subset(RESULTS,select=c("chromosome","start","end","Gene","Exon","Sample",
                                     "RC","medRC","RC.norm","medRC.norm","lowQual","CN" ))
  
  fwrite(RESULTS,file = output_file,sep = "\t")
  
}



# develop and test WES113
# args <- character(3)
# args[1] <- "/home/rj/4TB/CNV_DATA/WES113/AT-PRO-07krev.bam"
# args[2] <- "HyperExome_GRCh38"
# args[3] <- "AT-PRO-07krev"
# args[4] <- "/home/rj/4TB/CNV_DATA/WES113/AT-PRO-07krev.tsv"
#run as Rscript
# 
script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)