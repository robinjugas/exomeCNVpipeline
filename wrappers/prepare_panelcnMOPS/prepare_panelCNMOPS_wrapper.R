suppressMessages(library(data.table))
suppressMessages(library(panelcn.mops))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")


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
  exons <- subset(exons, select = -c(width))
  exons <- exons[exons$chromosome!="X",]
  exons <- exons[exons$chromosome!="Y",]
  exons <- unique(exons,by=c("chromosome", "start", "end"))
  exons <- as.data.table(exons)
  
  
  # MUTYH.E16.chr1.45794947.45795140
  exons$No <- 1
  exons[, rowNo := cumsum(No), by = gene]
  exons[, name := paste(gene,paste0("E",rowNo),chromosome,start,end,sep = ".")]
  exons[, exon := paste0("E",rowNo)]
  
  exons[,c("No","rowNo"):=NULL]
  exons <- exons[,c("chromosome", "start", "end", "name", "gene", "exon")]
  exons <- as.data.frame(exons)
  
  
  # ##
  # bed <- system.file("extdata/Genes_part.bed", package = "panelcn.mops")
  # countWindows <- getWindows(bed)
  
  ################################################################################
  ## load COHORT BAMS
  
  BAMcohort <- BAMS
  
  X_cohort <- panelcn.mops::countBamListInGRanges(countWindows = exons,
                                    bam.files = BAMcohort, read.width = FALSE)
  
  save(X_cohort,exons,file=Rdata)
 
  
}



# develop and test
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("panelcn.mops")
# args <- c("/home/rj/4TB/CEITEC/CNV_EXOM_BEDs/OvarianCancer_GRCh38.bed","cohort_data/panelcnMOPS_customCohort.Rdata","mapped/control_1.bam","mapped/control_2.bam","mapped/control_3.bam","mapped/control_4.bam","mapped/control_5.bam","mapped/control_6.bam","mapped/control_7.bam","mapped/control_8.bam")
# setwd("/media/rj/Exos8TB/CNV_OVARIA/FFPE/")
#run as Rscript
# 
script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)

