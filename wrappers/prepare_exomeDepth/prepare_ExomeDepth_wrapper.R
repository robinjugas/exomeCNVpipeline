suppressMessages(library(data.table))
suppressMessages(library(ExomeDepth))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")


################################################################################
## RUN
run_all <- function(args){
  #read all params as variables
  BED <- args[1]
  FASTA <- args[2]
  Rdata <- args[3]
  BAMS <- args[4:length(args)]
  
  
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
  
  # gene name
  exons$No <- 1
  exons[, rowNo := cumsum(No), by = gene]
  # exons[, name := paste(gene,paste0("E",rowNo),chromosome,start,end,sep = ".")]
  exons[, name := paste0(gene,"_E",rowNo)]
  
  exons[,c("No","rowNo","gene"):=NULL]
  exons <- exons[,c("chromosome", "start", "end", "name")]
  setnames(exons,"name","gene")
  exons <- as.data.frame(exons)
  
  ################################################################################
  ## load COHORT BAMS
  
  BAMcohort <- BAMS
  
  X_cohort <- ExomeDepth::getBamCounts(bed.frame = exons,
                           bam.files = BAMcohort,
                           include.chr = FALSE,
                           referenceFasta = FASTA)
  
  save(X_cohort,exons,file = Rdata)

  
}



# develop and test
# install.packages(c('aod', 'VGAM'))
# install.packages("/home/rj/Downloads/ExomeDepth_1.1.16.tar.gz", repos = NULL, type="source")
# args <- c("/home/rj/4TB/CEITEC/CNV_EXOM_BEDs/OvarianCancer_GRCh38.bed","/home/rj/4TB/CEITEC/homsap/GRCh38/seq/GRCh38.fa","cohort_data/exomeDepth_customCohort.Rdata","mapped/control_1.bam","mapped/control_2.bam","mapped/control_3.bam","mapped/control_4.bam","mapped/control_5.bam","mapped/control_6.bam","mapped/control_7.bam","mapped/control_8.bam")
# setwd("/media/rj/Exos8TB/CNV_OVARIA/FFPE/")

#run as Rscript
#
script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)

