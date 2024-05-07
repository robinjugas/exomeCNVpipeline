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
  
  X_cohort <- getBamCounts(bed.frame = exons,
                           bam.files = BAMcohort,
                           include.chr = FALSE,
                           referenceFasta = FASTA)
  
  save(X_cohort,exons,file = Rdata)

  
}



# develop and test WES113
# args <- character(3)
# args[1] <- "/home/rj/4TB/CNV_DATA/WES113/AT-PRO-07krev.bam"
# args[2] <- "HyperExome_GRCh38"
# args[3] <- "AT-PRO-07krev"
# args[4] <- "/home/rj/4TB/CEITEC/homsap/GRCh38-p10/seq/GRCh38-p10.fa"
# args[5] <- "/home/rj/4TB/CNV_DATA/WES113/AT-PRO-07krev.tsv"
#run as Rscript
#
script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)