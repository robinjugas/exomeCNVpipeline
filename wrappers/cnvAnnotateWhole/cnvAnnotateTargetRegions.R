
library(data.table)

#nemelo by byt napevno... ale s detekci UCSC/Ensembl
listofCHR<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
             "chr19","chr20","chr21","chr22","chrY","chrX")

listofCHR2<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
              "19","20","21","22","X","Y")
# develop and test
setwd("/media/rj/Exos8TB/CNV_OVARIA/FFPE/")
args <- c("CNV_TargetRegions/9466-2022_FFPE.callers_merged.tsv","CNV_TargetRegions/9466-2022_FFPE.classified.txt","CNV_TargetRegions/9466-2022_FFPE_called_TargetRegions.tsv","")



run_all <- function(args){
  # arguments
  tsv <- args[1]
  txt <- args[2]
  output <- args[3]


  
  classifyCNV_res <- fread(file=txt, sep='\t', header = TRUE)
  classifyCNV_res <- classifyCNV_res[,c(1,6)]
  classifyCNV_res$VariantID <- gsub("^chr","",classifyCNV_res$VariantID)
  
  exonsDF <- fread(file=tsv, sep='\t', header = TRUE)
  exonsDF$VariantID <- paste(exonsDF$chromosome,exonsDF$BED_start,exonsDF$BED_end,exonsDF$TYPE,sep = "_") # chr1_970657_970704_DEL
  
  
  # pozor na chr ve sloupci chromosome a VariantID
  exonsDF <- merge(exonsDF,classifyCNV_res,by="VariantID")
  exonsDF <- subset(exonsDF, select = -c(VariantID))
  #######################################################################################################################
  dir.create(file.path(getwd(),dirname(output)), recursive = TRUE)
  write.table(exonsDF, file=output, sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE )
  
}


# setwd("/home/rj/4TB/CNV_DATA_WES/resultsHANKA_TRUSEQEXOME/")
# 
# args <- c("CNV_exon_merged/DB9905tumor2017.merged.bed",
#           "CNV_exon_merged/DB9905tumor2017.merged.tsv",
#           "/home/rj/4TB/CEITEC/GTFs_GRCh37/gencode.v41lift37.basic.annotation.gtf",
#           "variant_calls/DB9905tumor2017/cnMOPS/cnMOPS_CNV_DB9905tumor2017.tsv",
#           "variant_calls/DB9905tumor2017/exomeDepth/exomeDepth_CNV_DB9905tumor2017.tsv",
#           "variant_calls/DB9905tumor2017/panelcnMOPS/panelcnMOPS_CNV_DB9905tumor2017.tsv")



#run as Rscript

args <- commandArgs(trailingOnly = T)
run_all(args)