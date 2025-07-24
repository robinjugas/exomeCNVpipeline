
library(data.table)

#nemelo by byt napevno... ale s detekci UCSC/Ensembl
listofCHR<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
             "chr19","chr20","chr21","chr22","chrY","chrX")

listofCHR2<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
              "19","20","21","22","X","Y")

# develop and test
# setwd("/media/rj/SSD_500GB/CNV_OVARIA/PLAZMY/")
# args <- c("CNV_Whole/942-23_FFPE.callers_merged.tsv","CNV_Whole/942-23_FFPE.classified.txt","CNV_Whole/942-23_FFPE_called_CNVs.tsv","")
# args <- c("CNV_Whole/1456-22_FFPE.callers_merged.tsv","CNV_Whole/1456-22_FFPE.classified.txt","CNV_Whole/1456-22_FFPE_called_CNVs.tsv","")

run_all <- function(args){
  # arguments
  tsv <- args[1]
  txt <- args[2]
  output <- args[3]

  ##############################################################################
  exonsDF <- fread(file=tsv, sep='\t', header = TRUE)
  names(exonsDF)
  exonsDF$VariantID <- paste(exonsDF$CHR,exonsDF$START,exonsDF$STOP,exonsDF$CNVtype,sep = "_") # chr1_970657_970704_DEL
  
  ##############################################################################
  ## DETECT USCS or ENSEMBL
  if( grepl( "^chr",     exonsDF[1,1], fixed = TRUE) ){
    UCSCorENSEMBL <- "ENS"
  }else{UCSCorENSEMBL <- "UCSC"}
  
  ##############################################################################
  classifyCNV_res <- fread(file=txt, sep='\t', header = TRUE)
  if("Chromosome" %in% names(classifyCNV_res)){setnames(classifyCNV_res,c("Chromosome"),c("CHR"))}
  if("Start" %in% names(classifyCNV_res)){setnames(classifyCNV_res,c("Start"),c("START"))}
  if("End" %in% names(classifyCNV_res)){setnames(classifyCNV_res,c("End"),c("STOP"))}
  if("Type" %in% names(classifyCNV_res)){setnames(classifyCNV_res,c("Type"),c("CNVtype"))}
  
 
  
  if(UCSCorENSEMBL=="UCSC"){
    classifyCNV_res$CHR <- gsub("^chr","",classifyCNV_res$CHR) # carefull
    classifyCNV_res <- classifyCNV_res[CHR %in% listofCHR2,]
    classifyCNV_res$VariantID <- gsub("^chr","",classifyCNV_res$VariantID)
  }else{
    classifyCNV_res <- classifyCNV_res[CHR %in% listofCHR,]
  }
  
  classifyCNV_res <- classifyCNV_res[,-c(8:42)] # delete categories
  classifyCNV_res[,c("CHR","START","STOP","CNVtype") :=NULL] # delete cause merge columns
  
  ##############################################################################
  # pozor na chr ve sloupci chromosome a VariantID
  exonsDF <- merge(exonsDF,classifyCNV_res,by="VariantID")
  exonsDF[,c("VariantID") :=NULL] # delete cause merge columns
  
  exonsDF$IGV <- paste0(exonsDF$CHR,":",exonsDF$START,"-",exonsDF$STOP) #chr1:144,874-969,268)
  
  #######################################################################################################################
  dir.create(file.path(getwd(),dirname(output)), recursive = TRUE, showWarnings=FALSE)
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