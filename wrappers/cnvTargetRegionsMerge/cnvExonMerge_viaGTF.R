
library(data.table)
library(rtracklayer)

#nemelo by byt napevno... ale s detekci UCSC/Ensembl
listofCHR<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
             "chr19","chr20","chr21","chr22","chrY","chrX")

listofCHR2<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
              "19","20","21","22","X","Y")


run_all <- function(args){
  # arguments
  bed <- args[1]
  tsv <- args[2]
  bed_file <- args[3]
  gtf_file <- args[4]
  tsv_files <- args[5:length(args)]
  
  ## DEFAULT VALUES
  numCallers <- length(tsv_files)
  # if(is.na(distanceThreshold)){distanceThreshold <- 100} #default
  # if(is.na(min_sv_length)){min_sv_length <- 50} #default
  # max_sv_length depends on chromosome length
  
  ## DETECT USCS or ENSEMBL
  input_BED_DF <- fread(file=bed_file, sep='\t', header = F)
  if( grepl( "^chr",     input_BED_DF[1,1], fixed = TRUE) ){
    UCSCorENSEMBL <- "ENS"
  }else{UCSCorENSEMBL <- "UCSC"}
    

  ############################################################################################################################################################
  ## read GTF EXONS OR GENES????
  gtf_GENCODE <- rtracklayer::import(gtf_file)
  gtf_df_all <- as.data.table(gtf_GENCODE)
  gtf_df_all <- gtf_df_all[gtf_df_all$type=="exon",]   #gtf_df_all <- gtf_df_all[gtf_df_all$type=="exon",]
  gtf_df_all <- gtf_df_all[gtf_df_all$gene_type=="protein_coding",]
  # gtf_df_all <- gtf_df_all[gtf_df_all$source=="ENSEMBL",]
  setnames(gtf_df_all,"seqnames","chromosome")
  gtf_df_all <- subset(gtf_df_all,select=c("chromosome","start","end","width","strand","source",
                                           "gene_id","gene_type","gene_name","exon_id","hgnc_id"))
  gtf_df_all <- unique(gtf_df_all, by=c("chromosome","start","end","gene_name"))
  gtf_df_all$chromosome <- as.character(gtf_df_all$chromosome)
  
  if(UCSCorENSEMBL=="UCSC"){
    gtf_df_all$chromosome <- gsub("^chr","",gtf_df_all$chromosome) # carefull
    gtf_df_all <- gtf_df_all[chromosome %in% listofCHR2,]
  }else{
    gtf_df_all <- gtf_df_all[chromosome %in% listofCHR,]
  }
  
  gtf_df_all <- subset(gtf_df_all,select=c("chromosome","start","end","width","strand",
                                           "gene_id","gene_type","gene_name","exon_id","hgnc_id"))
  data.table::setDT(gtf_df_all)
  data.table::setkey(gtf_df_all, chromosome,start,end)
  rm(gtf_GENCODE)
  ############################################################################################################################################################
  ## detect files by caller
  cnvkit <- tsv_files[which(grepl( "/cnvkit", tsv_files, fixed = TRUE))]
  cnmops <- tsv_files[which(grepl( "/cnMOPS", tsv_files, fixed = TRUE))]
  exomedepth <- tsv_files[which(grepl( "/exomeDepth", tsv_files, fixed = TRUE))]
  panelcnmops <- tsv_files[which(grepl( "/panelcnMOPS", tsv_files, fixed = TRUE))]
  
  ##############################################################################
  ## LOAD TSV FILES
  
  #cnvkit
  if(!identical(cnvkit, character(0))){
    cnvkitDF <- fread(file=cnvkit, sep='\t', header = TRUE)
    if(nrow(cnvkitDF)){
      
      cnvkitDF[,gene:=NULL]
      cnvkitDF <- cnvkitDF[cnvkitDF$cn!=2,]
      cnvkitDF <- cnvkitDF[cnvkitDF$chromosome %in% listofCHR2,]
      cnvkitDF$CALLER <- "cnvkit"
      
    }else{cnvkitDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character() )}
  }else{
    cnvkitDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character() )
  }
  
  # cnmops
  if(!identical(cnmops, character(0))){
    cnmopsDF <- fread(file=cnmops, sep='\t', header = TRUE)
    if(nrow(cnmopsDF)){
      
      if("seqnames" %in% colnames(cnmopsDF)){setnames(cnmopsDF,c("seqnames"),c("chromosome"))}
      
      cnmopsDF <- cnmopsDF[cnmopsDF$chromosome %in% listofCHR2,]
      cnmopsDF$CALLER <- "cnmops"
      
      cnmopsDF <- subset(cnmopsDF,select=c("chromosome","start","end","CALLER","CN"))
      cnmopsDF$CN <- as.numeric(gsub("CN","",cnmopsDF$CN))
      
    }else{cnmopsDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character() )}
  }else{
    cnmopsDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character() )
  }
  
  # exomedepth
  if(!identical(exomedepth, character(0))){
    exomedepthDF <- fread(file=exomedepth, sep='\t', header = TRUE)
    if(nrow(exomedepthDF)){
      
      exomedepthDF <- exomedepthDF[exomedepthDF$chromosome %in% listofCHR2,]
      exomedepthDF$CALLER <- "exomedepth"
      
      
      exomedepthDF <- subset(exomedepthDF,select=c("chromosome","start","end","CALLER","type"))
      exomedepthDF$CN <- 0
      exomedepthDF[exomedepthDF$type=="duplication","CN"] <- 10
      exomedepthDF <- subset(exomedepthDF,select=c("chromosome","start","end","CALLER","CN"))
      
      
    }else{exomedepthDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character() )}
  }else{
    exomedepthDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character() )
  }
  
  # panelcnmops
  if(!identical(panelcnmops, character(0))){
    panelcnmopsDF <- fread(file=panelcnmops, sep='\t', header = TRUE)
    if(nrow(panelcnmopsDF)){
      
      panelcnmopsDF <- panelcnmopsDF[panelcnmopsDF$chromosome %in% listofCHR2,]
      panelcnmopsDF$CALLER <- "panelcnmops"
      
      panelcnmopsDF <- subset(panelcnmopsDF,select=c("chromosome","start","end","CALLER","lowQual","CN"))
      panelcnmopsDF <- panelcnmopsDF[panelcnmopsDF$lowQual!="lowQual"]
      panelcnmopsDF <- subset(panelcnmopsDF,select=c("chromosome","start","end","CALLER","CN"))
      panelcnmopsDF$CN <- as.numeric(gsub("CN","",panelcnmopsDF$CN))
      
      
      
    }else{panelcnmopsDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character() )}
  }else{
    panelcnmopsDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character() )
  }
  
  
  ############################################################################################################################################################
  ## MERGE WITH EXONS
  
  # deletions - homo/heterezog
  # duplications
  # separately
  
  cnvkitDF$chromosome <- as.character(cnvkitDF$chromosome)
  data.table::setDT(cnvkitDF)
  data.table::setkey(cnvkitDF, chromosome,start,end)
  
  cnmopsDF$chromosome <- as.character(cnmopsDF$chromosome)
  data.table::setDT(cnmopsDF)
  data.table::setkey(cnmopsDF, chromosome,start,end)
  
  exomedepthDF$chromosome <- as.character(exomedepthDF$chromosome)
  data.table::setDT(exomedepthDF)
  data.table::setkey(exomedepthDF, chromosome,start,end)
  
  panelcnmopsDF$chromosome <- as.character(panelcnmopsDF$chromosome)
  data.table::setDT(panelcnmopsDF)
  data.table::setkey(panelcnmopsDF, chromosome,start,end)
  
  #new ALL CNVs together
  if(nrow(cnmopsDF)){  overlapCNMOPS <- data.table::foverlaps(gtf_df_all,cnmopsDF, by.x = c("chromosome","start","end"), mult="all",type="any",nomatch=NULL, which=FALSE)
  }else{overlapCNMOPS <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"exon_id"=character(),"CN"=numeric() ) }
  if(nrow(panelcnmopsDF)){  overlapPANELCNMOPS <- data.table::foverlaps(gtf_df_all,panelcnmopsDF, by.x = c("chromosome","start","end"), mult="all",type="any",nomatch=NULL, which=FALSE)
  }else{overlapCNMOPS <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"exon_id"=character(),"CN"=numeric() ) }
  if(nrow(exomedepthDF)){  overlapEXOMEDEPTH <- data.table::foverlaps(gtf_df_all,exomedepthDF, by.x = c("chromosome","start","end"), mult="all",type="any",nomatch=NULL, which=FALSE)
  }else{overlapCNMOPS <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"exon_id"=character(),"CN"=numeric() ) }
  
  
  
  #######################################################################################################################
  #udelat PARFOR a projit exony a dat ABSENCE/PRESENCE pokud je exon id pritomno v overlap dataframech... 
  gtf_df_all$panCNMOPS <- NA
  gtf_df_all$CNMOPS <- NA
  gtf_df_all$exomeDepth <- NA
  
  
  all_exons <- c(overlapCNMOPS$exon_id,overlapPANELCNMOPS$exon_id,overlapEXOMEDEPTH$exon_id)
  EXONS <- gtf_df_all[gtf_df_all$exon_id %in% all_exons,]
  
  for(i in 1:nrow(EXONS)){
    exonID <- EXONS$exon_id[i]
    #cnMOPS
    if(length(overlapCNMOPS[overlapCNMOPS$exon_id==exonID,CN]>0)){
      EXONS$CNMOPS[i] <- overlapCNMOPS[overlapCNMOPS$exon_id==exonID,CN]
    }
    # panelCNMOPS
    if(length(overlapPANELCNMOPS[overlapPANELCNMOPS$exon_id==exonID,CN]>0)){
      EXONS$panCNMOPS[i] <- overlapPANELCNMOPS[overlapPANELCNMOPS$exon_id==exonID,CN]
    }
    # exomedepth
    if(length(overlapEXOMEDEPTH[overlapEXOMEDEPTH$exon_id==exonID,CN]>0)){
      EXONS$exomeDepth[i] <- overlapEXOMEDEPTH[overlapEXOMEDEPTH$exon_id==exonID,CN]
    }
    
  }
  
  
  EXONS$Count_Detected <- rowSums(!is.na(EXONS[, 11:13]))
  
  #######################################################################################################################
  # pozor na exomedpth CN numbers
  EXONSfiltered <- EXONS[EXONS$Count_Detected>=2,]
  EXONSfiltered$exomeDepth <- as.character(EXONSfiltered$exomeDepth)
  EXONSfiltered[EXONSfiltered$exomeDepth==0,"exomeDepth"] <- "DEL"
  EXONSfiltered[EXONSfiltered$exomeDepth==10,"exomeDepth"] <- "DUP"
  
  EXONSfiltered <- EXONSfiltered[order(EXONSfiltered$chromosome, EXONSfiltered$start),]
  EXONSfiltered <- unique(EXONSfiltered, by="exon_id")
  
  EXONSfiltered$type <- EXONSfiltered$exomeDepth
  EXONSfiltered[is.na(EXONSfiltered$type) & EXONSfiltered$panCNMOPS<2 & EXONSfiltered$CNMOPS<2  ,"type"] <- "DEL"
  EXONSfiltered[is.na(EXONSfiltered$type) & EXONSfiltered$panCNMOPS>2 & EXONSfiltered$CNMOPS>2  ,"type"] <- "DUP"
  EXONSfiltered <- EXONSfiltered[!is.na(EXONSfiltered$type),]
  
  EXONSbed <- EXONSfiltered[,c("chromosome","start","end","type","gene_name")]
  
  #######################################################################################################################
  dir.create(file.path(getwd(),dirname(bed)), recursive = TRUE, showWarnings = FALSE)
  write.table(EXONSbed, file=bed, sep = "\t", quote = FALSE,row.names = FALSE, col.names = FALSE )
  write.table(EXONSfiltered, file=tsv, sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE )
  
  
  #######################################################################################################################
  
    
    
}


# develop and test

setwd("/media/rj/Exos8TB/CNV_OVARIA/FFPE/")
args <- c("CNV_exon_merged/9466-2022_FFPE.callers_merged.bed","CNV_exon_merged/9466-2022_FFPE.callers_merged.tsv","/home/rj/4TB/CEITEC/CNV_EXOM_BEDs/OvarianCancer_GRCh38.bed","/home/rj/4TB/CEITEC/GENCODE_v47_GRCh38.p14/gencode.v47.primary_assembly.annotation.gtf","variant_calls/9466-2022_FFPE/cnMOPS/cnMOPS_CNV_9466-2022_FFPE.tsv","variant_calls/9466-2022_FFPE/exomeDepth/exomeDepth_CNV_9466-2022_FFPE.tsv","variant_calls/9466-2022_FFPE/panelcnMOPS/panelcnMOPS_CNV_9466-2022_FFPE.tsv")


#run as Rscript

args <- commandArgs(trailingOnly = T)
run_all(args)