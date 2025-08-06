
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer))

#nemelo by byt napevno... ale s detekci UCSC/Ensembl
listofCHR<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
             "chr19","chr20","chr21","chr22","chrY","chrX")

listofCHR2<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
              "19","20","21","22","X","Y")


# develop and test
# setwd("/media/rj/Exos8TB/CNV_OVARIA/FFPE/")
# args <- c("CNV_exon_merged/9466-2022_FFPE.callers_merged.bed","CNV_exon_merged/9466-2022_FFPE.callers_merged.tsv","/home/rj/4TB/CEITEC/CNV_EXOM_BEDs/OvarianCancer_GRCh38.bed","/home/rj/4TB/CEITEC/GENCODE_v47_GRCh38.p14/gencode.v47.primary_assembly.annotation.gtf","variant_calls/9466-2022_FFPE/cnMOPS/cnMOPS_CNV_9466-2022_FFPE.tsv","variant_calls/9466-2022_FFPE/exomeDepth/exomeDepth_CNV_9466-2022_FFPE.tsv","variant_calls/9466-2022_FFPE/panelcnMOPS/panelcnMOPS_CNV_9466-2022_FFPE.tsv")

run_all <- function(args){
  # arguments
  output_bed <- args[1]
  output_tsv <- args[2]
  bed_file <- args[3]
  gtf_file <- args[4]
  tsv_files <- args[5:length(args)]
  
  ## DEFAULT VALUES
  numCallers <- length(tsv_files)
  # if(is.na(distanceThreshold)){distanceThreshold <- 100} #default
  # if(is.na(min_sv_length)){min_sv_length <- 50} #default
  # max_sv_length depends on chromosome length
  
  ##############################################################################
  ## DETECT USCS or ENSEMBL
  input_BED_DF <- fread(file=bed_file, sep='\t', header = F)
  if( grepl( "^chr",     input_BED_DF[1,1], fixed = TRUE) ){
    UCSCorENSEMBL <- "ENS"
  }else{UCSCorENSEMBL <- "UCSC"}
  
 
  ##############################################################################
  # PROCESS GTF
  gtf_GENCODE <- rtracklayer::import(gtf_file)
  gtf_df_all <- as.data.table(gtf_GENCODE)
  gtf_df_all <- gtf_df_all[gtf_df_all$type=="gene",]   #gtf_df_all <- gtf_df_all[gtf_df_all$type=="exon",]
  gtf_df_all <- gtf_df_all[gtf_df_all$gene_type=="protein_coding",]
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
  setDT(gtf_df_all)
  setkey(gtf_df_all, chromosome,start,end)
  rm(gtf_GENCODE)

  
  ##############################################################################
  # PROCESS BED
  input_BED_DF <- fread(file=bed_file, sep='\t', header = F)
  setDT(input_BED_DF)
  setnames(input_BED_DF,c("chromosome","BED_start","BED_end","BEDxxx"), skip_absent=TRUE)
  setkey(input_BED_DF, chromosome,BED_start,BED_end)
  
  ############################################################################################################################################################
  ## overlap BED and GTF
  overlapBEDxGTF <- foverlaps(input_BED_DF, gtf_df_all, by.x = c("chromosome","BED_start","BED_end"),
                                          by.y = c("chromosome","start","end"), mult="all", type="any", nomatch=NULL, which=FALSE)
  
  overlapBEDxGTF <- unique(overlapBEDxGTF, by = c("chromosome","BED_start","BED_end","gene_name"))
  setkey(overlapBEDxGTF, chromosome, BED_start,BED_end)
  overlapBEDxGTF[, duplicate := .N > 1, by = key(overlapBEDxGTF)]
  # remove those with ENSG in gene name if there are also duplicate TRUE
  overlapBEDxGTF <- overlapBEDxGTF[!((gene_name %like% "ENSG") & (duplicate==TRUE)),]
  overlapBEDxGTF <- unique(overlapBEDxGTF, by = c("chromosome","BED_start","BED_end","gene_name"))
  
  
  names(overlapBEDxGTF)
  overlapBEDxGTF <- subset(overlapBEDxGTF,select=c("chromosome","BED_start","BED_end","gene_id","gene_type","gene_name","hgnc_id","duplicate"))
  setnames(overlapBEDxGTF,
           c("chromosome","BED_start","BED_end","gene_id","gene_type","gene_name","hgnc_id","duplicate"),
           c("chromosome","BED_start","BED_end","gene_id","gene_type","gene_name","hgnc_id","BEDduplicate"),
           skip_absent=TRUE
           )
  
  overlapBEDxGTF$BED_IGV <- paste0(overlapBEDxGTF$chromosome,":",overlapBEDxGTF$BED_start,"-",overlapBEDxGTF$BED_end) #chr1:144,874-969,268)
  setDT(overlapBEDxGTF)
  overlapBEDxGTF[, gene_name_collapsed := paste(gene_name, collapse = ", "), by = BED_IGV]
  overlapBEDxGTF[, gene_id_collapsed := paste(gene_id, collapse = ", "), by = BED_IGV]
  overlapBEDxGTF[, hgnc_id_collapsed := paste(hgnc_id, collapse = ", "), by = BED_IGV]
  
  overlapBEDxGTF <- unique(overlapBEDxGTF, by = c("chromosome","BED_start","BED_end","gene_name_collapsed"))
  
  overlapBEDxGTF <- subset(overlapBEDxGTF,select=c("chromosome","BED_start","BED_end","gene_name_collapsed","gene_id_collapsed","hgnc_id_collapsed"))
  
  setnames(overlapBEDxGTF,
           c("chromosome","BED_start","BED_end","gene_name_collapsed","gene_id_collapsed","hgnc_id_collapsed"),
           c("chromosome","BED_start","BED_end","gene_name","gene_id","hgnc_id"),
           skip_absent=TRUE
  )
  
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

      if(UCSCorENSEMBL=="UCSC"){
        cnvkitDF$chromosome <- gsub("^chr","",cnvkitDF$chromosome) # carefull
        cnvkitDF <- cnvkitDF[chromosome %in% listofCHR2,]
      }else{
        cnvkitDF <- cnvkitDF[chromosome %in% listofCHR,]
      }
      
      cnvkitDF$CALLER <- "cnvkit"
      
    }else{cnvkitDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character(),"CN"=numeric(),"TYPE"=character() )}
  }else{
    cnvkitDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character(),"CN"=numeric(),"TYPE"=character() )
  }
  
  # cnmops
  if(!identical(cnmops, character(0))){
    cnmopsDF <- fread(file=cnmops, sep='\t', header = TRUE)
    if(nrow(cnmopsDF)){
      if("seqnames" %in% colnames(cnmopsDF)){setnames(cnmopsDF,c("seqnames"),c("chromosome"))}
      #handle USCS or ENSEMBL
      if(UCSCorENSEMBL=="UCSC"){
        cnmopsDF$chromosome <- gsub("^chr","",cnmopsDF$chromosome) # carefull
        cnmopsDF <- cnmopsDF[chromosome %in% listofCHR2,]
      }else{
        cnmopsDF <- cnmopsDF[chromosome %in% listofCHR,]
      }
      # HANDLE CN AND CNV TYPE
      cnmopsDF$CN <- as.numeric(gsub("CN","",cnmopsDF$CN))
      cnmopsDF$TYPE <- "NA"
      cnmopsDF[cnmopsDF$CN<2,"TYPE"] <- "DEL"
      cnmopsDF[cnmopsDF$CN>2,"TYPE"] <- "DUP"
      
      cnmopsDF$CALLER <- "cnmops"
      
      cnmopsDF <- subset(cnmopsDF,select=c("chromosome","start","end","CALLER","CN","TYPE"))
      
    }else{cnmopsDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character(),"CN"=numeric(),"TYPE"=character() )}
  }else{
    cnmopsDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character(),"CN"=numeric(),"TYPE"=character() )
  }
  
  # exomedepth
  if(!identical(exomedepth, character(0))){
    exomedepthDF <- fread(file=exomedepth, sep='\t', header = TRUE)
    if(nrow(exomedepthDF)){
      if("seqnames" %in% colnames(exomedepthDF)){setnames(exomedepthDF,c("seqnames"),c("chromosome"))}
      #handle USCS or ENSEMBL
      if(UCSCorENSEMBL=="UCSC"){
        exomedepthDF$chromosome <- gsub("^chr","",exomedepthDF$chromosome) # carefull
        exomedepthDF <- exomedepthDF[chromosome %in% listofCHR2,]
      }else{
        exomedepthDF <- exomedepthDF[chromosome %in% listofCHR,]
      }
      # HANDLE CN AND CNV TYPE
      exomedepthDF$CN <- NA
      exomedepthDF$TYPE <- "NA"
      exomedepthDF[exomedepthDF$type=="deletion","TYPE"] <- "DEL"
      exomedepthDF[exomedepthDF$type=="duplication","TYPE"] <- "DUP"
      
      exomedepthDF$CALLER <- "exomedepth"

      exomedepthDF <- subset(exomedepthDF,select=c("chromosome","start","end","CALLER","CN","TYPE"))
      
    }else{exomedepthDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character(),"CN"=numeric(),"TYPE"=character() )}
  }else{
    exomedepthDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character(),"CN"=numeric(),"TYPE"=character() )
  }
  
  # panelcnmops
  if(!identical(panelcnmops, character(0))){
    panelcnmopsDF <- fread(file=panelcnmops, sep='\t', header = TRUE)
    if(nrow(panelcnmopsDF)){
      if("seqnames" %in% colnames(panelcnmopsDF)){setnames(panelcnmopsDF,c("seqnames"),c("chromosome"))}
      #handle USCS or ENSEMBL
      if(UCSCorENSEMBL=="UCSC"){
        panelcnmopsDF$chromosome <- gsub("^chr","",panelcnmopsDF$chromosome) # carefull
        panelcnmopsDF <- panelcnmopsDF[chromosome %in% listofCHR2,]
      }else{
        panelcnmopsDF <- panelcnmopsDF[chromosome %in% listofCHR,]
      }
      # delete lowQual
      panelcnmopsDF <- panelcnmopsDF[panelcnmopsDF$lowQual!="lowQual",]
      # HANDLE CN AND CNV TYPE
      panelcnmopsDF$CN <- as.numeric(gsub("CN","",panelcnmopsDF$CN))
      panelcnmopsDF$TYPE <- "NA"
      panelcnmopsDF[panelcnmopsDF$CN<2,"TYPE"] <- "DEL"
      panelcnmopsDF[panelcnmopsDF$CN>2,"TYPE"] <- "DUP"
      
      panelcnmopsDF$CALLER <- "panelcnmops"

      panelcnmopsDF <- subset(panelcnmopsDF,select=c("chromosome","start","end","CALLER","CN","TYPE"))

    }else{panelcnmopsDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character(),"CN"=numeric(),"TYPE"=character() )}
  }else{
    panelcnmopsDF <- data.table("chromosome"=character(),"start"=numeric(),"end"=numeric(),"CALLER"=character(),"CN"=numeric(),"TYPE"=character() )
  }
  
  
  ############################################################################################################################################################
  ## MERGE WITH EXONS
  
  # deletions - homo/heterezog
  # duplications
  # separately
  
  cnvkitDF$chromosome <- as.character(cnvkitDF$chromosome)
  setDT(cnvkitDF)
  setkey(cnvkitDF, chromosome,start,end)
  setnames(cnvkitDF,c("TYPE","CN"),c("TYPE_cnvkit","CN_cnvkit"))
  
  cnmopsDF$chromosome <- as.character(cnmopsDF$chromosome)
  setDT(cnmopsDF)
  setkey(cnmopsDF, chromosome,start,end)
  setnames(cnmopsDF,c("TYPE","CN"),c("TYPE_cnmops","CN_cnmops"))
  
  exomedepthDF$chromosome <- as.character(exomedepthDF$chromosome)
  setDT(exomedepthDF)
  setkey(exomedepthDF, chromosome,start,end)
  setnames(exomedepthDF,c("TYPE","CN"),c("TYPE_exomedepthDF","CN_exomedepthDF"))
  
  panelcnmopsDF$chromosome <- as.character(panelcnmopsDF$chromosome)
  setDT(panelcnmopsDF)
  setkey(panelcnmopsDF, chromosome,start,end)
  setnames(panelcnmopsDF,c("TYPE","CN"),c("TYPE_pcnmops","CN_pcnmops"))
  
  ############################################################################################################################################################
  
  CNVS <- rbind(cnvkitDF,cnmopsDF,exomedepthDF,panelcnmopsDF, fill=TRUE )
  CNVS$IGV <- paste0(CNVS$chromosome,":",CNVS$start,"-",CNVS$end) #chr1:144,874-969,268)
  CNVS$IGV_CALLER <- paste0(CNVS$chromosome,":",CNVS$start,"-",CNVS$end,"_",CNVS$CALLER) #chr1:144,874-969,268)
  overlapBEDxGTF$BED_IGV <- paste0(overlapBEDxGTF$chromosome,":",overlapBEDxGTF$BED_start,"-",overlapBEDxGTF$BED_end) #chr1:144,874-969,268)
  
  setDT(CNVS)
  setkey(CNVS, chromosome,start,end)
  
  setDT(overlapBEDxGTF)
  setkey(overlapBEDxGTF, chromosome,BED_start,BED_end)
  
  BED_CNVs <- foverlaps(overlapBEDxGTF, CNVS, by.x = c("chromosome","BED_start","BED_end"),
                                          by.y = c("chromosome","start","end"), mult="all", type="any", nomatch=NULL, which=FALSE)
  
  names(BED_CNVs)
  setnames(BED_CNVs, c("start","end"),c("CNV_start","CNV_end"),skip_absent=TRUE)
 
  
  # setkey(BED_CNVs, chromosome, BED_start, BED_end, BED_IGV)
  # BED_CNVs[, CNVduplicate := .N > 1, by = key(BED_CNVs)]
  
  #######################################################################################################################
  # 
  setkey(BED_CNVs, BED_IGV, IGV_CALLER)
  BED_CNVs[, collapsed_CNVs := paste(IGV_CALLER, collapse = ", "), by = BED_IGV]
  BED_CNVs[, collapsed_CALLERS := paste(CALLER, collapse = ", "), by = BED_IGV]
  BED_CNVs[, CallersNo := .N, by = list(BED_IGV)]
  
  names(BED_CNVs)

  
  
  # Define mode function
  get_mode <- function(x) {
    x <- na.omit(x)
    ux <- unique(x)
    if (length(ux) == 0) return(NA)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # Apply mode row-wise to selected columns
  BED_CNVs[, CN := apply(.SD, 1, get_mode), .SDcols = c("CN_cnvkit","CN_cnmops","CN_pcnmops","CN_exomedepthDF")]
  BED_CNVs[, TYPE := apply(.SD, 1, get_mode), .SDcols = c("TYPE_cnvkit","TYPE_cnmops","TYPE_exomedepthDF","TYPE_pcnmops")]
  
  
  setcolorder(BED_CNVs, c("chromosome","BED_start","BED_end","BED_IGV","collapsed_CNVs","collapsed_CALLERS","CallersNo",
                          "CN","TYPE",
                          "IGV","IGV_CALLER","gene_name","gene_id","hgnc_id",
                          "CNV_start","CNV_end", "CALLER",
                          "CN_cnvkit","CN_cnmops","CN_pcnmops","CN_exomedepthDF",
                          "TYPE_cnvkit","TYPE_cnmops","TYPE_exomedepthDF","TYPE_pcnmops"
                          ), skip_absent=TRUE)
    
  BED_CNVs <- BED_CNVs[, (c("IGV_CALLER","CNV_start","CNV_end","CALLER")) := NULL]
  
  BED_CNVs_BED <- BED_CNVs[,c("chromosome","BED_start","BED_end","TYPE","gene_name","CN","BED_IGV")]
  
  #######################################################################################################################
  dir.create(file.path(getwd(),dirname(output_bed)), recursive = TRUE, showWarnings = FALSE)
  write.table(BED_CNVs_BED, file=output_bed, sep = "\t", quote = FALSE,row.names = FALSE, col.names = FALSE )
  write.table(BED_CNVs, file=output_tsv, sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE )
  
  
  #######################################################################################################################
  
    
    
}


#run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)