
suppressPackageStartupMessages(library(data.table))

#nemelo by byt napevno... ale s detekci UCSC/Ensembl
listofCHR<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
             "chr19","chr20","chr21","chr22","chrY","chrX")

listofCHR2<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
              "19","20","21","22","X","Y")

# develop and test
setwd("/media/rj/Exos8TB/CNV_OVARIA/FFPE/")
args <- c("CNV_Whole/1922-22_FFPE.callers_merged.bed","CNV_Whole/1922-22_FFPE.callers_merged.tsv","/home/rj/4TB/CEITEC/CNV_EXOM_BEDs/OvarianCancer_GRCh38.bed","variant_calls/1922-22_FFPE/cnMOPS/cnMOPS_CNV_1922-22_FFPE.tsv","variant_calls/1922-22_FFPE/exomeDepth/exomeDepth_CNV_1922-22_FFPE.tsv","variant_calls/1922-22_FFPE/panelcnMOPS/panelcnMOPS_CNV_1922-22_FFPE.tsv")


run_all <- function(args){
  # arguments
  output_bed <- args[1]
  output_tsv <- args[2]
  bed_file <- args[3]
  tsv_files <- args[4:length(args)]
  
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
  
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  ## MERGE WITH EACH OTHER
  
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
  
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  
  OVERLAPS <- rbind(cnvkitDF,cnmopsDF,exomedepthDF,panelcnmopsDF,fill=TRUE)
  setDT(OVERLAPS)
  setnames(OVERLAPS,c("chromosome", "start", "end"),c("CHR", "x.start", "x.stop"))
  
  OVERLAPS <- subset(OVERLAPS,select=c("CHR", "x.start", "x.stop"))
  
  OVERLAPS <- unique(OVERLAPS, by=c("CHR", "x.start", "x.stop"))
  setDT(OVERLAPS)
  data.table::setkey(OVERLAPS, CHR, x.start, x.stop)
  
  ######################################################################################################################################################################
  # using for loop
  usedCallers <- c()
  DTlist <- rbindlist(list())
  if(length(cnvkit)!=0){usedCallers <- c(usedCallers,"cnvkit"); DTlist <- rbindlist(list(DTlist,cnvkitDF))}
  if(length(cnmops)!=0){usedCallers <- c(usedCallers,"cnmops");DTlist <- rbindlist(list(DTlist,cnmopsDF))}
  if(length(exomedepth)!=0){usedCallers <- c(usedCallers,"exomedepth");DTlist <- rbindlist(list(DTlist,exomedepthDF))}
  if(length(panelcnmops)!=0){usedCallers <- c(usedCallers,"panelcnmops");DTlist <- rbindlist(list(DTlist,panelcnmopsDF))}

  caller <- usedCallers[1]
  for(caller in usedCallers){
    
    tempDF <- DTlist[CALLER==caller,]
    setDT(tempDF)
    data.table::setkey(tempDF, chromosome, start, end)
    
    OVERLAPS <- data.table::foverlaps(OVERLAPS,tempDF, by.x = c("CHR", "x.start", "x.stop"), by.y = c("chromosome", "start", "end"), mult="all",type="any",nomatch=NA, which=FALSE)
    
    colnames(OVERLAPS) <- gsub("^start$",sprintf("cnmops_START",usedCallers),colnames(OVERLAPS))
    colnames(OVERLAPS) <- gsub("^end$","cnmops_STOP",colnames(OVERLAPS))
    # colnames(OVERLAPS) <- gsub("^CALLER$","cnmops_CALLER",colnames(OVERLAPS))
    
    data.table::setkey(OVERLAPS, CHR, x.start, x.stop)
    
    
    
    
  }
  
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  ## overlap 1
  OVERLAPS <- data.table::foverlaps(OVERLAPS,cnmopsDF, by.x = c("CHR", "x.start", "x.stop"), by.y = c("chromosome", "start", "end"), mult="all",type="any",nomatch=NA, which=FALSE)
  
  colnames(OVERLAPS) <- gsub("^start$","cnmops_START",colnames(OVERLAPS))
  colnames(OVERLAPS) <- gsub("^end$","cnmops_STOP",colnames(OVERLAPS))
  # colnames(OVERLAPS) <- gsub("^CALLER$","cnmops_CALLER",colnames(OVERLAPS))
  
  data.table::setkey(OVERLAPS, CHR, x.start, x.stop)
  
  
  ## overlap 2 
  OVERLAPS <- data.table::foverlaps(OVERLAPS,exomedepthDF, by.x = c("CHR", "x.start", "x.stop"), by.y = c("chromosome", "start", "end"), mult="all",type="any",nomatch=NA, which=FALSE)
  
  colnames(OVERLAPS) <- gsub("^start$","exomedepth_START",colnames(OVERLAPS))
  colnames(OVERLAPS) <- gsub("^end$","exomedepth_STOP",colnames(OVERLAPS))
  # colnames(OVERLAPS) <- gsub("^CALLER$","exomedepth_CALLER",colnames(OVERLAPS))
  
  data.table::setkey(OVERLAPS, CHR, x.start, x.stop)
  
  
  ## overlap 3 
  OVERLAPS <- data.table::foverlaps(OVERLAPS,panelcnmopsDF, by.x = c("CHR", "x.start", "x.stop"), by.y = c("chromosome", "start", "end"), mult="all",type="any",nomatch=NA, which=FALSE)
  
  colnames(OVERLAPS) <- gsub("^start$","panelcnmops_START",colnames(OVERLAPS))
  colnames(OVERLAPS) <- gsub("^end$","panelcnmops_STOP",colnames(OVERLAPS))
  # colnames(OVERLAPS) <- gsub("^CALLER$","panelcnmops_CALLER",colnames(OVERLAPS))
  
  data.table::setkey(OVERLAPS, CHR, x.start, x.stop)
  
  ## ## ## 
  
  ## select NEW START as maximal (the most rigth) STOP as minimal (the most left) coordinates
  OVERLAPS <- OVERLAPS[, START := do.call(pmax, c(.SD, list(na.rm=TRUE)) ), .SDcols = grep("_START$", names(OVERLAPS))]
  OVERLAPS <- OVERLAPS[, STOP := do.call(pmin, c(.SD, list(na.rm=TRUE)) ), .SDcols = grep("_STOP$", names(OVERLAPS))]
  OVERLAPS <- OVERLAPS[, LENGTH := STOP-START]
  
  # treat not mutually overlapped regions
  OVERLAPS[, START := ifelse(LENGTH<=0,  x.start , START)]
  OVERLAPS[, STOP := ifelse(LENGTH<=0,  x.stop , STOP)]
  OVERLAPS[, mutualOverlap := ifelse(LENGTH>0,  "yes" , "no")]
  OVERLAPS <- OVERLAPS[, LENGTH := STOP-START]
  
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################

  
  OVERLAPS$Callers <- ""
  OVERLAPS$CallersNum <- 0
  for (j in 1:nrow(OVERLAPS)){
    callers=""
    if(!is.na(OVERLAPS$cnmops_CALLER[j])){callers <- paste(callers,"cnmops",sep = ";")} # append to callers string
    if(!is.na(OVERLAPS$exomedepth_CALLER[j])){callers <- paste(callers,"exomedepth",sep = ";")} # append to callers string
    if(!is.na(OVERLAPS$panelcnmops_CALLER[j])){callers <- paste(callers,"panelcnmops",sep = ";")} # append to callers string
    
    callers <- gsub("^;","",callers)
    OVERLAPS$Callers[j] <- callers
    
    OVERLAPS$CallersNum[j] <- sapply(strsplit(callers, ";"), length)
    
  }
  
  OVERLAPS <- OVERLAPS[mutualOverlap=="yes",]
  
  setDT(OVERLAPS)
  OVERLAPS <- OVERLAPS[,c("cnmops_CALLER","exomedepth_CALLER","panelcnmops_CALLER",
                          "panelcnmops_START","panelcnmops_STOP","exomedepth_START","exomedepth_STOP",
                          "cnmops_START","cnmops_STOP",
                          "x.start","x.stop","mutualOverlap") :=NULL 
                       ]
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  # CN processing
  
  # exomedepth CN numbers
  OVERLAPS$exomedepth_Type <- ""
  OVERLAPS[OVERLAPS$exomedepth_CN==0,"exomedepth_Type"] <- "DEL"
  OVERLAPS[OVERLAPS$exomedepth_CN==10,"exomedepth_Type"] <- "DUP"
  OVERLAPS[,exomedepth_CN:=NULL]
  
  OVERLAPS$cnmops_Type <- ""
  OVERLAPS[OVERLAPS$cnmops_CN<2,"cnmops_Type"] <- "DEL"
  OVERLAPS[OVERLAPS$cnmops_CN>2,"cnmops_Type"] <- "DUP"
  OVERLAPS[,cnmops_CN:=NULL]
  
  OVERLAPS$panelcnmops_Type <- ""
  OVERLAPS[OVERLAPS$panelcnmops_CN<2,"panelcnmops_Type"] <- "DEL"
  OVERLAPS[OVERLAPS$panelcnmops_CN>2,"panelcnmops_Type"] <- "DUP"
  OVERLAPS[,panelcnmops_CN:=NULL]
  
  j=1
  #
  OVERLAPS$CNVtype <- ""
  for (j in 1:nrow(OVERLAPS)){
    
    # vector of cnv types
    types <- paste(OVERLAPS$exomedepth_Type[j],OVERLAPS$cnmops_Type[j],OVERLAPS$panelcnmops_Type[j],sep=";")
    x <- unlist(strsplit(types, ";"))
    x <- x[!x==""]
    # return the most frequent one
    OVERLAPS$CNVtype[j] <- names(sort(table(x),decreasing=TRUE,useNA = "no")[1])
    
    #
    # OVERLAPS[is.na(OVERLAPS$CNVtype) & OVERLAPS$panCNMOPS<2 & OVERLAPS$CNMOPS<2  ,"CNVtype"] <- "DEL"
    # OVERLAPS[is.na(OVERLAPS$CNVtype) & OVERLAPS$panCNMOPS>2 & OVERLAPS$CNMOPS>2  ,"CNVtype"] <- "DUP"
    
  }
  
  # names(OVERLAPS)
  setDT(OVERLAPS)
  OVERLAPS <- OVERLAPS[,c("exomedepth_Type","cnmops_Type","panelcnmops_Type") :=NULL ]
  BED <- OVERLAPS[,c("CHR","START","STOP","CNVtype")]
  #######################################################################################################################
  dir.create(file.path(getwd(),dirname(bed)), recursive = TRUE)
  write.table(BED, file=bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
  write.table(OVERLAPS, file=tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  
}


# # develop and test
# setwd("/mnt/ssd/ssd_1/workspace/ROBIN2/zewen_WES/")
# args <- c("CNV_exon_merged/14394_17_Tu_red.callers_merged.bed",
#           "CNV_exon_merged/14394_17_Tu_red.callers_merged.tsv",
#           "gencode.v42.basic.annotation.gtf.gz",
#           "variant_calls/14394_17_Tu_red/cnMOPS/cnMOPS_CNV_14394_17_Tu_red.tsv",
#           "variant_calls/14394_17_Tu_red/exomeDepth/exomeDepth_CNV_14394_17_Tu_red.tsv",
#           "variant_calls/14394_17_Tu_red/panelcnMOPS/panelcnMOPS_CNV_14394_17_Tu_red.tsv")


#run as Rscript

args <- commandArgs(trailingOnly = T)
run_all(args)