
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))

#nemelo by byt napevno... ale s detekci UCSC/Ensembl
listofCHR<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
             "chr19","chr20","chr21","chr22","chrY","chrX")

listofCHR2<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
              "19","20","21","22","X","Y")

# develop and test
# setwd("/media/rj/SSD_500GB/CNV_OVARIA/FFPE/")
# setwd("/media/rj/SSD_500GB/CNV_OVARIA/PLAZMY/")
# args <- c("CNV_Whole/1922-22_FFPE.callers_merged.bed","CNV_Whole/1922-22_FFPE.callers_merged.tsv","/home/rj/4TB/CEITEC/CNV_EXOM_BEDs/OvarianCancer_GRCh38.bed",
#           "variant_calls/1922-22_FFPE/cnMOPS/cnMOPS_CNV_1922-22_FFPE.tsv",
#           "variant_calls/1922-22_FFPE/exomeDepth/exomeDepth_CNV_1922-22_FFPE.tsv",
#           "variant_calls/1922-22_FFPE/panelcnMOPS/panelcnMOPS_CNV_1922-22_FFPE.tsv",
#           "variant_calls/1922-22_FFPE/cnvkit/cnvkit_CNV_1922-22_FFPE.tsv"
#           )
# args <- c("CNV_Whole/1850-23_plazma.callers_merged.bed","CNV_Whole/1850-23_plazma.callers_merged.tsv","/home/rj/4TB/CEITEC/CNV_EXOM_BEDs/OvarianCancer_capture_targets_GRCh38.bed","variant_calls/1850-23_plazma/cnMOPS/cnMOPS_CNV_1850-23_plazma.tsv","variant_calls/1850-23_plazma/exomeDepth/exomeDepth_CNV_1850-23_plazma.tsv","variant_calls/1850-23_plazma/panelcnMOPS/panelcnMOPS_CNV_1850-23_plazma.tsv")


run_all <- function(args){
  # arguments
  output_bed <- args[1]
  output_tsv <- args[2]
  bed_file <- args[3]
  tsv_files <- args[4:length(args)]
  threshold <- 2
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
    cnvkitDF <- fread(file=cnvkit, sep='\t', header = FALSE)
    setnames(cnvkitDF,c("chromosome","start","end","xx","CN"))
    if(nrow(cnvkitDF)){
      
      cnvkitDF[,xx:=NULL]
      cnvkitDF <- cnvkitDF[cnvkitDF$CN!=2,]
      
      if(UCSCorENSEMBL=="UCSC"){
        cnvkitDF$chromosome <- gsub("^chr","",cnvkitDF$chromosome) # carefull
        cnvkitDF <- cnvkitDF[chromosome %in% listofCHR2,]
      }else{
        cnvkitDF <- cnvkitDF[chromosome %in% listofCHR,]
      }
      # HANDLE CN AND CNV TYPE
      cnvkitDF$CN <- as.numeric(gsub("CN","",cnvkitDF$CN))
      cnvkitDF$TYPE <- "NA"
      cnvkitDF[cnvkitDF$CN<2,"TYPE"] <- "DEL"
      cnvkitDF[cnvkitDF$CN>2,"TYPE"] <- "DUP"
      
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
  
  
  ##############################################################################
  ##############################################################################
  ## MERGE WITH EACH OTHER
  
  # deletions - homo/heterezog
  # duplications
  # separately
  
  cnvkitDF$chromosome <- as.character(cnvkitDF$chromosome)
  setDT(cnvkitDF)
  setkey(cnvkitDF, chromosome,start,end)
  setnames(cnvkitDF,c("CN"),c("CN_cnvkit"))
  
  cnmopsDF$chromosome <- as.character(cnmopsDF$chromosome)
  setDT(cnmopsDF)
  setkey(cnmopsDF, chromosome,start,end)
  setnames(cnmopsDF,c("CN"),c("CN_cnmops"))
  
  exomedepthDF$chromosome <- as.character(exomedepthDF$chromosome)
  setDT(exomedepthDF)
  setkey(exomedepthDF, chromosome,start,end)
  setnames(exomedepthDF,c("CN"),c("CN_exomedepthDF"))
  
  panelcnmopsDF$chromosome <- as.character(panelcnmopsDF$chromosome)
  setDT(panelcnmopsDF)
  setkey(panelcnmopsDF, chromosome,start,end)
  setnames(panelcnmopsDF,c("CN"),c("CN_pcnmops"))
  
  
  ##############################################################################
  # COMBINE INTO SINGLE DT
  OVERLAPS <- rbind(cnvkitDF,cnmopsDF,exomedepthDF,panelcnmopsDF,fill=TRUE)
  setDT(OVERLAPS)
  names(OVERLAPS)
  
  ##############################################################################
  # Assuming TYPE column exists in cnv_data
  sv_types <- unique(OVERLAPS$TYPE)
  
  # Create an empty list to collect final results
  final_list <- list()
  
  for (sv in sv_types) {
    message("Processing SV type: ", sv)
    
    # Subset CNVs of current TYPE
    subset_data <- OVERLAPS[TYPE == sv]
    
    # Convert to GRanges
    gr <- GRanges(seqnames = subset_data$chromosome,
                  ranges = IRanges(start = subset_data$start, end = subset_data$end),
                  caller = subset_data$CALLER,
                  cnCNVkit = subset_data$CN_cnvkit,
                  cnCNMOPS = subset_data$CN_cnmops,
                  cnPCNMOPS = subset_data$CN_pcnmops
    )
    
    # Reduce overlapping CNVs
    reduced_gr <- reduce(gr)
    
    # Overlap between reduced and original GRanges
    overlaps <- findOverlaps(reduced_gr, gr)
    
    # Build a data.table of overlapping metadata
    overlap_dt <- data.table(
      merged_id = queryHits(overlaps),
      caller = mcols(gr)$caller[subjectHits(overlaps)],
      cnCNVkit = mcols(gr)$cnCNVkit[subjectHits(overlaps)],
      cnCNMOPS = mcols(gr)$cnCNMOPS[subjectHits(overlaps)],
      cnPCNMOPS = mcols(gr)$cnPCNMOPS[subjectHits(overlaps)]
    )
    
    # Aggregate per merged_id
    summary_dt <- overlap_dt[, .(
      callers = list(unique(caller)),
      n_callers = uniqueN(caller),
      CN_cnvkit = paste(na.omit(cnCNVkit), collapse = ","),
      CN_cnmops = paste(na.omit(cnCNMOPS), collapse = ","),
      CN_pcnmops = paste(na.omit(cnPCNMOPS), collapse = ",")
    ), by = merged_id]
    
    # Apply threshold filter
    filtered <- summary_dt[n_callers >= threshold]
    
    # Get final merged regions
    final_regions <- reduced_gr[filtered$merged_id]
    
    if(nrow(filtered)>0){
    # Create final dataframe for this TYPE
      final_df <- data.frame(
        seqnames = seqnames(final_regions),
        start = start(final_regions),
        end = end(final_regions),
        n_callers = filtered$n_callers,
        callers = sapply(filtered$callers, paste, collapse = ","),
        CN_cnvkit = filtered$CN_cnvkit,
        CN_cnmops = filtered$CN_cnmops,
        CN_pcnmops = filtered$CN_pcnmops,
        svtype = sv
      )
      
      final_list[[sv]] <- final_df
    }else{
      final_df <- data.frame(
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        n_callers = integer(0),
        callers = character(0),
        CN_cnvkit = integer(0),
        CN_cnmops = integer(0),
        CN_pcnmops = integer(0),
        svtype =  character(0)
      )
      final_list[[sv]] <- final_df
    }
    
  }
  
  # Combine all results into one dataframe
  all_results <- rbindlist(final_list)
  # length
  all_results$length <- all_results$end-all_results$start
  # rename
  names(all_results)
  setnames(all_results,
           c("seqnames","start","end","n_callers","callers","CN_cnvkit","CN_cnmops","CN_pcnmops","svtype","length"),
           c("CHR","START","STOP","n_CALLERS","CALLERS","CN_cnvkit","CN_cnmops","CN_pcnmops","CNVtype","CNVlength"))


  setDT(OVERLAPS)
  OVERLAPS <- all_results[,c("CHR","START","STOP","CNVlength","CNVtype","n_CALLERS","CALLERS","CN_cnvkit","CN_cnmops","CN_pcnmops")]
  BED <- OVERLAPS[,c("CHR","START","STOP","CNVtype")]
  
  ##############################################################################
  dir.create(file.path(getwd(), dirname(output_bed)), recursive = TRUE, showWarnings = FALSE)
  write.table(BED, file=output_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )
  write.table(OVERLAPS, file=output_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  
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