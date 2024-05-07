suppressMessages(library(data.table))
suppressMessages(library(ExomeDepth))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")


################################################################################
## RUN
run_all <- function(args){
  #read all params as variables
  BAM_file <- args[1]
  lib_ROI <- args[2]
  sampleNameWildCard <- args[3]
  fastaFile <- args[4]
  output_file <- args[5]
  
  if(lib_ROI=="HyperExome_GRCh38"){
    Rdata_file <- paste0(script_dir,"/COHORT_exonsHyperExome_exomeDepth.Rdata")
    load(Rdata_file)
    exons <- exons.GRCh38
  }
  if(lib_ROI=="TruSeq_Exome"){
    Rdata_file <- paste0(script_dir,"/COHORT_exonsTruSeqExome_exomeDepth.Rdata")
    load(Rdata_file)
  }
  if(lib_ROI=="HyperExome_w_CNVbb_GRCh38"){
    Rdata_file <- paste0(script_dir,"/COHORT_exonsHyperExome_w_CNVbb_exomeDepth.Rdata")
    load(Rdata_file)
    exons <- exons.GRCh38
  }
  # if(lib_ROI=="XXX"){
  #   load("COHORT_exonsXXXX_cn.MOPS.Rdata")
  #   exons <- exons.GRCh38
  # }
  
  
  ##############################################################################
  ## LOAD SAMPLE
  X_sample <- getBamCounts(bed.frame = exons,
                           bam.files = BAM_file,
                           include.chr = FALSE,
                           referenceFasta = fastaFile)
  
  
  # rename column
  # sampleNameWildcard2 <- gsub("-",".",sampleNameWildCard)
  # idx <- grep(paste0("^",sampleNameWildcard2), colnames(X_sample))
  idx <- grep(".bam$", colnames(X_sample))
  names(X_sample)[idx] <- "SAMPLENAME"
  
  # merge with X_cohort
  X_cohort <- merge(X_cohort,X_sample,by=c("chromosome", "start", "end","exon","GC"))
  
  ##############################################################################
  # Build the most appropriate reference set
  # for selected sample
  all.samples <- colnames(X_cohort)[5:length(colnames(X_cohort))]
  my.ref.samples <- setdiff(all.samples,"SAMPLENAME")
  
  # subset dataframes
  my.test <- X_cohort[,colnames(X_cohort)=="SAMPLENAME"]
  my.reference.set <- as.matrix(X_cohort[, my.ref.samples])
  my.choice <- select.reference.set (test.counts = my.test,
                                     reference.counts = my.reference.set,
                                     bin.length = (X_cohort$end - X_cohort$start)/1000,
                                     n.bins.reduced = 10000)
  
  # selected correlated samples as reference used
  #print(my.choice[[1]])
  
  # Create the aggregate reference set for this sample
  my.matrix <- as.matrix( X_cohort[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(X = my.matrix,
                                 MAR = 1,
                                 FUN = sum)
  
  
  all.exons <- new('ExomeDepth',
                   test = my.test,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')
  
  ##############################################################################
  ## CNV calling
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = X_cohort$chromosome,
                        start = X_cohort$start,
                        end = X_cohort$end,
                        name = X_cohort$exon)
  
  
  DF <- all.exons@CNV.calls
  # length column 
  DF$length <- DF$end-DF$start
  DF <- subset(DF, select=c("chromosome","start","end", "length",  "id", "type" ,        
                            "nexons","BF",
                            "reads.expected","reads.observed","reads.ratio"
  ))
  

  if(nrow(DF)>0){
    dir.create(file.path(dirname(output_file)), showWarnings = FALSE)
    fwrite(DF,file = output_file,sep = "\t")
  }else{
    DF <- data.frame(chromosome=character(0),start=character(0),end=character(0),length=character(0),id=character(0),type=character(0),
                     nexons=character(0),BF=character(0),reads.expected=character(0),reads.observed=character(0),reads.ratio=character(0))
    dir.create(file.path(dirname(output_file)), showWarnings = FALSE)
    fwrite(DF,file = output_file,sep = "\t")
    
  }

  
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