
library(data.table)
library(rtracklayer)

#nemelo by byt napevno... ale s detekci UCSC/Ensembl
listofCHR<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
             "chr19","chr20","chr21","chr22","chrY","chrX")

listofCHR2<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
              "19","20","21","22","X","Y")

# develop and test
# setwd("/media/rj/SSD_500GB/CNV_OVARIA/FFPE/")
# args <- c("variant_calls/4434-2022_FFPE/cnvkit/cnvkit_CNV_preannot_4434-2022_FFPE.tsv","variant_calls/4434-2022_FFPE/cnvkit/4434-2022_FFPE.classified.txt","/home/rj/4TB/CEITEC/homsap/GRCh38/annot/r111/GRCh38.gtf","cnvkit_final/4434-2022_FFPE_final_CNVs_annotated.tsv","cnvkit_final/4434-2022_FFPE_final_CNVs_annotated.xlsx","")

run_all <- function(args){
  # arguments
  gtf_file <- args[1]
  output_tsv <- args[2]

  ##############################################################################
  # PROCESS GTF
  gtf_GENCODE <- rtracklayer::import(gtf_file)
  gtf_df_all <- as.data.table(gtf_GENCODE)

  if(any(grepl( "^type$",    names(gtf_df_all))) ){
    gtf_df_all <- gtf_df_all[gtf_df_all$type=="gene",]   #gtf_df_all <- gtf_df_all[gtf_df_all$type=="exon",]
  }
  
  if(any(grepl( "^gene_type$",    names(gtf_df_all))) ){
    gtf_df_all <- gtf_df_all[gtf_df_all$gene_type=="protein_coding",]
  }
  
  if(any(grepl( "^gene_biotype$",    names(gtf_df_all))) ){
    gtf_df_all <- gtf_df_all[gtf_df_all$gene_biotype=="protein_coding",]
  }
  
  
  if(all(c("seqnames","start","end","width","gene_name") %in% names(gtf_df_all))){
    setnames(gtf_df_all,"seqnames","chromosome")
    gtf_df_all <- subset(gtf_df_all,select=c("chromosome","start","end","gene_name"))
    gtf_df_all <- unique(gtf_df_all, by=c("chromosome","start","end","gene_name"))
    gtf_df_all$chromosome <- as.character(gtf_df_all$chromosome)
  }

  setDT(gtf_df_all)


  #######################################################################################################################
  write.table(gtf_df_all, file=output_tsv, sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE )
  
}

#run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)