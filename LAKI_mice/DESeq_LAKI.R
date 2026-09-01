#///////////////////////////////////////////////////////////////////////////////
############################## DESeq Analysis ##################################
#///////////////////////////////////////////////////////////////////////////////

setwd('/home/antotartier/data/20231017_RNASeq_LAKI/')

# Variables
out_path <- 'results/transposons/subfamily/'
out_object <- 'R_files/transposons/subfamily/'
counts_path <- Sys.glob('TE_count/count_files_unstr/*cntTable')
sample_info <- read.delim('R_files/meta/Sample_info_LAKI_corrected.tsv')

# Libraries
require(DESeq2)
require(readr)
require(dplyr)
require(tibble)

#///////////////////////////////////////////////////////////////////////////////
############################## DESeq Dataset ##################################
#///////////////////////////////////////////////////////////////////////////////

  #Select count data present in sample info
  counts_path<-counts_path[grepl(paste0(sample_info$SampleName,'\\.',collapse = '|'),counts_path)]
  
  #Read count data
  countData<-read_delim(counts_path, id = "sample", col_names = F, skip = 1,show_col_types = F)
  countData<-pivot_wider(countData,names_from = sample, values_from = X2)
  colnames(countData)<-gsub(".cntTable","",basename(colnames(countData)))
  countData<-data.frame(countData, row.names = 1)
  rownames(countData)<-gsub("\\.[0-9]*","",rownames(countData))
  countData<-countData[,order(as.numeric(gsub("Tube_","",colnames(countData))))]
  
  #Create colData
  colData<-sample_info |>
    column_to_rownames(var = 'SampleName') |>
    dplyr::select(-SampleID,-InferedGeno)
  
  #Convert in factor and determine the passage order 
  colData$Sex<-factor(colData$Sex)
  colData$Genotype<-factor(colData$Genotype,levels = c('WT','KO'))
  
  #Save colData
  saveRDS(colData,paste0(out_object,'colData.RDS'))
  
  #Generate DESeqDataSet
  dds_data<-DESeqDataSetFromMatrix(countData,colData,design = ~Sex + Genotype)
  
  #Add annotation to rowdata
  annotation_data<-readRDS('/data/antotartier/annotation/gencode_vM34_primary_annot.RDS')
  rnames<-data.frame(gene_id = rownames(dds_data))
  row_anot<-left_join(rnames,annotation_data,by = 'gene_id') |>
    dplyr::select(gene_name)
  rowData(dds_data)<-row_anot
  
  #separate tissues
  tissues<-unique(dds_data$Tissue)
  dds_data_tis<-list()
  for(x in tissues){
    subset_data<-dds_data[,dds_data$Tissue==x]
    dds_data_tis[[x]]<-subset_data
  }
  
  #Filtering those genes with at least 10 lectures in the smaller discrete group (KO/WT)
  dds_data_tis_filt<-list()
  for(name in names(dds_data_tis)){
    y <-dds_data_tis[[name]]
    smallestGroupSize <- min(c(sum(y$Genotype=="KO"),sum(y$Genotype=="WT")))
    y <- y[rowSums(counts(y) >= 10) >= smallestGroupSize,]
    dds_data_tis_filt[[name]]<-y
  }
  
  
# Sample 86 from colon has low quantity of reads and seems an outlier so we eliminate it from the analysis
  
  #Eliminate sample 86 and save dds dataset
  dds_data_tis$Colon<-dds_data_tis$Colon[,!colnames(dds_data_tis$Colon)=='Tube_86']
  write_rds(dds_data_tis, paste0(out_object,'dds_data.rds'))
  
  #Filter
  dds_data_tis_filt<-list()
  for(name in names(dds_data_tis)){
    y <-dds_data_tis[[name]]
    smallestGroupSize <- min(c(sum(y$Genotype=="KO"),sum(y$Genotype=="WT")))
    y <- y[rowSums(counts(y) >= 10) >= smallestGroupSize,]
    dds_data_tis_filt[[name]]<-y
  }


#///////////////////////////////////////////////////////////////////////////////
############################## DESeq Experiment ################################
#///////////////////////////////////////////////////////////////////////////////
  
  #DESeq experiment
  dds_exp<-lapply(dds_data_tis_filt,DESeq)
  saveRDS(dds_exp,file = paste0(out_object,'dds_exp.rds'))
  

#///////////////////////////////////////////////////////////////////////////////
############################## DESeq Result ################################
#///////////////////////////////////////////////////////////////////////////////
  
  #Generating results and saving
  res<-lapply(dds_exp,results,saveCols = 'gene_name')
  
  saveRDS(res,file = paste0(out_object,'res.rds'))
  res_df<-lapply(res,as.data.frame)
  write.xlsx(res_df,paste0(out_path,"res.xlsx"),rowNames=T)