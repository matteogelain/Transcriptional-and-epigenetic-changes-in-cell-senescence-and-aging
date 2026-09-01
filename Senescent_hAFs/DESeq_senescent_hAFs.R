#///////////////////////////////////////////////////////////////////////////////
############################# DESEq Analysis ###################################
#///////////////////////////////////////////////////////////////////////////////

# We start the analysis from RNAseq counts previousy generated in the laboratory

#Libraries
require(readr)
require(DESeq2)
require(tidyr)
require(openxlsx)
require(textclean)
require(dplyr)

#Variables
setwd('data/fib_senescent_RNAseq/')
out_path<-'results/transposons/subfamily/'
out_object<-'R_files/transposons/subfamily/'
counts_path<- Sys.glob(paste('transposons/TEcount/count_files_unstr', '/*cntTable', sep = ''))

#///////////////////////////////////////////////////////////////////////////////
############################# DESEq Dataset ###################################
#///////////////////////////////////////////////////////////////////////////////

#Read count data
countData<-read_delim(counts_path, id = "sample", col_names = F, skip = 1,show_col_types = F)
countData<-pivot_wider(countData,names_from = sample, values_from = X2)
colnames(countData)<-gsub(".cntTable","",basename(colnames(countData)))
countData<-data.frame(countData, row.names = 1)
rownames(countData)<-gsub("\\.[0-9]*","",rownames(countData))
countData<-countData[,order(as.numeric(gsub("Tube_","",colnames(countData))))]

#Create colData
sample_info<-read.xlsx("sample_info.xlsx")
colData<-sample_info[c(1:9,19:27),2:3]
colData$Passage<-mgsub(colData$Passage, c("p14","p20","p24","p21","p27","p34"), c("Early","Mid","Late","Early","Mid","Late"))
rownames(colData)<-colnames(countData)

#Convert in factor and determine the passage order
colData$Cells<-factor(colData$Cells)
colData$Passage<-factor(colData$Passage,levels = c("Early","Mid","Late"))

#Generate DESeqDataSet
dds_data<-DESeqDataSetFromMatrix(countData,colData,design = ~Cells + Passage)

#Add annotation to rowdata
annotation_data<-readRDS('/data/antotartier/annotation/gencode_v44_primary_annot.RDS')
rnames<-data.frame(gene_id = rownames(dds_data))
row_anot<-left_join(rnames,annotation_data,by = 'gene_id') |>
  dplyr::select(gene_name)
rowData(dds_data)<-row_anot

#Separate cell lines
dds_data_indep<-list(GM05565 = dds_data[, dds_data$Cells == "GM05565"],
                     GM00038 = dds_data[, dds_data$Cells == "GM00038"])
design(dds_data_indep$GM05565)<-formula(~Passage)
design(dds_data_indep$GM00038)<-formula(~Passage)

#Filtering
smallestGroupSize <- 3 #3 samples/passage
for (name in names(dds_data_indep)){
  keep <- rowSums(counts(dds_data_indep[[name]]) >= 10) >= smallestGroupSize
  dds_data_indep[[name]] <- dds_data_indep[[name]][keep,]
}


#///////////////////////////////////////////////////////////////////////////////
############################# DESEq Experiment #################################
#///////////////////////////////////////////////////////////////////////////////

#DEseq experiment
dds_exp_indep<-lapply(dds_data_indep,DESeq,test = 'LRT', reduced = ~1)
saveRDS(dds_exp_indep,file = paste0(out_object,'dds_exp_indep.rds'))

#///////////////////////////////////////////////////////////////////////////////
############################# DESEq Results #################################
#///////////////////////////////////////////////////////////////////////////////

#Generating results
res_indep<-lapply(dds_exp_indep,results,saveCols = 'gene_name')

saveRDS(res_indep,file = paste0(out_object,'res_indep.rds'))
res_indep_df<-lapply(res_indep,as.data.frame)
write.xlsx(res_indep_df,paste0(out_path,"DESeq_res_indep.xlsx"),rowNames=T)
