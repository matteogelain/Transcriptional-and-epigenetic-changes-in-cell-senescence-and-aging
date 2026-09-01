#///////////////////////////////////////////////////////////////////////////////
######################### Histone Expression Analysis ##########################
#///////////////////////////////////////////////////////////////////////////////

# Once analyzed gene expression with DESeq we proceed to explore the results

# Libraries
require(dplyr)
require(ggplot2)
require(clusterProfiler)

# Variables
out_path <- 'results/histones/'
out_object <- 'R_files/histones/'
res <- readRDS('R_files/transposons/subfamily/res_indep.rds')
sample_info <- read.xlsx("sample_info.xlsx")


#///////////////////////////////////////////////////////////////////////////////
######################### Histone Differential Expression ######################
#///////////////////////////////////////////////////////////////////////////////

#histone genes list from HGNC
histones<-read.xlsx("R_files/histones/gene_set_histonas.xlsx")

#classify histone
#define which are the replication dependent histones
dpts<-c("cluster","H2AX","H4C16")
histones <- histones |>
  mutate(type = ifelse(grepl("pseudogene", name) | !is.na(histones$pseudogene.org), "pseudogene", 'gene'),
         regulation = ifelse(grepl(paste(dpts, collapse = "|"),name) & !grepl("H3-4", name),"RDH",'RIH'),
         family = case_when(grepl('H1', name) ~ 'H1',
                           grepl('H2A', name) ~ 'H2A',
                           grepl('H2B', name) ~ 'H2B',
                           grepl('H3', name) ~ 'H3',
                           grepl('H4', name) ~ 'H4'
                           ))
histclas <- dplyr::select(histones, name, symbol, ensembl_gene_id, type, regulation, family)
write_rds(histclas,'R_files/histones/classified_histones.RDS')

#save ID and names
histone_name_id<-histones[,c(2,20)]
saveRDS(histone_name_id,file = paste0(out_object,"histone_name_id.RDS"))

#extract histone expression
res_df<-lapply(res, as.data.frame) |> 
  lapply(tibble::rownames_to_column, var = 'ensembl_gene_id') |> 
  bind_rows(.id = 'cell_line')
res_hist<-inner_join(res_df, histclas, by = 'ensembl_gene_id')

#Filter histone genes
res_hist_genes<-res_hist[res_hist$type=="gene",]

#compare expression changes between RDH and RIH and histone families
#Changes between RDH RIH
ggplot(res_hist_genes, aes(regulation, log2FoldChange, colour= regulation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, position = 'jitter') +
  theme_classic() +
  facet_wrap(~cell_line)

# Changes between families
ggplot(res_hist_genes, aes(family, log2FoldChange, colour= family)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, position = 'jitter') +
  theme_classic() +
  facet_wrap(~cell_line)

#///////////////////////////////////////////////////////////////////////////////
######################### Histone GSEA ######################
#///////////////////////////////////////////////////////////////////////////////

# To confirm histone downregulation we perform a GSEA of histone genes

#Set seed
seed(10)

# First we separate the coding genes
cod <- lapply(res, \(x) filter(as.data.frame(x), !grepl(':', rownames(x))))

# Generate ranking for gsea using pvalue
rank <- lapply(cod, \(x) na.omit(reframe(x, ranking = -log10(pvalue+min(pvalue[pvalue!=0], na.rm= T))*sign(log2FoldChange), ensembl_gene_id = rownames(x)))) |> 
  lapply(\(x) setNames(x$ranking, x$ensembl_gene_id)) |> 
  lapply(sort, decreasing = T)

# Generate term to gene
his <- filter(histclas, type == 'gene') |>
  select(ensembl_gene_id) |> 
  rename(gene = ensembl_gene_id) |> 
  mutate(term = 'RDH')
his <- his[, c("term", "gene")]
  

# Execute GSEA
hisGSEA <- lapply(rank, GSEA, TERM2GENE = his)

# Plots
p <- lapply(hisGSEA, \(x) gseaplot(x, geneSetID = 'RDH', by = 'runningScore', title = paste0('P.adjusted = ', as.data.frame(x)$p.adjust)))
