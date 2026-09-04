#///////////////////////////////////////////////////////////////////////////////
###################### Histone Expression Analysis #############################
#///////////////////////////////////////////////////////////////////////////////

# Variables
out_path <- 'results/histones/'
out_object <- 'R_files/histones/'
res <- readRDS('R_files/transposons/subfamily/res.rds')
term_to_gene <- readRDS('/data/antotartier/gene_sets/Histone_genes_GRCm39_GENCODE.rds')

# Libraries
require(openxlsx)
require(dplyr)
require(tibble)
require(ggsignif)

#///////////////////////////////////////////////////////////////////////////////
###################### Histone Differential Expression #########################
#///////////////////////////////////////////////////////////////////////////////

# I will inspect DESeq results to see if LAKI mice present changes in histone expression

#histone genes list from HGNC
histones<-read.xlsx("R_files/histones/gene_set_histonas_raton.xlsx")

#classify histones
dpts<-c("cluster","H2ax","H4c16")
histclas <- histones |>
  mutate(type = case_when(grepl('pseudogene', Marker_Name) ~ 'pseudogene',
                          grepl('lnc', BioTypes) ~ 'lncRNA',
                          !grepl('pseudogene', Marker_Name) ~ 'gene'),
         regulation = ifelse(grepl(paste(dpts, collapse = "|"),Marker_Name) & !grepl("H3-4", Marker_Name),"RDH",'RIH'),
         family = case_when(grepl('H1', Marker_Name) ~ 'H1',
                            grepl('H2A', Marker_Name) ~ 'H2a',
                            grepl('H2B', Marker_Name) ~ 'H2b',
                            grepl('H3', Marker_Name) ~ 'H3',
                            grepl('H4', Marker_Name) ~ 'H4'
         )) |> 
  dplyr::select(Marker_Symbol, Marker_Name, Ensembl_Accession_ID, type, regulation, family)

#save ID and names
histone_name_id<- dplyr:::select(histones, Marker_Symbol, Ensembl_Accession_ID)
saveRDS(histone_name_id,file = paste0(out_object,"histone_name_id.RDS"))

#extract histone expression
res_df<-lapply(res,as.data.frame) |> 
  lapply(rownames_to_column, var = 'Ensembl_Accession_ID')
res_hist<-lapply(res_df, \(x) inner_join(x, histclas, by = 'Ensembl_Accession_ID'))

#eliminate possible NA
res_hist<-lapply(res_hist,na.omit)

#Fuse data
res_hist<- bind_rows(res_hist,.id = 'tissue')

#filter histone genes
res_hist_genes<-res_hist[res_hist$type=="gene",]

#compare expression changes between RDH and RIH
ggplot(res_hist_genes, aes(regulation, log2FoldChange, fill= regulation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, position = 'jitter') +
  geom_abline(slope = 0, intercept = 0, linetype = 'dashed' ) +
  geom_signif(test = 't.test', comparisons = list(c('RDH', 'RIH'))) +
  theme_bw() +
  facet_wrap(~tissue)

#compare expression changes between histone families
ggplot(res_hist_genes, aes(family, log2FoldChange, fill= family)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, position = 'jitter') +
  geom_abline(slope = 0, intercept = 0, linetype = 'dashed' ) +
  theme_bw() +
  facet_wrap(~tissue)

#///////////////////////////////////////////////////////////////////////////////
################################# Histone GSEA ################################
#///////////////////////////////////////////////////////////////////////////////

# To confirm histone downregulation we perform a GSEA of histone genes

#set seed
set.seed(10)

#select coding genes
cod <- lapply(res, \(x) filter(as.data.frame(x), !grepl(':', rownames(x))))

#generating rank
rank <- lapply(cod, \(x) na.omit(reframe(x, ranking = -log10(pvalue+min(pvalue[pvalue!=0], na.rm= T))*sign(log2FoldChange), ensembl_gene_id = rownames(x)))) |> 
  lapply(\(x) setNames(x$ranking, x$ensembl_gene_id)) |> 
  lapply(sort, decreasing = T)

#execute GSEA
hisGSEA <- lapply(rank, GSEA, TERM2GENE = term_to_gene, pvalueCutoff = 1)

#plot
GSEA_plot<-mapply(
    function(x, plotTitle){
      gseaplot(x,
               geneSetID = "Histone genes",
               by = 'runningScore',
               title = paste(plotTitle,x[1][1],'|', 'P.adjusted =',x[1][7],sep = ' ')
      )
    },
    hisGSEA, names(rank), SIMPLIFY = F)
patchwork::wrap_plots(GSEA_plot)
