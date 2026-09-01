#///////////////////////////////////////////////////////////////////////////////
###################### Histone specific genes Analysis #########################
#///////////////////////////////////////////////////////////////////////////////

# Histone qPCR and Western in liver did not showed the clear downregulation of histones that we thought to see in RNAseq, therefore I am going to carefully explore the results, as it might be explained bu the contribution of each specific histone gene to the overall reads of its histone family 

# Libraries
require(dplyr)
require(ggplot2)

# Charge DESeq results and histone gene set
res <- readRDS('R_files/transposons/subfamily/res.rds')
hist_genes <- readRDS('R_files/histones/clasified_histones.rds')

# Extract histone results, add metadata and fuse them in a single dataframe and stablish up and downregulated genes
hist_res <- lapply(res, function(x)
{as.data.frame(x)[rownames(x) %in% hist_genes$Ensembl_Accession_ID,] |> 
    tibble:::rownames_to_column(var = 'Ensembl_Accession_ID') |> 
    left_join(hist_genes, by = 'Ensembl_Accession_ID')}) |> 
  bind_rows(.id = 'tissue') |> 
  mutate(
    status = case_when(
      log2FoldChange > 0 & padj < 0.05 ~ 'Upregulated',
      log2FoldChange < 0 & padj < 0.05 ~ 'Downregulated',
      padj > 0.05 ~ 'ns'),
    facet = interaction(tissue, family, sep = " - ")
  )

# Extract RDH histones
RDH <- filter(hist_res, reg == 'RDH' & type == 'gene') |> 
  group_by(facet) |> 
  mutate(perc_bmean = baseMean/sum(baseMean)) |> 
  ungroup()

# Plot of percentage of reads provided by each gene
ggplot(RDH |> filter(facet %in% c('Liver - H2b', 'Liver - H3')), aes(gene_name, perc_bmean, fill = log2FoldChange)) +
  geom_bar(aes(fill = log2FoldChange), stat = 'identity') +
  scale_fill_gradient2(
    low = '#e51f1f',
    mid = 'grey',
    high = '#44ce1b',
    midpoint = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315, hjust = 0),
        axis.title.x = element_blank()) +
  ylab('% of reads') +
  facet_wrap(~facet, ncol = 6, scales = 'free')

# Calculate and plot mean LFC for each histone family 
dds <- readRDS('R_files/transposons/subfamily/dds_exp.rds')

cts <- lapply(dds, DESeq2::counts, normalized = T) |> 
  lapply(as.data.frame) |> 
  lapply(tibble:::rownames_to_column, var = 'Ensembl_Accession_ID') |> 
  purrr:::reduce(full_join, by = 'Ensembl_Accession_ID')
hist_cts <- filter(cts, Ensembl_Accession_ID %in% hist_genes$Ensembl_Accession_ID) |> 
  left_join(hist_genes, by = 'Ensembl_Accession_ID')
RDH_cts <- filter(hist_cts, reg == 'RDH')

smpinf <- data.table::fread('R_files/meta/Sample_info_LAKI_corrected.tsv')

fam_normcts <- RDH_cts |> 
  tidyr::pivot_longer(starts_with('Tube'), names_to = 'SampleName', values_to = 'Norm_cts') |> 
  left_join(smpinf, by = 'SampleName') |> 
  group_by(family, SampleName) |> 
  summarise(
    Total_cts = sum(Norm_cts, na.rm = T),
    Genotype = unique(Genotype),
    Sex = unique(Sex),
    Tissue = unique(Tissue),
    Facet = interaction(unique(family), unique(Tissue), sep = '-'),
    .groups = 'drop')

fam_LFC <- fam_normcts|> 
  group_by(family, Genotype, Tissue) |> 
  summarise(Mean_cts = mean(Total_cts),
            .groups = 'drop') |> 
  tidyr::pivot_wider(names_from = Genotype, values_from = Mean_cts) |> 
  reframe(family = family,
          Tissue = Tissue,
          LFC = log2(KO/WT))

ggplot(fam_LFC |> filter(Tissue == 'Liver'), aes(family, LFC, fill = LFC)) +
  geom_col() +
  scale_fill_gradient2(high = 'darkgreen', low = '#e51f1f', mid = 'grey', midpoint = 0) +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        legend.title = element_text(size = 30),
        legend.text  = element_text(size = 15),
        legend.key.size = unit(1.2, "cm")) +
  geom_hline(yintercept = c(-1, 1), linetype = 'dashed')