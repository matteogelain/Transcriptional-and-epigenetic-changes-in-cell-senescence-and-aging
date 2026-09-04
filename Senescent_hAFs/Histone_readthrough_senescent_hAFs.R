#///////////////////////////////////////////////////////////////////////////////
######################### Histone Readthrough  Analysis ########################
#///////////////////////////////////////////////////////////////////////////////

# Some studies have described aberrant polyadenylation of histones as a driver of their downregulation both at RNA and protein level. However, the specific role of this process in cellular senescence remains unknown. To this end, we will analyze histone readthrough as an indirect measure of this phenomenom in RNAseq data from senescent human fibroblasts.

# Libraries
require(GenomicAlignments)
require(GenomicRanges)
require(rtracklayer)
require(ggplot2)
require(msigdbr)

## Functions
load_reads <- function(bafile, grs){
  bafile <- BamFile(bafile,asMates = T)
  
  # Define the params
  params <- ScanBamParam(
    what = c('qname'),
    which= grs,
    flag = scanBamFlag(
      isPaired = T,
      isProperPair = T,
      isSecondaryAlignment = F),
    mapqFilter = 5
  )
  # Load the reads
  reads <- readGAlignmentPairs(bafile, param = params,strandMode = 2)
  # Remove duplicated reads
  qnames <- mcols(first(reads))$qname
  reads <- reads[!duplicated(qnames)]
  # REmove reads with junctions
  # reads <- reads[njunc(reads) == 0]
  return(reads)
}

process_file <- function(bamfile, cdsgr, downgr){
  # Extend regions to get complete pairs, also those were read 1 is very far way
  extgr <- reduce(resize(c(cdsgr,downgr),fix='center',width=2E6), ignore.strand= TRUE)
  rds <- load_reads(bamfile, extgr) 
  # Pairs where mate2 overlaps cds and the pair range span downstream of cds, surely demonstrate readthrough
  m2incds <- findOverlaps(second(rds),cdsgr) # mate 2 in cds
  overDown <- findOverlaps(granges(rds),resize(downgr, 500, fix = 'start')) # Here I don't  need to use the annotation, since forcing mate2 to be in cds ensures that the read comes from the histone and not other gene
  cdstodown <- intersect(m2incds,overDown) # the intersection into these sets gives the pairs with mate 2 in cds and mate 1 in downstream region
  # Pairs where the mate2 falls in the downstream region defined by the annotation
  m2indown <- findOverlaps(second(rds), downgr)
  # Remove pairs where mate2 also overlap the cds (those are cdstodown pairs)
  onlyindown <- m2indown[!(from(m2indown) %in% from(cdstodown))]
  # Count number of reads with splice junctions in each read
  sj_cds2down <- njunc(rds[from(cdstodown)]) 
  sj_onlydown <- njunc(rds[from(onlyindown)])
  # determiner the number of reads with at least one splice junction in each gene
  sj_cds2down <- aggregate(sj_cds2down, by = list(gene = names(cdsgr[to(cdstodown)])), \(x) sum(x>0))
  sj_onlydown <- aggregate(sj_onlydown, by = list(gene = names(cdsgr[to(onlyindown)])), \(x) sum(x>0)) 
  # Return
  df <- data.frame(
    'gene' = histgenes, 
    'cds' = countOverlaps(cdsgr,rds),
    'mate2incds' = countRnodeHits(m2incds),
    'cds2down' = countRnodeHits(cdstodown),
    'onlydown' = countRnodeHits(onlyindown),
    'downstot' = countRnodeHits(cdstodown)+countRnodeHits(onlyindown)
  )
  df <- merge(df, sj_cds2down, by = 'gene', all.x = T) |> 
    merge(sj_onlydown, by = 'gene', all.x = T)
  colnames(df) <- c("gene", "cds", "mate2incds", "cds2down", "onlydown", "downstot", "sj_cds2down", "sj_onlydown")
  df[is.na(df)] <- 0# If there are no cds to down or only down genes the merge generates NAs that correspond to 0
  df$sj_downstot <- df$sj_cds2down + df$sj_onlydown
  return(df)
}


#///////////////////////////////////////////////////////////////////////////////
######################### Histone Readthrough Measure ########################
#///////////////////////////////////////////////////////////////////////////////

# Load the annotation
histgtf <- import('R_files/histones/Hs_Histones_with_annotated_SL_and_UTR_UCSCcoords.gtf',format='gtf')

# Extract CDS ranges and Downstream regions with SAME ORDER (important for overlap intersections)
names(histgtf) <- histgtf$gene
histgenes <- unique(histgtf$gene)
cdsgr <- subset(histgtf, type == 'CDS')[histgenes]
downgr <- subset(histgtf, type == 'DownsReg')[histgenes]

# Senescence dataset from the lab
myfiles <- list.files('STAR_aligns/',pattern='.*sorted.bam$',full.names=T)
fnames <- lapply(myfiles, basename) |> 
  gsub(pattern = '_Aligned.*$', replacement = '')
names(myfiles) <- fnames

# Execute function
res <- lapply(myfiles, process_file, cdsgr, downgr) |> 
  dplyr::bind_rows(.id = 'ID')

# Save results
write.table(res, 'results/histones/Histone_polyadenilation_table_results.tsv', sep = '\t')

# Load results
res <- read.delim('results/histones/Histone_polyadenilation_table_results.tsv')

# Load sample info
smpinf <- readxl::read_xlsx('sample_info.xlsx')
smpinf$ID <- paste0('Tube_', smpinf$ID)
smpinf <- smpinf[c(1:9,19:27),c('ID','Cells','Passage')]
smpinf$Passage <- factor(rep(c('early','mid','late'),each=3, times = 2), levels = c('early','mid','late'))

# Left join sample info to results
res <- dplyr::left_join(res, smpinf, by = 'ID')

# Apply some filters of expression and some filters of readthrough (if there is no readthrough when you add a count you will be only analyzing expression, since you normalize readthrough by it)
res <- dplyr::group_by(res, gene, Cells) |> 
  dplyr::mutate(
    mean_cds = mean(mate2incds),
    min_cds = min(mate2incds),
    min_cds2down = sum(cds2down > 0))

# First I will be filtering those genes with a minimum expression and also a minmum evidence of readthrough (if there is no cds2down reads in any sample, you polyA measures is just a measure of gene expression changes)
res_filt <- subset(res, mean_cds > 5 & min_cds2down > 2) # at least evidence of polyA in 3 samples

# After inspecting which genes had nocds2down reads in some samples; I could see that it was always mid or late passage samples with very few cds reads. Hence I will not exclude them from the analysis

# Inspect those genes were there are very low cds to downstream reads (direct evidence of readthorugh) but higher only downstream reads (indirect evidence of readthrough). As it is not the same to have 100 onlydown reads in a downstream region of 200 bp or 15000 bp, I include down region length info
# downsize <- as.data.frame(lengths(downgr)) |> 
#   tibble::rownames_to_column(var = 'gene')
# res_filt<- merge(res_filt, downsize, by = 'gene')
# res_filt$ratio <- res_filt$onlydown/res_filt$cds2down

# The majority of genes that have high only down reads and low cds to down reads seem due to reads that very likely do not proceed from aberrant readthrough. Hence I will use cds to down reads as they provide direct evidence of readthorough and in theory they should be proportional to those onlydown


# Calculate the polyA ratio
res$pAratio <- (res$cds2down + 1)/res$mate2incds
res_filt$pAratio <- (res_filt$cds2down + 1)/res_filt$mate2incds

# Plot general filtered results
ggplot(res_filt, aes(Passage, pAratio, colour = Passage, group = ID)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_point(shape = 16, position = position_jitterdodge()) +
  ylab('Readthrough ratio') +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        strip.text = element_text(size = 20)
  ) +
  facet_wrap(~ Cells)


#///////////////////////////////////////////////////////////////////////////////
###################### Histone Readthrough and Expression #####################
#///////////////////////////////////////////////////////////////////////////////

# To see if there is a relationship between aberrant readthorugh and expression I will cross readthrough and DE data

# First I will compute the LFC of readthrough ratio between early and late passages
LFCpA <- res_filt |> 
  subset(Passage != 'mid') |> 
  dplyr::group_by(gene, Cells) |> 
  dplyr::reframe(LFCpA = log2(mean(pAratio[Passage == 'late'])/mean(pAratio[Passage == 'early']))) |> 
  dplyr::ungroup()

# Plot LFCpA
ggplot(LFCpA, aes(Cells, LFCpA, colour = Cells)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, position = position_jitter(width = 0.2), size = 2.5) +
  theme_classic() +
  ylab('LFC readthrough ratio') +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.position = 'none')

# Load DESeq results
DEres <- readRDS('R_files/transposons/subfamily/res_indep.rds')
DEres <- DEres |> 
  lapply(dplyr::as_tibble) |> 
  dplyr::bind_rows(.id = 'Cells')

# Extract histone results
DEh <- DEres |>
  dplyr::semi_join(
    LFCpA,
    by = c("gene_name" = "gene", "Cells" = "Cells")
  ) |> 
  dplyr::rename(gene = gene_name)

# Merge both
DEpA <- merge(DEh, LFCpA, by = c('gene', 'Cells'))

# Correlation between DE LFC and LFCpA
ggplot(DEpA, aes(LFCpA, log2FoldChange)) +
  geom_point(shape = 16,  size = 2.5) +
  geom_smooth(method = 'lm') +
  ggpmisc::stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
                        formula = x ~ y,
                        method = 'lm',
                        parse = T,
                        size = 6
  ) +
  ylab('LFC expression') +
  xlab('LFC readthrough') +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20)) +
  facet_wrap(~ Cells, scales = 'free')



# Test if histones with no evidence of readthrough are less downregulated than those with evidence of readthrough
DEres <- dplyr::rename(DEres, gene = gene_name)
DEh_nofilt <- subset(DEres, gene %in% unique(res$gene))
DEh_nofilt <- DEh_nofilt |> 
  dplyr::mutate(Rdthr = ifelse(gene %in% DEh$gene, 'Evidence', 'No Evidence'))
ggplot(DEh_nofilt, aes(Rdthr, log2FoldChange, colour = Cells)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_point(shape = 16, position = position_jitterdodge()) +
  ylab('LFC expression') +
  xlab('Readthrough') +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25)) 

# It seems it is not the case, so this histones were there is no evidence of readthrough could be downregulated due to other reasons or have very low expression in senescence samples so readthrough cannot be technically detected. I will test this possibility
res <- res |> 
  dplyr::left_join(DEh_nofilt, by = c('gene', 'Cells'))
res_ev <- na.omit(res) |> 
  subset(Passage == 'late')
ggplot(res_ev, aes(Rdthr, log2(mate2incds), colour = Cells)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_point(shape = 16, position = position_jitterdodge()) +
  ylab('log2(Reads late Passage)') +
  xlab('Readthrough') +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25)) 


#///////////////////////////////////////////////////////////////////////////////
###################### Histone Readthrough and Splicing #####################
#///////////////////////////////////////////////////////////////////////////////

# We also saw that some of the histones not only displayed aberrant readthrough in senescent cells but also an alternative splicing after the stem loop. Hence I will see if this alternative splicing is somehow related with readthrough level and with histone expression

# First of all I will extract the splice junctions using the SJ files generated by STAR
SJf <- list.files('STAR_aligns/', pattern = '*._SJ.out.tab', full.names = T)
names(SJf) <- gsub('_SJ.*$', '', basename(SJf))
SJ <- lapply(SJf, read.delim, header = F)
SJgr <- lapply(SJ, \(sj) GRanges(seqnames = sj$V1, ranges = IRanges(sj$V2, sj$V3), strand = c("*", "+", "-")[sj$V4 + 1]))

# Extract de splice junctions from the cds2down reads
SJgr_cds2down <- lapply(myfiles, function(bamfile){
  # Extend regions to get complete pairs, also those were read 1 is very far way
  extgr <- reduce(resize(c(cdsgr,downgr),fix='center',width=2E6), ignore.strand= TRUE)
  rds <- load_reads(bamfile, extgr) 
  # Pairs where mate2 overlaps cds and the pair range span downstream of cds, surely demonstrate readthrough
  m2incds <- findOverlaps(second(rds),cdsgr) # mate 2 in cds
  overDown <- findOverlaps(granges(rds),resize(downgr, 500, fix = 'start')) # Here I don't  need to use the annotation, since forcing mate2 to be in cds ensures that the read comes from the histone and not other gene
  cdstodown <- intersect(m2incds,overDown)
  cds2down_rds <- rds[from(cdstodown)]
  names(cds2down_rds) <- names(cdsgr[to(cdstodown)])
  SJ <- unlist(junctions(cds2down_rds))
})

# Intersect the cds2down splice junctions with the ones detected by STAR
SJ_ovlap <- mapply(findOverlaps, SJgr_cds2down, SJgr, MoreArgs = list(type = 'equal'))

# !! Not all junctions of cds2down reads can be found in STAR SJ
lapply(SJgr_cds2down, length) |> dplyr::bind_rows(.id = 'sample') == lapply(SJ_ovlap, length) |> dplyr::bind_rows(.id = 'sample')
rbind(lapply(SJgr_cds2down, length) |> dplyr::bind_rows(.id = 'sample'), lapply(SJ_ovlap, length) |> dplyr::bind_rows(.id = 'sample'))

# Inspect what are the ones not present in STAR SJ
SJgr_cds2down$Tube_23[!(1:length(SJgr_cds2down$Tube_23) %in% from(SJ_ovlap$Tube_23))]

# Select the STAR SJ present in the cds2down reads
SJ_cds2down <- mapply(\(x,y) x[to(y),], SJ, SJ_ovlap, SIMPLIFY = F) |> 
  dplyr::bind_rows(.id = 'ID')
colnames(SJ_cds2down) <- c(
  "ID",
  "chr",
  "start",
  "end",
  "strand",
  "motif",
  "annotated",
  "unique_reads",
  "multi_reads",
  "max_overhang"
)

# Change annotation and strand columns
SJ_cds2down <- SJ_cds2down |> 
  dplyr::mutate(annotated = dplyr::case_when(annotated == 0 ~ 'unannotated',
                                             annotated == 1 ~ 'annotated'),
                strand = dplyr::case_when(strand == 0 ~ '*',
                                          strand == 1 ~ '+',
                                          strand == 2 ~ '-'))

# Asign histone genes to STAR splice junctions
SJgr_cds2down_filt <- mapply(\(x,y) x[from(y)], SJgr_cds2down, SJ_ovlap, SIMPLIFY = F) |> 
  GRangesList() |> 
  unlist(,use.names = F)
identical(start(SJgr_cds2down_filt), SJ_cds2down$start) # Test that Granges filtered object has the same order of the SJ STAR file
SJ_cds2down$gene <- names(SJgr_cds2down_filt)

# Eliminate duplicate STAR junctions, as we overlap reads with recognized splice junctions, one junction will be supported by several reads, generating duplicates, hence I will take only the unique values
SJ_cds2down <- unique(SJ_cds2down)

# Fuse this data with metadata
SJ_cds2down <- merge(SJ_cds2down, smpinf, by = 'ID')

# Complete data with samples where splice junction were not annotated adding the junction info but with 0 counts
# Extract all detected junctions
junction_info <- SJ_cds2down |>
  dplyr::ungroup() |>
  dplyr::distinct(
    chr, start, end, strand,
    motif, annotated, gene
  )
# Generate a combination of all detected junctions with all samples
complete_grid <- tidyr::crossing(
  junction_info,
  smpinf
)
# Add the samples were splice junction was not detected as rows with 0 counts
SJ_complete <- complete_grid |>
  dplyr::left_join(
    SJ_cds2down |>
      dplyr::ungroup() |>
      dplyr::select(
        ID, Cells,
        chr, start, end, strand,
        unique_reads, multi_reads, max_overhang
      ),
    by = c(
      "ID", "Cells",
      "chr", "start", "end", "strand"
    )
  ) |>
  dplyr::mutate(
    dplyr::across(
      c(unique_reads, multi_reads, max_overhang),
      ~ tidyr::replace_na(.x, 0L)
    )
  )

# Test I have done it correctly
# Test I have not changed original data
cols <- sort(names(SJ_cds2down))

df1 <- SJ_cds2down |>
  dplyr::ungroup() |>
  dplyr::select(dplyr::all_of(cols)) |>
  dplyr::arrange(dplyr::across(dplyr::everything()))

df2 <- SJ_complete |>
  dplyr::filter(unique_reads > 0 | multi_reads > 0) |>
  dplyr::select(dplyr::all_of(cols)) |>
  dplyr::arrange(dplyr::across(dplyr::everything()))

identical(as.data.frame(df1), as.data.frame(df2))
# Test that zero rows are all different
nrow(subset(SJ_complete, unique_reads == 0 & multi_reads == 0)) == nrow(dplyr::distinct(subset(SJ_complete, unique_reads == 0 & multi_reads == 0), dplyr::across(dplyr::all_of(cols))))

# Explore in how many samples is present each splice junction
SJ_complete <- SJ_complete |> 
  dplyr::group_by(Cells, start, end, strand) |> 
  dplyr::mutate(Nsamp = sum(unique_reads > 0)) |> 
  dplyr::ungroup()

# Filter out splice junctions that are not present in any sample of a cell line 
SJ_complete <- subset(SJ_complete, Nsamp > 0)

# I compute the splicing ratio to normalize the splicing by gene expression
# FIrst I add the expression data for each sample
SJ_complete <- dplyr::left_join(SJ_complete, res_filt[,c(1:2,4)], by = c('ID','gene')) |> 
  na.omit()
SJ_complete <- SJ_complete |> 
  dplyr::group_by(chr, start, end, strand, ID, Cells) |> 
  dplyr::mutate(SJratio = (unique_reads + 1)/mate2incds)

# Plot SJ ratio by passage
# Those splice junctions present in all samples
ggplot(SJ_complete |> subset(Nsamp == 9), aes(Passage, SJratio, colour = Passage)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_point(shape = 16, position = position_jitterdodge()) +
  ylab('Splice Ratio') +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25)) +
  facet_wrap(~Cells)

# I compute the LFC of late vs early for each splice junction
LFC_SJ <- SJ_complete |> 
  dplyr::group_by(Cells, start, end, strand) |> 
  dplyr::summarise(LFC_SJ = log2(mean(SJratio[Passage == 'late'])/mean(SJratio[Passage == 'early'])),
                   Nsamp = unique(Nsamp),
                   mean_rds = mean(unique_reads),
                   gene = unique(gene),
                   annotated = unique(annotated),
                   .groups = "drop")


# I will start with an anlysis of all splice junctions to see if they are correlated with readthrough (they should be as they are after the stem loop) and with expression changes. But first I will check if there are differences in number of samples and LFC of splice junction depending on if they are or not annotated

# As for correlating expression readthrough and splicing I must have one splice junction per gene, for those with more than one splice junction I will select the one with the greatest mean of unique reads, if there are ties the one present in more smaples and if still peak arbitraretely the first
LFC_SJ_filt <- LFC_SJ |> 
  dplyr::group_by(Cells, gene) |> 
  dplyr::filter(mean_rds == max(mean_rds, na.rm = TRUE)) |>
  dplyr::filter(Nsamp == max(Nsamp, na.rm = TRUE)) |>
  dplyr::slice_head(n = 1) |> # if still ties peak the first
  dplyr::ungroup()

# Merge to SJ data the expression and readthrough data
SJpAexpr <- dplyr::left_join(LFC_SJ_filt, DEres, by = c('gene', 'Cells')) |> 
  dplyr::left_join(LFCpA, by =c('gene', 'Cells'))

# Correlation splicing and readthrough
ggplot(SJpAexpr, aes(LFCpA, LFC_SJ)) +
  geom_point(shape = 16, size = 3) +
  geom_smooth(method = 'lm') +
  ggpmisc::stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
                        formula = x ~ y,
                        method = 'lm',
                        parse = T,
                        size = 6) +
  ylab('LFC Splicing') +
  xlab('LFC Readthrough') +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20)) +
  facet_wrap(~Cells, scales = 'free')

# As expected there is a correlation between readthrough LFC and splicing LFC.
# I will divide LFCpA data into histones with splicing detected and without and check the readthroufh LFC
LFCpA <- dplyr::mutate(LFCpA, splicing = ifelse(gene %in% LFC_SJ_filt$gene, 'Evidence', 'No Evidence'))
ggplot(LFCpA, aes(splicing, LFCpA, colour = Cells)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, position = position_jitterdodge()) +
  theme_classic() +
  ylab('LFC Readthrough') +
  xlab('Splicing') +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))

# Correlation splicing and expression
ggplot(SJpAexpr, aes(log2FoldChange, LFC_SJ)) +
  geom_point(shape = 16, size = 3) +
  geom_smooth(method = 'lm') +
  ggpmisc::stat_poly_eq(aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
                        formula = x ~ y,
                        method = 'lm',
                        parse = T,
                        size = 6) +
  ylab('LFC Splicing') +
  xlab('LFC Expression') +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20)) +
  facet_wrap(~ Cells, scales = 'free')
# It seems that those histones with higher splicing in senescent cells are the ones that suffer less downregulation
# Check if in histones with evidence of readthrough the ones with splicing are less downregulated
DEpA <- merge(DEh, LFCpA, by = c('gene', 'Cells'))
DEpA <- dplyr::mutate(DEpA, splicing = ifelse(gene %in% LFC_SJ_filt$gene, 'Evidence', 'No Evidence'))
ggplot(DEpA, aes(splicing, log2FoldChange, colour = Cells)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, position = position_jitterdodge()) +
  theme_classic() +
  ylab('LFC Expression') +
  xlab('Splicing') +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))

# It seems that histones with evidence of splicing are a little bit less downregulated, but lest check those where splicing is detected in all samples
DEpA <- dplyr::mutate(DEpA, splicing = ifelse(gene %in% LFC_SJ_filt[LFC_SJ_filt$Nsamp == 9,]$gene, 'Evidence', 'No Evidence'))
ggplot(DEpA, aes(splicing, log2FoldChange, colour = Cells)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, position = position_jitterdodge()) +
  theme_classic() +
  ylab('LFC Expression') +
  xlab('Splicing') +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))

# In this case it is clear that they are clearly less downregulated, reinforcing the idea that splicing could somehow avoid the degradation of the mRNA


#///////////////////////////////////////////////////////////////////////////////
###################### Histone Readthrough Regulators #####################
#///////////////////////////////////////////////////////////////////////////////

# I will look for possible candidates that could explain histone aberrant readthrough

# First of all I generate gene sets that summarise genes regulating histone mRNA processing

# Get GO gene sets
go <- msigdbr(
  species = "Homo sapiens",
  collection = "C5"
) |> 
  dplyr::filter(gs_subcollection %in% c("GO:BP", "GO:CC", "GO:MF"))

# Terms we want, the core terms strictly related to histone mRNA processing and the extended related also with other histone processes
core_GO <- c(
  "GO:0008334",  # histone mRNA metabolic process
  "GO:0035363"   # histone locus body
)

extended_GO <- c(
  core_GO,
  "GO:0006335",  # DNA replication-dependent chromatin assembly
  "GO:0140713"   # histone chaperone activity
)

# Make gene lists
HISTONE_MRNA_CORE <- go |> 
  dplyr::filter(gs_exact_source %in% core_GO) |> 
  dplyr::pull(gene_symbol) |> 
  unique()

HISTONE_MRNA_EXTENDED <- go |> 
  dplyr::filter(gs_exact_source %in% extended_GO) |> 
  dplyr::pull(gene_symbol) |> 
  unique()

# Number of genes
length(HISTONE_MRNA_CORE)
length(HISTONE_MRNA_EXTENDED)

# Check expression of these genes

# Load DESeq experiment
dds <- readRDS('R_files/transposons/subfamily/dds_exp_indep.rds')

# Extract stabilized varianze counts
vsd <- lapply(dds, DESeq2::vst)
cts <- lapply(vsd, assay) |> 
  lapply(as.data.frame) |> 
  lapply(tibble::rownames_to_column, var = 'ensmbleID')
cts <- dplyr::full_join(cts$GM05565, cts$GM00038, by = 'ensmbleID')

# Associate gene names to ensemble ID
gnames <- lapply(dds, rowData) |> 
  lapply(as.data.frame) |> 
  lapply(tibble::rownames_to_column, var = 'ensmbleID') |> 
  lapply(dplyr::select, ensmbleID, gene_name)
gnames <- dplyr::full_join(gnames$GM05565, gnames$GM00038)

cts <- merge(cts, gnames, by = 'ensmbleID')

# Filter out genes from datasets
histcore <- subset(cts, gene_name %in% HISTONE_MRNA_CORE) |> 
  dplyr::select(-ensmbleID)
rownames(histcore) <- histcore$gene_name
histcore$gene_name <- NULL
histext <- subset(cts, gene_name %in% HISTONE_MRNA_EXTENDED) |> 
  dplyr::select(-ensmbleID)
rownames(histext) <- histext$gene_name
histext$gene_name <- NULL
# Plot a heatmap of the genes
ann <- rbind(
  as.data.frame(colData(dds$GM05565)),
  as.data.frame(colData(dds$GM00038))
)
# Heatmap of core regulators
pheatmap::pheatmap(histcore, scale = 'row',
                   color = colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(25),
                   annotation_col = ann[,1:2],
                   fontsize = 10,
                   annotation_names_col = F)

# Heatmap of extended regulators
pheatmap::pheatmap(histext, scale = 'row',
                   color = colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(25),
                   annotation_col = ann[,1:2],
                   fontsize = 10,
                   annotation_names_col = F)

# Extended regulators cluster late and mid passage samples separated from early samples
