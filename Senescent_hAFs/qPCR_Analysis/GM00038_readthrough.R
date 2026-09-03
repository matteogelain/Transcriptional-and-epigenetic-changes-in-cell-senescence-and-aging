# libraries
require('dplyr')
require('data.table')
require('ggplot2')
require('ggh4x')
require('readr')

#20260722 Histone analysis GM00038 senescent vs proliferative  aberrant polyadenilation#

# Read file, and modify column containing senescent and proliferative, by eliminating
# the number of the sample
data <- fread("20260722_Cq_GM38_Hist_Prol_Sen_Aberrant_H3.csv") |>
  mutate(Condition = ifelse(grepl("Prol", Sample), "Prol", "Sen"), Cq = as.numeric(Cq))

dt_path <- '20260722_Cq_GM38_Hist_Prol_Sen_Aberrant_H3.csv'
control <- 'Prol'
pattern <- "GM38 | [0-9]+"
filt_outlier <- F
out_path <- ''

# Load results
data <- fread(dt_path) |>
  mutate(Condition = gsub(pattern, '', Sample))
data$Cq <- as.numeric(data$Cq)

# Define the endogenous control for each target
endogen_map <- tibble(Target = c("Aberrant H3C2","Aberrant H3C8","Aberrant H3C11"),endogen = c("H3C2","H3C8","H3C11"))

# Eliminate outliers if wanted
if (filt_outlier == T) {
  clean_data <- qpcr_clean(data,cq = Cq,threshold = threshold,Target,Sample,Condition)
  
  # Check eliminated values
  outliers <- data[
    !data$`Well Position` %in% clean_data$`Well Position`,
  ]
  
  outlier_context <- data |>
    mutate(outlier = ifelse(data$`Well Position` %in% outliers$`Well Position`, T, F))
  
  write.csv(outlier_context,
    paste0(out_path, gsub('_.*', '', basename(dt_path)),'_qPCR_outliers.csv'))
  
  # Calculate mean Cq for each Sample and Target
  data_noNA <- clean_data |>
    filter(!is.na(Cq)) |>
    group_by(Sample, Target) |>
    summarise(Condition = unique(Condition),Cq = mean(Cq),.groups = "drop")
  
} else {
  
  # Calculate median Cq to reduce the effect of outliers
  data_noNA <- data |>
    filter(!is.na(Cq)) |> 
    group_by(Sample, Target) |>
    summarise(Condition = unique(Condition),median_Cq = median(Cq),.groups = "drop")
}

# For the outlier-filtered data, rename the Cq column, so that both branches have the same structure
if (filt_outlier == T) {
  data_noNA <- data_noNA |>
    rename(median_Cq = Cq)
}


# Extract the endogenous controls
endogen_data <- data_noNA |>
  filter(Target %in% endogen_map$endogen) |>
  select(Sample, Target,median_Cq) |>
  rename(endogen = Target,
    endogen_Cq = median_Cq)

# Get target data and associate each target with its endogenous control
dCq_data <- data_noNA |>
  filter(Target %in% endogen_map$Target) |>
  left_join(endogen_map, by = "Target") |>
  left_join(endogen_data,by = c("Sample", "endogen")) |>
  mutate(dCq = median_Cq - endogen_Cq)

# Set control condition and calculate delta-delta Cq
control_data <- dCq_data |>
  filter(Condition == control) |>
  group_by(Target) |> 
  summarise(dCq_control = mean(dCq),.groups = "drop")

# Calculate ddCq and RQ
RQ_data <- dCq_data |>
  left_join(control_data,by = "Target") |>
  group_by(Target, Condition) |> 
  mutate(ddCq = dCq - dCq_control, RQ = 2^(-ddCq),mean_RQ = 2^(-mean(ddCq))) |> 
  ungroup()


# Write results
write_rds(RQ_data, paste0(
    out_path,
    gsub(
      '.csv',
      '_multiple_endogen.RDS',
      basename(dt_path)
)))

p<-ggplot(RQ_data, aes(Condition, RQ, fill = Condition)) +
  stat_summary(position = position_dodge(width = 0.9), fun = mean, geom = 'bar') +
  stat_summary(position = position_dodge(width = 0.9),fun.data = mean_se, geom = 'errorbar',width = 0.2) +
  geom_point(position = position_jitter(width = 0.1),alpha = 0.8) +
  theme_classic() +
  facet_wrap(~Target, scales = 'free_y')

ggsave('20260727_Human_Prolif_vs_Senesc_Aberrant_H3_Normalized_GM38.pdf', p, width = 10, height = 5)

# Statistical analysis
res_GAPDH <- readRDS("20260722_Cq_GM38_Hist_Prol_Sen_Aberrant_H3_multiple_endogen.RDS")

# check normality, only Mm_SINE1 is not normal, the other ones don't follow
# a normal distribution
shap_t <- by(res_GAPDH$dCq, res_GAPDH$Target, shapiro.test)

# Check variances
F_t <- by(res_GAPDH, res_GAPDH$Target, \(df) var.test(dCq ~ Condition, data = df))

# T-test
T_t <- by(res_GAPDH, res_GAPDH$Target, \(df) t.test(dCq ~ Condition, data = df, var.equal = TRUE))
pvals <- data.frame(
  Target = names(T_t),
  p = sapply(T_t, \(x) x$p.value)
)

# plot the results
yplot <- res_GAPDH |> 
  group_by(Target) |> 
  summarise(y = max(RQ, na.rm = TRUE) * 1.1)

pvals <- left_join(pvals, yplot, by = "Target")
pvals$signif <- dplyr::case_when(
  pvals$p < 0.001 ~ "***",
  pvals$p < 0.01  ~ "**",
  pvals$p < 0.05  ~ "*",
  TRUE            ~ "ns"
)
