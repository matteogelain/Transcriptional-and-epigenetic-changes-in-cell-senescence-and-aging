setwd('C:/Users/Lenovo/Desktop/Erasmus/Results/qPCR results/qPCR_R_Analysis')
# libraries
require('dplyr')
require('data.table')
require('ggplot2')
require('ggh4x')

######### 20260619 Histone analysis GM00038 senescent vs proliferative #########

# Read file, and modify column containing senescent and proliferative, by eliminating
# the number of the sample
data <- fread("20260619_Human_Prolif_vs_Senesc_Cq.csv") |>
  mutate(Condition = gsub(" [0-9]", "", Sample))
data$Cq <- as.numeric(data$Cq)

# Variables creation
dt_path <- '20260619_Human_Prolif_vs_Senesc_Cq.csv'
endogen <- 'GAPDH'
control <- 'Proliferative'
pattern <- ' +[0-9]'
filt_outlier <- F
out_path <- ''

source('qPCR_analysis.R')
p<-ggplot(RQ_data, aes(x = Target, y = RQ, fill = Condition)) +
  stat_summary(position = position_dodge(width = 0.9), fun = mean, geom = 'bar') +
  stat_summary(position = position_dodge(width = 0.9), fun.data = mean_se, geom = 'errorbar', width = 0.2) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
  theme_classic()

ggsave('20260619_Human_Prolif_vs_Senesc.pdf', p, width = 10, height = 5)

# Statistical analysis
res_GAPDH <- readRDS("20260619_Human_Prolif_vs_Senesc_Cq_GAPDH.RDS")

# check normality, they all follow the normal distribution
shap_t <- by(res_GAPDH$dCq, res_GAPDH$Target, shapiro.test)

# check variances, they have all the same variances, except for H3C11
F_t <- by(res_GAPDH, res_GAPDH$Target, \(df) var.test(dCq ~ Condition, data = df))

# Compare means using t-test, assuming equal variances
T_t <- by(res_GAPDH, res_GAPDH$Target, \(df) t.test(dCq ~ Condition, data = df, var.equal = TRUE))
pvals <- data.frame(Target = names(T_t),p = sapply(T_t, \(x) x$p.value))

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
