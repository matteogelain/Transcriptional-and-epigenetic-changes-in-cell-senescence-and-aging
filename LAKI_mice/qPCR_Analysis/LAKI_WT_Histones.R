# libraries
require('dplyr')
require('data.table')
require('ggplot2')

### 20260421 LAKI histones liver ####

# Separate tissues: liver and kidney
res <- fread('20260421_LAKI_liver_kidney_hist_Matteo_Results_20260422_132441.csv') |>  
  mutate(Tissue = gsub('.* ', '', Sample))
res_liv <- filter(res, Tissue == 'LIVER')

write.csv(res_liv, '20260421_LAKI_liver_hist.csv')

# Liver 
dt_path <- '20260421_LAKI_liver_hist.csv'
endogen <- '18S'
control <- 'WT LIVER'
pattern <- '[0-9]+ '
filt_outlier <- F
out_path <- ''

source('C:/Users/Lenovo/Desktop/Erasmus/Results/qPCR results/qPCR_R_Analysis/qPCR_analysis.R')

p<-ggplot(RQ_data, aes(Condition, RQ, fill = Condition)) +
  stat_summary(position = position_dodge(width = 0.9), fun = mean, geom = 'bar') +
  stat_summary(position = position_dodge(width = 0.9),fun.data = mean_se, geom = 'errorbar',width = 0.2) +
  geom_point(position = position_jitter(width = 0.1),alpha = 0.8) +
  theme_classic() +
  facet_wrap(~Target, scales = 'free_y')

# Statistical analysis
# check normality
res_18S <- readRDS("20260421_LAKI_liver_hist_18S.RDS")
shap_t <- by(res_18S$dCq, res_18S$Target, shapiro.test)

# Check the variances
F_t <- by(res_18S, res_18S$Target, \(df) var.test(dCq ~ Condition, data = df))

# T-test
T_t <- by(res_18S, res_18S$Target, \(df) t.test(dCq ~ Condition, data = df, var.equal = TRUE))
pvals <- data.frame(Target = names(T_t), p = sapply(T_t, \(x) x$p.value))

# plot the results
yplot <- res_18S |> 
  group_by(Target) |> 
  summarise(y = max(RQ, na.rm = TRUE) * 1.1)

pvals <- left_join(pvals, yplot, by = "Target")

pvals$signif <- dplyr::case_when(
  pvals$p < 0.001 ~ "***",
  pvals$p < 0.01  ~ "**",
  pvals$p < 0.05  ~ "*",
  TRUE            ~ "ns"
)
