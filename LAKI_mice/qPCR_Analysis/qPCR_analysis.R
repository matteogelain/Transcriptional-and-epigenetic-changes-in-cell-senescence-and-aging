# libraries
require('dplyr')
require('data.table')
require('readr')

#load results
  data <- fread(dt_path) |>
    mutate(Condition = gsub(pattern, '', Sample))
  data$Cq <- as.numeric(data$Cq)

#eliminate outliers if wanted
if(filt_outlier == T){
  clean_data <- qpcr_clean(data, cq = Cq, threshold = threshold, Target, Sample, Condition)
  #check the eliminated values
  outliers <- data[!data$`Well Position` %in% clean_data$`Well Position`,]
  outlier_context <- data |>
    mutate(outlier = ifelse(data$`Well Position` %in% outliers$`Well Position`, T, F))
  
  write.csv(outlier_context, paste0(out_path,gsub('_.*','',basename(dt_path)),'_qPCR_outliers.csv'))
  
  #calculate mean Cq
  clean_data <- clean_data |>
    group_by(Sample, Target) |>
    mutate(mean_Cq = mean(Cq))
  
  # Calculate mean_Cq for the endogen control per condition and sample
  endogen_data <- clean_data |>
    filter(Target == endogen) |>
    group_by(Sample) |>
    summarise(endogen_Cq = mean(mean_Cq, na.rm = TRUE),.groups = "drop")
} else {

  #Calculate median Cq, to avoid outliers effect
  data_noNA <- data |>
    filter(!is.na(Cq)) |> 
    group_by(Sample, Target) |>
    reframe(Condition = Condition, Cq = Cq, median_Cq = median(Cq))
  
  #Calculate median_Cq for the endogen control per condition and sample
  endogen_data <- data_noNA |>
    filter(Target == endogen) |>
    group_by(Sample) |>
    summarise(endogen_Cq = unique(median_Cq, na.rm = TRUE), .groups = "drop")
}

#Get the target data (non-endogen) and join with the endogen
dCq_data <- data_noNA |>
  filter(Target != endogen) |>
  group_by(Sample, Target) |>
  summarise(median_Cq = unique(median_Cq, na.rm = TRUE),
    Condition = unique(Condition),.groups = "drop") |> 
  left_join(endogen_data, by = "Sample") |> 
#Calculate deltaCq for each target measurement.
  mutate(dCq = median_Cq - endogen_Cq)

#Set control condition and calculate deltadeltaCq
control_data <- dCq_data |>
  filter(Condition == control) |>
  group_by(Target) |> 
  summarise(dCq_control = mean(dCq))

#join deltaCq data of targets with control data and calculate deltadeltaCq and RQ
RQ_data <- dCq_data |>
  left_join(control_data, by = 'Target') |>
  group_by(Target, Condition) |> 
  mutate(ddCq = dCq - dCq_control, RQ = 2^(-ddCq), mean_RQ = 2^(-mean(ddCq))) |> 
  ungroup()

#write results to an object
write_rds(RQ_data,paste0(out_path,gsub('.csv',paste0('_', endogen, '.RDS'),basename(dt_path))))
