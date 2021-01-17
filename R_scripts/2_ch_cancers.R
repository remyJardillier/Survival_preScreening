# data frame with all the characteristics
clinical_feat <- c("Name", "n patients", "Censoring rate",
                   "Survival - 3 years", "C-index")
data_ch <- data.frame(matrix(ncol = length(clinical_feat), nrow = length(cancers_all)))
row.names(data_ch) <- cancers_all
colnames(data_ch) <- clinical_feat

for(cancer in cancers_all){
  
  print(paste("*** start learning for", cancer, "***"))
  
  # load the data ---
  source(file = "load_data/load_data_final.R")
  
  # number of patients
  data_ch[cancer, "n patients"] <- nrow(clin)

  # censoring rate
  data_ch[cancer, "Censoring rate"] <- round(sum(clin$status == 0)/nrow(clin), digits = 2)

  # survival at 3 years
  km_fit <- survfit(y_cox~1)
  km_fit_summary <- summary(km_fit, times = 3)
  data_ch[cancer, "Survival - 3 years"] <- signif(km_fit_summary$surv[1],2)
  
  # C-index médian
  load(file = paste0("data_fit/", cancer, "/ridge/pred_clin_mRNA.RData"))
  data_ch[cancer, "C-index"] <- signif(median(C_df_final$mRNA, na.rm = T),2)
}

print(data_ch )

write.xlsx(data_ch, file = "tables/data_ch.xlsx", row.names = T)
