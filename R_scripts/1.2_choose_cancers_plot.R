# choose cancers - plot ---------------------------------------------------

# build the data frame
C_df_all_cancers <- data.frame(matrix(ncol=length(cancers_all), 
                                      nrow = n_rep * K_folds))
colnames(C_df_all_cancers) <- cancers_all

for(cancer in cancers_all){
  
  load(file = paste0("data_fit/", cancer, "/ridge/pred_clin_mRNA.RData"))
  
  C_df_all_cancers[, cancer] <- C_df_final$mRNA
}

C_df_ggplot <- stack(C_df_all_cancers)

# p-values and colors
stars_p_val_BH <- my_stars.pval_vect(p_val_test_BH)
col_p_val = col_cancer_name <- rep("black", length(cancers_all))
col_p_val[p_val_test_BH <= 0.01] = col_cancer_name[p_val_test_BH <= 0.01] <- "red"
names(col_p_val) = names(col_cancer_name) <- cancers_all

# plot
C_all_cancers_plot <- ggplot() + 
  geom_boxplot(data = C_df_ggplot, aes(x = ind,# factor(ind, levels = names(sort(med_C_sorted, decreasing = T))), 
                                       y = values), position = position_dodge(0.75), 
               fill = "#386cb0") + 
  xlab("Cancer") + ylab("C-index") + ylim(NA, 1) +
  theme_Publication() + 
  geom_abline(slope = 0, intercept = 0.6, col = "red", lwd = 1, lty = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                   colour = col_cancer_name, 
                                   size = 12)) +
  geom_text(aes(x = 1:length(cancers_all), y = rep(1, length(cancers_all)), 
                label = stars_p_val_BH), 
            col = col_p_val[names(med_C_sorted)], size = 4)
print(C_all_cancers_plot)

ggsave(C_all_cancers_plot, filename = "pdf/C_ridge_all_cancers.pdf")
