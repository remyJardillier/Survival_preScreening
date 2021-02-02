cancers_used <- names(p_val_test_BH[p_val_test_BH <= s_level])
cancers_not_used <- names(p_val_test_BH[p_val_test_BH > s_level])

# Clinical data available -------------------------------------------------

print("*** Checking which clinical data is available for each cancer ***")

# data frame with all the clinical features
clin_available_mRNA <- data.frame(matrix(ncol = 6, nrow = length(cancers_all)))
colnames(clin_available_mRNA) <- c("age", "gender", "grade", "T_stage", "N_stage", "M_stage")
row.names(clin_available_mRNA) <- cancers_all

# vector for the plot
clin_text_plot <- rep(NA, length(cancers_all))
names(clin_text_plot) <- cancers_all

for(cancer in cancers_all){
  
  print(paste("*** start learning for", cancer, "***"))
  
  # load the data ---
  source(file = "load_data/load_data_final.R")
  
  n_patients <- nrow(clin_data_cox)
  n_patients
  
  clin_available_mRNA[cancer,] <- as.numeric(c("age", "gender", "grade", 
                                               "T_stage", "N_stage", "M_stage") %in%
                                               colnames(clin_data_cox))
  
  text_tmp <- ""
  if("grade" %in% colnames(clin_data_cox)){
    text_tmp <- paste0(text_tmp, "G")
  }
  if("T_stage" %in% colnames(clin_data_cox)){
    text_tmp <- paste0(text_tmp, "T")
  }
  if("N_stage" %in% colnames(clin_data_cox)){
    text_tmp <- paste0(text_tmp, "N")
  }
  if("M_stage" %in% colnames(clin_data_cox)){
    text_tmp <- paste0(text_tmp, "M")
  }
  
  clin_text_plot[cancer] <- text_tmp
  
}

print(clin_available_mRNA)

write.xlsx(clin_available_mRNA, file = "tables/clin_available_mRNA.xlsx",
           row.names = T)
print("Table saved in: 'tables/clin_available_mRNA.xlsx'")

save(clin_text_plot, file = "clin_avalilable.RData")


# p-values - Mix > clinical ? ---------------------------------------------

p_val_vect_C = p_val_vect_IBS  <- rep(NA, length(cancers_all))
names(p_val_vect_C) = names(p_val_vect_IBS)  <- cancers_all

for(cancer in cancers_all){
  
  load(file = paste0("data_fit/", cancer, "/ridge/pred_clin_mRNA.RData"))
  
  p_val_vect_C[cancer] <- wilcox.test(C_df_final$both, C_df_final$clin, 
                                      alternative = "greater")$p.value
  p_val_vect_IBS[cancer] <- wilcox.test(IBS_df_final$both, IBS_df_final$clin, 
                                        alternative = "less")$p.value
}

p_val_vect_C_BH <- p.adjust(p_val_vect_C, method = "BH")
stars_C <- sapply(p_val_vect_C, my_stars.pval)
length(cancers_all)
sum(p_val_vect_C < 0.05)

p_val_vect_IBS_BH <- p.adjust(p_val_vect_IBS, method = "BH")
stars_IBS <- sapply(p_val_vect_IBS, my_stars.pval)
length(cancers_all)
sum(p_val_vect_IBS < 0.05)

# Boxplot - cancers not used -----------------------------------------------------------------

# build the data frame
C_list = IBS_list <- list()

for(cancer in cancers_not_used){
  
  load(file = paste0("data_fit/", cancer, "/ridge/pred_clin_mRNA.RData"))
  
  C_list[[cancer]] <- stack(C_df_final)
  IBS_list[[cancer]] <- stack(IBS_df_final)
}

C_final_ggplot <- melt(C_list)
IBS_final_ggplot <- melt(IBS_list)

load(file = "clin_avalilable.RData")

# ggplot - C-index
plot_clin_mRNA_C <- ggplot() + 
  geom_boxplot(data = C_final_ggplot, aes(x = factor(L1, levels = cancers_not_used), 
                                          y = value, fill = ind),
               position = position_dodge(0.75)) + 
  theme_Publication_legend_right() + xlab("Cancer") + ylab("C-index")+
  theme(legend.title=element_blank(), legend.text = element_text(size = 15)) + ylim(NA, 1) +
  geom_text(aes(x = 1:length(cancers_not_used), y = rep(1, length(cancers_not_used)),
                label = stars_C[cancers_not_used]), col = "darkorchid", size = 5) +
  scale_fill_manual(values=c("brown3","#386cb0", "orchid"),
                    labels = c("Clinical", "mRNA", "Mixed")) +
  geom_text(aes(x = 1:length(cancers_not_used), y = rep(0.15, length(cancers_not_used)), 
                label = clin_text_plot[cancers_not_used]), col = "brown3", size = 4) +
  #theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1)) +
  geom_bracket(xmin = length(cancers_not_used) - 0.25, 
               xmax = length(cancers_not_used) + 0.25, y.position = 0.92, 
               label = "", col = "orchid", size = 1)

plot_clin_mRNA_C


# ggplot - IBS
plot_clin_mRNA_IBS <- ggplot() + 
  geom_boxplot(data = IBS_final_ggplot, aes(x = factor(L1, levels = cancers_not_used), 
                                            y = value, fill = ind),
               position = position_dodge(0.75)) + 
  theme_Publication() + xlab("Cancer") + ylab("IBS")+ ylim(NA, 0.6) +
  theme(legend.title=element_blank(), legend.text = element_text(size = 15)) + 
  geom_text(aes(x = 1:length(cancers_not_used), y = rep(0.5, length(cancers_not_used)),
                label = stars_IBS[cancers_not_used]), col = "darkorchid", size = 5) +
  scale_fill_manual(values=c("brown3","#386cb0", "orchid"),
                    labels = c("Clinical", "mRNA", "Mixed")) +
  geom_text(aes(x = 1:length(cancers_not_used), y = rep(0, length(cancers_not_used)), 
                label = clin_text_plot[cancers_not_used]), col = "brown3", size = 4) +
  # theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1)) +
  geom_bracket(xmin = length(cancers_not_used) - 0.25, 
               xmax = length(cancers_not_used) + 0.25, y.position = 0.45, 
               label = "", col = "orchid", size = 1)

plot_clin_mRNA_IBS

plot_clin_mRNA <- ggarrange(plot_clin_mRNA_C, plot_clin_mRNA_IBS, nrow = 2, ncol = 1,
                            labels = c("A", "B"))
plot_clin_mRNA <- annotate_figure(plot_clin_mRNA, top = text_grob("Cancers not used"))
print(plot_clin_mRNA)

ggsave(plot_clin_mRNA, filename = "pdf/suppl_clin_mRNA_can_NotUsed.pdf")
print("Figure saved in 'pdf/suppl_clin_mRNA_can_NotUsed.pdf'")

# Boxplot - cancers used -----------------------------------------------------------------

# build the data frame
C_list = IBS_list <- list()

for(cancer in cancers_used){
  
  load(file = paste0("data_fit/", cancer, "/ridge/pred_clin_mRNA.RData"))
  
  C_list[[cancer]] <- stack(C_df_final)
  IBS_list[[cancer]] <- stack(IBS_df_final)
}

C_final_ggplot <- melt(C_list)
IBS_final_ggplot <- melt(IBS_list)

load(file = "clin_avalilable.RData")

# ggplot - C-index
plot_clin_mRNA_C <- ggplot() + 
  geom_boxplot(data = C_final_ggplot, aes(x = factor(L1, levels = cancers_used), 
                                          y = value, fill = ind),
               position = position_dodge(0.75)) + 
  theme_Publication_legend_right() + xlab("Cancer") + ylab("C-index")+
  theme(legend.title=element_blank(), legend.text = element_text(size = 15)) + ylim(NA, 1) +
  geom_text(aes(x = 1:length(cancers_used), y = rep(1, length(cancers_used)),
                label = stars_C[cancers_used]), col = "darkorchid", size = 5) +
  scale_fill_manual(values=c("brown3","#386cb0", "orchid"),
                    labels = c("Clinical", "mRNA", "Mixed")) +
  geom_text(aes(x = 1:length(cancers_used), y = rep(0.15, length(cancers_used)), 
                label = clin_text_plot[cancers_used]), col = "brown3", size = 4) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1)) +
  geom_bracket(xmin = length(cancers_used) - 0.25, 
               xmax = length(cancers_used) + 0.25, y.position = 0.92, 
               label = "", col = "orchid", size = 1)

plot_clin_mRNA_C


# ggplot - IBS
plot_clin_mRNA_IBS <- ggplot() + 
  geom_boxplot(data = IBS_final_ggplot, aes(x = factor(L1, levels = cancers_used), 
                                          y = value, fill = ind),
               position = position_dodge(0.75)) + 
  theme_Publication() + xlab("Cancer") + ylab("IBS")+
  theme(legend.title=element_blank(), legend.text = element_text(size = 15)) + 
  geom_text(aes(x = 1:length(cancers_used), y = rep(0.5, length(cancers_used)),
                label = stars_IBS[cancers_used]), col = "darkorchid", size = 5) +
  scale_fill_manual(values=c("brown3","#386cb0", "orchid"),
                    labels = c("Clinical", "mRNA", "Mixed")) +
  geom_text(aes(x = 1:length(cancers_used), y = rep(0, length(cancers_used)), 
                label = clin_text_plot[cancers_used]), col = "brown3", size = 4) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1)) +
  geom_bracket(xmin = length(cancers_used) - 0.25, 
               xmax = length(cancers_used) + 0.25, y.position = 0.45, 
               label = "", col = "orchid", size = 1)

plot_clin_mRNA_IBS

plot_clin_mRNA <- ggarrange(plot_clin_mRNA_C, plot_clin_mRNA_IBS, nrow = 2, ncol = 1,
                             labels = c("A", "B"))
plot_clin_mRNA <- annotate_figure(plot_clin_mRNA, top = text_grob("Cancers used"))
print(plot_clin_mRNA)

ggsave(plot_clin_mRNA, filename = "pdf/suppl_clin_mRNA_can_used.pdf")
print("Figure saved in'pdf/suppl_clin_mRNA_can_used.pdf'")

