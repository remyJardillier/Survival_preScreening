
# main - C-index ----------------------------------------------------------

C_list_can <- list()

for(cancer in cancers_all){
  
  C_meth_list <- list()
  
  # load the data
  load(file = paste0("data_fit/", cancer, "/C_nested.RData"))
  
  for(method in methods_vect){
    
    C_meth_list[[method]] <- C_df_nested_flt[, method]
    
  }
  
  C_list_can[[cancer]] <- stack(C_meth_list)
}

C_df_bench <- melt(C_list_can)


pred_C_flt <- ggplot() + 
  geom_boxplot(data = C_df_bench, aes(x= factor(L1, levels = cancers_all),
                                      y = value, fill = ind),
               position = position_dodge(0.75)) +
  theme_Publication_legend_right() + xlab("Cancer") + ylab("C-index") +
  scale_fill_Publication() + theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
                                   legend.text = element_text(size = 13),
                                   legend.title = element_blank())

print(pred_C_flt)

ggsave(pred_C_flt, filename = "Figures/main_pred_C_flt.pdf")
print("Figures saved in 'Figures/main_pred_C_flt.pdf' file.")

# Suppl - IBS ----------------------------------------------------------

IBS_list_can <- list()

for(cancer in cancers_all){
  
  IBS_meth_list <- list()
  
  # load the data
  load(file = paste0("data_fit/", cancer, "/IBS_nested.RData"))
  
  for(method in methods_vect){
    
    IBS_meth_list[[method]] <- IBS_df_nested_flt[, method]
    
  }
  
  IBS_list_can[[cancer]] <- stack(IBS_meth_list)
}

IBS_df_bench <- melt(IBS_list_can)


pred_IBS_suppl_flt <- ggplot() + 
  geom_boxplot(data = IBS_df_bench, aes(x= factor(L1, levels = cancers_all),
                                        y = value, fill = ind),
               position = position_dodge(0.75)) + ylim(NA, 0.5)+ 
  theme_Publication_legend_right() + xlab("Cancer") + ylab("IBS") +
  scale_fill_Publication() + theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
                                   legend.text = element_text(size = 13),
                                   legend.title = element_blank())

print(pred_IBS_suppl_flt)

ggsave(pred_IBS_suppl_flt, filename = "pdf/suppl_pred_IBS_flt.pdf")
print("Figures saved in 'pdf/suppl_pred_IBS_flt.pdf' file.")
