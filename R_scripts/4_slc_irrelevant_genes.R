

# optimal thresholds
load(file = paste0("data_fit/", cancer, "/", method, "/pred_flt.RData"))
higher <- best_thrs(C_ary_final, "max")

# load the data
source(file = "load_data/load_data_final.R")

# normalize the data
logCPM <- log.cpm(count_mRNA)
logCPM_std <- scale(logCPM)

VST <- t(vst(t(round(count_mRNA))))

# learn a model
fit <- learn_models(logCPM_std, clin, method)
genes_slc <- genes_selected_func(fit, "lambda.min")

# p-values of the genes selected
p_val_all <- p_val_univCox_func(logCPM_std, y_cox)
p_val_all_BH <- p.adjust(p_val_all, method = "BH")
p_val_slc <- p_val_all_BH[genes_slc]
hist(p_val_slc)

plot_p_val <- ggplot() +
  geom_histogram(data = data.frame(p_val = p_val_slc), aes(x = p_val), 
                 col = "black", fill = "lightblue") + 
  theme_Publication() + xlab("Corrected p-values") +
  geom_vline(xintercept = higher[1], linetype = "dashed", 
             color = "red", size = 2)
plot_p_val

sum(p_val_slc > higher[1])
sum(p_val_slc <= higher[1])
length(genes_slc)

# IQR among patients
IQR_VST_slc <- apply(VST[,genes_slc], 2, IQR)

sum(IQR_VST_slc > higher[2])
sum(IQR_VST_slc <= higher[2])
length(genes_slc)

plot_IQR <- ggplot() +
  geom_histogram(data = data.frame(IQR = IQR_VST_slc), aes(x = IQR), 
                 col = "black", fill = "lightblue") + 
  theme_Publication() + xlab("IQR") +
  geom_vline(xintercept = higher[2], linetype = "dashed", 
             color = "red", size = 2)
plot_IQR

# genes selected
genes_slc_plot <- ggarrange(plot_p_val, plot_IQR, ncol = 2, nrow = 1,
                            labels = c("C", "D"))
genes_slc_plot <- annotate_figure(genes_slc_plot, 
                                  top = text_grob(cancer, face = "bold", size = 14))
print(genes_slc_plot)

ggsave(genes_slc_plot, filename = "pdf/suppl_irrelvant_genes.pdf")
print("Figure saved in 'pdf/suppl_irrelvant_genes.pdf' file")

