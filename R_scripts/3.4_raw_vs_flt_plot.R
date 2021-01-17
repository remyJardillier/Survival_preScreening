

# p-values with BH correction ---------------------------------------------

p_val_df_BH_list <- list()

med_inc_list <- list()
med_inc_all_list <- list()

med_n_genes_list <- list()
med_n_genes_all_list <- list()

num_NA_list <- list()

for(type in c("C", "IBS")){ 
  
  # median increase
  med_inc_df <- matrix(ncol = length(methods_vect), 
                       nrow = length(cancers_all))
  colnames(med_inc_df) <- methods_vect
  row.names(med_inc_df) <- cancers_all
  
  # median increase all cancers
  med_inc_all_tmp <- list()
  
  # p-value with BH correction
  p_val_df <- matrix(ncol = length(methods_vect), nrow = length(cancers_all))
  colnames(p_val_df) <- methods_vect
  row.names(p_val_df) <- cancers_all
  
  # median number of genes
  med_n_genes <- matrix(ncol = length(methods_vect), nrow = length(cancers_all))
  colnames(med_n_genes) <- methods_vect
  row.names(med_n_genes) <- cancers_all
  
  # median number of genes - all cancers
  med_n_genes_all_tmp <- list()
  
  # number of NA
  num_NA <- matrix(ncol = length(methods_vect), nrow = length(cancers_all))
  colnames(num_NA) <- methods_vect
  row.names(num_NA) <- cancers_all
  
  # loop over cancers and methods
  for(cancer in cancers_all){
    
    med_inc_all = med_n_genes_all <- matrix(ncol = length(methods_vect), nrow = K_folds * n_rep)
    colnames(med_inc_all) = colnames(med_n_genes_all) <- methods_vect
      
    for(method in methods_vect){
      
      load(file = paste0("data_fit/", cancer, "/", method, "/pred_flt.RData"))
      
      # threshold for the highest median C-index / p-value / IBS
      if(type == "C"){
        higher <- best_thrs(C_ary_final, "max")
      }else if(type == "p"){
        higher <- best_thrs(p_val_ary_final, "min")
      }else if(type == "IBS"){
        higher <- best_thrs(IBS_ary_final, "min")
      }
      
      id_p_val_best <- which(thrs_p_val %in% higher[1]) 
      id_IQR_best <- which(thrs_IQR %in% higher[2])
      
      med_n_genes[cancer,method] <- median(n_genes_ary_final[id_p_val_best, id_IQR_best,], 
                                 na.rm = T)
      med_n_genes_all[, method] <- n_genes_ary_final[id_p_val_best, id_IQR_best,]
      
      # wilcoxon test between the 2 (flt vs no flt)
      if(type == "C"){
        load(file = paste0("data_fit/", cancer, "/C_nested.RData"))
        df_tmp <- data.frame(raw = C_df_nested[,method], 
                             flt = C_df_nested_flt[,method])
        try(p_val_df[cancer,method] <- wilcox.test(df_tmp$flt, df_tmp$raw, 
                                         alternative = "greater",paired = T)$p.value)
        
        
        med_inc_df[cancer,method] <- median(df_tmp$flt - df_tmp$raw, na.rm = T)
        med_inc_all[, method] <- df_tmp$flt - df_tmp$raw
        
        num_NA[cancer,method] <- sum(is.na(df_tmp$flt - df_tmp$raw))
        
      }else if(type == "IBS"){
        
        load(file = paste0("data_fit/", cancer, "/IBS_nested.RData"))
        
        df_tmp <- data.frame(raw = IBS_df_nested[,method], 
                             flt = IBS_df_nested_flt[,method])
        
        try(p_val_df[cancer,method] <- wilcox.test(df_tmp$flt, df_tmp$raw, 
                                         alternative = "less",paired = T)$p.value)
        
        med_inc_df[cancer,method] <- median(df_tmp$raw - df_tmp$flt, na.rm = T)
        med_inc_all[, method] <- df_tmp$raw - df_tmp$flt
        
        num_NA[cancer,method] <- sum(is.na(df_tmp$flt - df_tmp$raw))
      }
      
      
    } 
    
    med_inc_all_tmp[[cancer]] <- med_inc_all
    med_n_genes_all_tmp[[cancer]] <- med_n_genes_all
  }
  
  p_val_df_BH_list[[type]] <- apply(p_val_df, 2, function(x) p.adjust(x, method = "BH"))
  
  med_inc_list[[type]] <- med_inc_df
  med_inc_all_tmp <- do.call(rbind,med_inc_all_tmp)
  med_inc_all_list[[type]] <- apply(med_inc_all_tmp, 2, function(x) median(x, na.rm = T))
  
  med_n_genes_list[[type]] <- med_n_genes
  med_n_genes_all_tmp <- do.call(rbind,med_n_genes_all_tmp)
  med_n_genes_all_list[[type]] <- apply(med_n_genes_all_tmp, 2, function(x) median(x, na.rm = T))
  
  num_NA_list[[type]] <- num_NA
}

p_val_df_BH_list
med_inc_list
med_n_genes_list 
num_NA_list 
med_inc_all_list
med_n_genes_all_list




# plot - all -----------------------------------------------------------------

if (!dir.exists("pdf_tmp"))
  dir.create("pdf_tmp", recursive = T)

for(cancer in cancers_all){
  
  print(paste0("Start producing the figure for: ", cancer))
  
  plot_C <- raw_vs_flt_plot_nested(cancer, p_val_df_BH_list, ylim = c(0.4,1), type = "C")
  #plot_p <- raw_vs_flt_plot_nested(cancer, p_val_df_BH_list, ylim = c(-0.5, 7), type = "p")
  plot_IBS <- raw_vs_flt_plot_nested(cancer, p_val_df_BH_list, ylim = c(-0.1,0.4), type = "IBS")
  
  fig_final <- ggarrange(plot_C, plot_IBS,
                         labels = c("A", "B"), ncol = 1, nrow = 2)
  fig_final <- annotate_figure(fig_final,
                               top = text_grob(cancer, face = "bold", size = 18))
  ggsave(fig_final, filename = paste0("pdf_tmp/", cancer, ".pdf"), 
         width = 210, height = 297, units = "mm")
}

pdf_names <- list.files("pdf_tmp", full.names = T)
pdf_combine(pdf_names, output = "pdf/raw_vs_flt_nested_CV.pdf")
unlink("pdf_tmp", recursive = T)

print("Figure saved in 'pdf/raw_vs_flt_nested_CV.pdf' file")

# one -----------------------------------------------------------------

cancer <- "ACC"
type <- "C"
ylim <- c(0.4,0.9)
main <- ""

flt_vs_raw_nested <- raw_vs_flt_plot_nested(cancer, p_val_df_BH_list, 
                                            ylim = ylim, type = type)

print(flt_vs_raw_nested)


# final plot - main figure ------------------------------------------------

# grid used to compute the optimal thresholds
method <- "EN"

load(file = paste0("data_fit/", cancer, "/", method, "/pred_flt.RData"))

grid_C <- grid_plot_one(C_ary_final, n_genes_ary_final, max_or_min = "max", name = "C-index")
venn_C <- venn_digramm(C_ary_final, n_genes_IQR_ary_final, n_genes_p_val_ary_final,
                       n_genes_ary_final, "max")

grid_venn_C <- ggarrange(grid_C, venn_C, ncol = 2, nrow = 1, 
                         widths = c(3,1), labels = c("A", "B"))

# Prediction obtained with and without filtering after nested cros-validation
ylim <- c(0.4,0.9)
main <- ""
flt_vs_raw_nested <- raw_vs_flt_plot_nested(cancer, p_val_df_BH_list, 
                                            ylim = ylim, type = type)
flt_vs_raw_nested <- ggarrange(flt_vs_raw_nested, labels = "C")

# final plot
final_plot <- ggarrange(grid_venn_C, flt_vs_raw_nested, ncol = 1, nrow = 2)

print(final_plot)

ggsave(final_plot, filename = "Figures/main_fig.pdf")
print("Figure saved in 'Figures/main_fig.pdf' file")

# Tables - sum up of the median increase for all cancers ------------------

# C-index ---
p_val_df_BH_C <- p_val_df_BH_list$C
stars_BH_C <- apply(p_val_df_BH_C, 2, my_stars.pval_vect)
med_inc_df_C <- med_inc_list$C
med_n_genes_C <- med_n_genes_list$C

final_df <- data.frame(matrix(ncol = length(methods_vect)+1, 
                              nrow = 3*length(cancers_all)))
colnames(final_df) <- c("cancer", methods_vect)
final_df$cancer <- rep(cancers_all, each = 3)

for(i in 1:length(cancers_all)){
  final_df[3*i-2,2:(length(methods_vect)+1)] <- signif(med_inc_df_C[i, ], 2)
  final_df[3*i-1,2:(length(methods_vect)+1)] <- stars_BH_C[i, ]
  final_df[3*i,2:(length(methods_vect)+1)] <- signif(med_n_genes_C[i,], 2)
}

print(final_df)
write.xlsx(final_df, file = "tables/sum_up_inc_C_df.xlsx")
print("Table saved in 'tables/sum_up_inc_C_df.xlsx' file")

# IBS ---
p_val_df_BH_IBS <- p_val_df_BH_list$IBS
stars_BH_IBS <- apply(p_val_df_BH_IBS, 2, my_stars.pval_vect)
med_inc_df_IBS <- med_inc_list$IBS
med_n_genes_IBS <- med_n_genes_list$IBS


final_df <- data.frame(matrix(ncol = (length(methods_vect)+1), 
                              nrow = 3*length(cancers_all)))
colnames(final_df) <- c("cancer", methods_vect)
final_df$cancer <- rep(cancers_all, each = 3)

for(i in 1:length(cancers_all)){
  final_df[3*i-2,2:(length(methods_vect)+1)] <- signif(med_inc_df_IBS[i, ], 2)
  final_df[3*i-1,2:(length(methods_vect)+1)] <- stars_BH_IBS[i, ]
  final_df[3*i,2:(length(methods_vect)+1)] <- signif(med_n_genes_IBS,2)
}

print(final_df)
write.xlsx(final_df, file = "tables/sum_up_inc_IBS_df.xlsx")
print("Table saved in 'tables/sum_up_inc_IBS_df.xlsx' file")

# significant increase of the metrics ---
inc_C <- apply(p_val_df_BH_C, 2, function(x) x<0.05)
inc_IBS <- apply(p_val_df_BH_IBS, 2, function(x) x<0.05)

apply(inc_IBS, 2, function(x) sum(x))
apply(inc_C, 2, function(x) sum(x))
apply(inc_C & inc_IBS, 2, function(x) sum(x))


