
# grid --------------------------------------------------------------------


grid_plot_one <- function(ary_final, n_genes_ary_final, max_or_min, name, main = ""){
  
  med_df <- apply(ary_final, c(1,2), function(x) median(x, na.rm = T))
  med_n_genes_df <- round(apply(n_genes_ary_final, c(1,2), function(x) median(x, na.rm = T)))
  
  # data frame for the plot
  med_df <- melt(med_df)
  med_n_genes_df <- melt(med_n_genes_df)
  
  # str(med_C_df)
  med_df[,c(1,2)] <- apply(med_df[,c(1,2)], 2, as.character)
  med_n_genes_df[,c(1,2)] <- apply(med_n_genes_df[,c(1,2)], 2, as.character)
  
  # plot ---
  higher <- best_thrs(ary_final, max_or_min)
  id_p_val_best <- which(thrs_p_val %in% higher[1]) 
  id_IQR_best <- which(thrs_IQR %in% higher[2]) 
  
  # inverse the colors or not
  if(max_or_min == "min"){
    trans <- 'reverse'
  }else{
    trans <- 'identity'
  }
  
  
  grid_plot <- ggplot() + 
    geom_tile(data = med_df, aes(x=Var1, y=Var2, fill=value))  + xlab("Threshold (p-value)") + 
    ylab("Threshold (IQR)") + theme_Publication_legend_right() +
    geom_text(data = med_n_genes_df, aes(x=Var1, y=Var2, label=value), color = "white", size = 5) +
    geom_rect(size=1., fill=NA, colour="grey",
              aes(xmin=length(thrs_p_val) - 0.5, xmax=length(thrs_p_val) + 0.5, ymin=1 - 0.5, ymax=1+ 0.5)) +
    geom_rect(size=1., fill=NA, colour="#386cb0",
              aes(xmin=id_p_val_best - 0.5, xmax=id_p_val_best + 0.5, 
                  ymin=id_IQR_best - 0.5, ymax=id_IQR_best+ 0.5)) +
    scale_fill_gradientn(colours = pal, name = name, trans = trans) +
    ggtitle(main)
  
  return(grid_plot)
}

grid_plot_one_IQR <- function(ary_final, n_genes_ary_final, max_or_min, name, main = ""){
  
  IQR_df <- apply(ary_final, c(1,2), function(x) IQR(x, na.rm = T))
  med_n_genes_df <- round(apply(n_genes_ary_final, c(1,2), function(x) median(x, na.rm = T)))
  
  # data frame for the plot
  IQR_df <- melt(IQR_df)
  med_n_genes_df <- melt(med_n_genes_df)
  
  # str(med_C_df)
  IQR_df[,c(1,2)] <- apply(IQR_df[,c(1,2)], 2, as.character)
  med_n_genes_df[,c(1,2)] <- apply(med_n_genes_df[,c(1,2)], 2, as.character)
  
  # plot ---
  higher <- best_thrs(ary_final, max_or_min)
  id_p_val_best <- which(thrs_p_val %in% higher[1]) 
  id_IQR_best <- which(thrs_IQR %in% higher[2]) 
  
  # inverse the colors or not
  if(max_or_min == "min"){
    trans <- 'reverse'
  }else{
    trans <- 'identity'
  }
  
  
  grid_plot <- ggplot() + 
    geom_tile(data = IQR_df, aes(x=Var1, y=Var2, fill=value))  + xlab("Threshold (p-value)") + 
    ylab("Threshold (IQR)") + theme_Publication_legend_right() +
    geom_text(data = med_n_genes_df, aes(x=Var1, y=Var2, label=value), color = "white", size = 4) +
    geom_rect(size=1., fill=NA, colour="grey",
              aes(xmin=6 - 0.5, xmax=6 + 0.5, ymin=1 - 0.5, ymax=1+ 0.5)) +
    geom_rect(size=1., fill=NA, colour="#386cb0",
              aes(xmin=id_p_val_best - 0.5, xmax=id_p_val_best + 0.5, 
                  ymin=id_IQR_best - 0.5, ymax=id_IQR_best+ 0.5)) +
    scale_fill_gradientn(colours = pal, name = name, trans = trans) +
    ggtitle(main)
  
  return(grid_plot)
}


# raw vs filtering --------------------------------------------------------

raw_vs_flt_plot_nested <- function(cancer, p_val_df_BH_list, ylim = c(0.4,1), 
                                   type = "C", main = ""){
  
  # p-values with BH correction ---
  p_val_BH <- p_val_df_BH_list[[type]][cancer,]
  stars_BH <- sapply(p_val_BH, my_stars.pval)
  
  # load and store the data from nested cross validation ... ---
  df_list <- list()
  
  # ... and number of genes retained by the pre-filtering ...
  n_genes_flt = n_genes_IQR = n_genes_p_val <- rep(NA, length(methods_vect))
  names(n_genes_flt) = names(n_genes_IQR) = names(n_genes_p_val) <- methods_vect
  
  # ... and threshold in the optimal situation
  thrs_opt_IQR = thrs_opt_p_val <- rep(NA, length(methods_vect))
  names(thrs_opt_IQR) = names(thrs_opt_p_val) <- methods_vect
  
  for(method in methods_vect){
    
    load(file = paste0("data_fit/", cancer, "/", method, "/pred_flt.RData"))
    
    if(type == "C"){
      higher <- best_thrs(C_ary_final, "max")
      ylab <- "C-index"
      
      y_n_genes <- ylim[1]+0.1
      y_n_genes2 <- ylim[1]+0.05
      
    }else if(type == "IBS"){
      higher <- best_thrs(IBS_ary_final, "min")
      ylab <- "IBS"
      
      y_n_genes <- ylim[1]+0.1
      y_n_genes2 <- ylim[1]+0.05
    }
    
    id_p_val_best <- which(thrs_p_val %in% higher[1]) 
    id_IQR_best <- which(thrs_IQR %in% higher[2])
    
    thrs_opt_IQR[method] <- higher[2]
    thrs_opt_p_val[method] <- higher[1]
    
    n_genes_flt[method] <- round(median(n_genes_ary_final[id_p_val_best, id_IQR_best,], na.rm = T))
    n_genes_IQR[method] <- round(median(n_genes_IQR_ary_final[id_p_val_best, id_IQR_best,], na.rm = T))
    n_genes_p_val[method] <- round(median(n_genes_p_val_ary_final[id_p_val_best, id_IQR_best,], na.rm = T))
    
    # save the information in a data frame ---
    if(type == "C"){
      
      load(file = paste0("data_fit/", cancer, "/C_nested.RData"))
      
      df_list[[method]] <- data.frame(raw = C_df_nested[,method], 
                           flt = C_df_nested_flt[,method])
      
    }else if(type == "IBS"){
      
      load(file = paste0("data_fit/", cancer, "/IBS_nested.RData"))
      
      df_list[[method]] <- data.frame(raw = IBS_df_nested[,method], 
                           flt = IBS_df_nested_flt[,method])
    }
    
  }
  
  
  
  df_final <- melt(df_list)
  df_final$L1 <- factor(df_final$L1, methods_vect)
  
  ggplot() +
    geom_boxplot(data = df_final, aes(x = L1, y = value, fill = variable), 
                 position = position_dodge(0.75)) + theme_Publication() +
    discrete_scale("fill","Publication",manual_pal(values = c("grey","#386cb0"))) +
    xlab("Method") + ylab(ylab) +
    geom_text(aes(x = 1:length(methods_vect), y = rep(ylim[2], length(methods_vect)), 
                  label = stars_BH), size = 5) +
    geom_text(aes(x = 1:length(methods_vect), y = rep(y_n_genes, length(methods_vect)), 
                  label = n_genes_flt), size = 5, col = "#386cb0") +
    geom_text(aes(x = 1:length(methods_vect), y = rep(y_n_genes2, length(methods_vect)), 
                  label = paste(n_genes_p_val, "/", n_genes_IQR)), size = 5) +
    geom_text(aes(x = 1:length(methods_vect), y = rep(ylim[1], length(methods_vect)), col = "darkred", fontface=2, 
                  label = paste(thrs_opt_p_val, "/", thrs_opt_IQR)), size = 5) +
   ggtitle(main)
    
}







# Venn diagramm -----------------------------------------------------------

venn_digramm <- function(ary_final, n_genes_IQR_ary_final, n_genes_p_val_ary_final,
                         n_genes_ary_final, max_or_min){
  
  # optimal threshold
  higher <- best_thrs(ary_final, max_or_min)
  id_p_val_best <- which(thrs_p_val %in% higher[1]) 
  id_IQR_best <- which(thrs_IQR %in% higher[2]) 
  
  # median number of genes for each filter
  med_n_genes_IQR <- round(median(n_genes_IQR_ary_final[id_p_val_best, id_IQR_best,], na.rm = T))
  med_n_genes_p_val <- round(median(n_genes_p_val_ary_final[id_p_val_best, id_IQR_best,], na.rm = T))
  med_n_genes_tot <- round(median(n_genes_ary_final[id_p_val_best, id_IQR_best,], na.rm = T))
  
  grid.newpage()
  draw.pairwise.venn(med_n_genes_IQR, med_n_genes_p_val, med_n_genes_tot, 
                     category = c("IQR", "p-value"), 
                     col = c("#440154ff", '#21908dff'),
                     fill = c("#440154ff", '#21908dff'), alpha = rep(0.3, 2), 
                     cat.pos = c(0, 0), cat.dist = rep(0.025, 2),
                     fontfamily = rep("Helvetica", 3), cat.fontfamily = rep("Helvetica", 2),
                     cex = rep(1.5, 3), cat.cex = rep(1.5, 2))
}



# SIS ---------------------------------------------------------------------

boxplot_SIS <- function(cancers_vect, cancers_all, type, method, ylim){
  
  
  df_ggplot <- list()
  n_cancer <- length(cancers_vect)
  
  n_genes_flt = n_genes_S <- rep(NA, n_cancer)
  names(n_genes_flt) = names(n_genes_S) <- cancers_vect
  
  for(cancer in cancers_vect){
    
    if(type == "C"){
      load(file = paste0("data_fit/", cancer, "/C_nested.RData"))
      ylab <- "C-index"
      y_text <- ylim[1] + 0.07
      
    }else if(type == "IBS"){
      load(file = paste0("data_fit/", cancer, "/IBS_nested.RData"))
      ylab <- "IBS"
      y_text <- ylim[1] + 0.05
    }
    
    # number of genes selected
    n_genes_flt[cancer] <- round(median(n_genes_slc_flt[, method], na.rm = T))
    
    load(file = paste0("data_fit/", cancer, "/ridge/SIS.RData"))
    n_genes_S[cancer] <- round(median(n_genes_SIS, na.rm = T))
    
    # save the information in a data frame ---
    if(type == "C"){
      df_ggplot[[cancer]] <- data.frame(raw = C_df_nested_flt[, method], 
                                        screen = C_df_nested[, method],
                                        ISIS = C_SIS)
      
    }else if(type == "IBS"){
      df_ggplot[[cancer]] <- data.frame(raw = IBS_df_nested_flt[, method], 
                                        screen = IBS_df_nested[, method],
                                        ISIS = IBS_SIS)
    }
    
  }
  
  # compute the p-values with BH correction
  p_val_wilcox = col_p_val <- rep(NA, length(cancers_all))
  names(p_val_wilcox) = names(col_p_val) <- cancers_all
  
  for(cancer in cancers_all){
    
    load(file = paste0("data_fit/", cancer, "/ridge/SIS.RData"))
    
    # compute the p-value (wilcoxon test)
    if(type == "C"){
      load(file = paste0("data_fit/", cancer, "/C_nested.RData"))
      
      df_tmp <- data.frame(flt = C_df_nested_flt[, method],
                           SIS = C_SIS)
      p_val_wilcox[cancer] <- wilcox.test(df_tmp$flt, df_tmp$SIS)$p.value
      
      if(median(df_tmp$flt, na.rm = T) >= median(df_tmp$SIS, na.rm = T)){
        col_p_val[cancer] <- "#7fc97f"
      }else{
        col_p_val[cancer] <- "#ef3b2c"
      }
      
      
    }else if(type == "IBS"){
      load(file = paste0("data_fit/", cancer, "/IBS_nested.RData"))
      
      df_tmp <- data.frame(flt = IBS_df_nested_flt[, method],
                           SIS = IBS_SIS)
      p_val_wilcox[cancer] <- wilcox.test(df_tmp$flt, df_tmp$SIS)$p.value
      
      if(median(df_tmp$flt, na.rm = T) >= median(df_tmp$SIS, na.rm = T)){
        col_p_val[cancer] <- "#ef3b2c"
      }else{
        col_p_val[cancer] <- "#7fc97f"
      }
    }
    
  }
  p_val_wilcox_BH <-p.adjust(p_val_wilcox, method = "BH")
  p_val_stars <- sapply(p_val_wilcox_BH, my_stars.pval)  
  
  # ggplot
  df_ggplot <- melt(df_ggplot)
  
    comp_SIS <- ggplot() +
      geom_boxplot(data = df_ggplot, aes(x = factor(L1, levels = cancers_vect),
                                         y = value, fill=variable),
                   position = position_dodge(0.75)) + 
      theme_Publication_legend_right() + xlab("Cancer") + ylab(ylab) + ylim(ylim[1],ylim[2]) +
      discrete_scale("fill","Publication",manual_pal(values = c("grey","#7fc97f", "#ef3b2c"))) +
      geom_text(aes(x = 1:n_cancer, y = rep(y_text, n_cancer), 
                    label = n_genes_flt), size = 5, col = "#7fc97f")+
      geom_text(aes(x = 1:n_cancer, y = rep(ylim[1], n_cancer), 
                    label = n_genes_S), size = 5, col = "#ef3b2c") +
      geom_text(aes(x = 1:n_cancer, y = rep(ylim[2], n_cancer), 
                    label = p_val_stars[cancers_vect]), size = 5, col = col_p_val[cancers_vect]) +
      theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
            legend.title=element_blank(), legend.text = element_text(size = 15)) +
      geom_bracket(xmin = length(cancers_all), 
                   xmax = length(cancers_all) + 0.25, y.position = 0.92, 
                   label = "", col = "black", size = 1)
  
  
  return(comp_SIS)
}

