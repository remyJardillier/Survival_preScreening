

# Learn models ------------------------------------------------------------

for(cancer in cancers_all){
  
  print(paste0("Start learning for: ", cancer))
  
  # load the data ---
  source(file = "load_data/load_data_final.R")
  
  # compute the C-index, IBS and p-value of the univariate Cox
  # for all pre-filtering thresholds 
  C_ary_final = IBS_ary_final =
  n_genes_ary_final = n_genes_IQR_ary_final = n_genes_p_val_ary_final <- 
    array(dim=c(length(thrs_p_val), length(thrs_IQR),  K_folds * n_rep),
                dimnames = list(thrs_p_val, thrs_IQR, 1:(K_folds * n_rep)))
  
  # i=1
  # k=1
  
  ind <- 1
  
  for(i in 1:n_rep){
    
    print(paste0("*** Start learning for repetition number: ", i, " (", cancer, ") ***"))
    
    flds <- createFolds(1:nrow(count_mRNA), k = K_folds, list = TRUE, returnTrain = FALSE)
    
    for(k in 1:K_folds){
      
      print(paste0("*** Start learning for fold: ", k, " ***"))
      
      # build training and testing dataset ---
      id_test <- flds[[k]]
      
      # id of the training dataset
      id_train <- 1:nrow(count_mRNA)
      id_train <- id_train[-id_test]
      
      # testing and training dataset
      count_train <- count_mRNA[id_train,]
      count_test <- count_mRNA[id_test,]
      
      clin_train <- clin[id_train,]
      y_cox_train <- Surv(clin_train$time, clin_train$status)
      clin_test <- clin[id_test,]
      y_cox_test <- Surv(clin_test$time, clin_test$status)
      
      dim(count_train)
      dim(count_test)
      
      # filtering (genes expressed or not) ang logCPM in the training data
      logCPM_list <- log.cpm.cv(count_train, count_test)
      
      logCPM_train <- logCPM_list$logCPM_train
      sd_train <- apply(logCPM_train, 2, sd)
      logCPM_train <- logCPM_train[, sd_train != 0]
      count_train <- count_train[, colnames(logCPM_train)]
      
      logCPM_test <- logCPM_list$logCPM_test
      logCPM_test <- logCPM_test[, colnames(logCPM_train)]
      
      # standardization
      logCPM_train_std <- std_train(logCPM_train)
      logCPM_test_std <- std_test(logCPM_train, logCPM_test)
      
      # p-values and IQR for the filtering ---
      # compute p-values
      p_val_cox_train <- p_val_univCox_func(logCPM_train_std, y_cox_train)
      p_val_cox_train_BH <- p.adjust(p_val_cox_train, method = "BH")
      
      id_NA <- is.na(p_val_cox_train_BH)
      logCPM_train_std <- logCPM_train_std[,!id_NA]
      count_train <- count_train[,!id_NA]
      logCPM_test_std <- logCPM_test_std[, !id_NA]
      p_val_cox_train_BH <-  p_val_cox_train_BH[!id_NA]
      
      dim(logCPM_train_std)
      dim(count_train)
      dim(logCPM_test_std)
      length(p_val_cox_train_BH)
      
      # compute IQR with variance stabilization
      vst_train <- t(vst(t(round(count_train))))
      vst_train[1:5, 1:5]
      
      IQR_vect_train <- apply(vst_train, 2, IQR)
      
      
      for(p in 1:length(thrs_p_val)){
        
        print(paste0("Start learning for p_val_univCox: ", thrs_p_val[p]))
        
        for(q in 1:length(thrs_IQR)){
          
          print(paste0("Start learning for thrs_IQR: ", thrs_IQR[q]))
          
          if(sum(p_val_cox_train_BH <= thrs_p_val[p] & 
                 IQR_vect_train >= thrs_IQR[q]) >= 2){
            
            # pre-filter the data and learn the models
            logCPM_train_std_tmp <- logCPM_train_std[, p_val_cox_train_BH <= thrs_p_val[p] & 
                                                         IQR_vect_train >= thrs_IQR[q]]
            logCPM_test_std_tmp <- logCPM_test_std[, p_val_cox_train_BH <= thrs_p_val[p] & 
                                                       IQR_vect_train >= thrs_IQR[q]]
            
            # learn a model
            res <- try(fit <- learn_models(logCPM_train_std_tmp, clin_train, method))
            
            # compute the C-index, IBS and p-value
            if(!inherits(res, "try-error")){
              
              # number of genes
              n_genes_ary_final[p,q,ind + k-1] <- sum(p_val_cox_train_BH <= thrs_p_val[p] & 
                                                      IQR_vect_train >= thrs_IQR[q])
              n_genes_IQR_ary_final[p,q,ind + k-1] <- sum(IQR_vect_train >= thrs_IQR[q])
              n_genes_p_val_ary_final[p,q,ind + k-1] <- sum(p_val_cox_train_BH <= thrs_p_val[p])
              
              # C-index, IBS and p-value
              C_IBS_pVal_tmp <- C_IBS_pVal_func(fit, clin_train, clin_test, logCPM_train_std_tmp, 
                                                logCPM_test_std_tmp, method)
             
              
              C_ary_final[p,q,ind + k-1] <- C_IBS_pVal_tmp[1]
              IBS_ary_final[p,q,ind + k-1] <- C_IBS_pVal_tmp[2]
            }
          }
        }
      }
    }
      
    ind <- ind + K_folds
  }
  
  save(C_ary_final, IBS_ary_final, 
       n_genes_ary_final, n_genes_IQR_ary_final, n_genes_p_val_ary_final,
       file = paste0("data_fit/", cancer, "/", method, "/pred_flt.RData"))
}



# grid + venn diagramm - all cancers --------------------------------------

method <- "EN"

if (!dir.exists("pdf_tmp"))
  dir.create("pdf_tmp", recursive = T)

for(cancer in cancers_all){
  
  print(paste0("Start producing the figure for: ", cancer))
  
  load(file = paste0("data_fit/", cancer, "/", method, "/pred_flt.RData"))
  
  # C-index
  grid_C <- grid_plot_one(C_ary_final, n_genes_ary_final, max_or_min = "max", name = "C-index")
  venn_C <- venn_digramm(C_ary_final, n_genes_IQR_ary_final, n_genes_p_val_ary_final,
                         n_genes_ary_final, "max")
  
  # IBS
  grid_IBS <- grid_plot_one(IBS_ary_final, n_genes_ary_final, max_or_min = "min", name = "IBS")
  venn_IBS <- venn_digramm(IBS_ary_final, n_genes_IQR_ary_final, n_genes_p_val_ary_final,
                           n_genes_ary_final, "min")
  
  
  final_plot <- ggarrange(grid_C, venn_C, grid_IBS, venn_IBS,
                          ncol = 2, nrow = 2, widths = rep(c(2,1),3))
  final_plot <- annotate_figure(final_plot,
                                top = text_grob(paste(cancer, "-", method), face = "bold", size = 14))
  # final_plot
  
  ggsave(final_plot, filename = paste0("pdf_tmp/", cancer, "_", method, ".pdf"),
         width = 210, height = 297, units = "mm")
}

pdf_names <- list.files("pdf_tmp", full.names = T)
pdf_combine(pdf_names, output = paste0("pdf/C_grid_", method, "_venn.pdf"))
unlink("pdf_tmp", recursive = T)

print(paste0("Figure saved in 'pdf/C_grid_", method, "_venn.pdf' file"))


# grid + venn diagramm - one cancer ---------------------------------------

cancer <- "ACC"
method <- "EN"

load(file = paste0("data_fit/", cancer, "/", method, "/pred_flt.RData"))

# C-index
grid_C <- grid_plot_one(C_ary_final, n_genes_ary_final, max_or_min = "max", name = "C-index")
venn_C <- venn_digramm(C_ary_final, n_genes_IQR_ary_final, n_genes_p_val_ary_final,
                       n_genes_ary_final, "max")

grid_venn_C <- ggarrange(grid_C, venn_C, ncol = 2, nrow = 1, 
                        widths = c(3,1), labels = c("A", "B"))

# IBS
grid_IBS <- grid_plot_one(IBS_ary_final, n_genes_ary_final, max_or_min = "min", name = "IBS")
venn_IBS <- venn_digramm(IBS_ary_final, n_genes_IQR_ary_final, n_genes_p_val_ary_final,
                       n_genes_ary_final, "min")

grid_venn_IBS <- ggarrange(grid_IBS, venn_IBS, ncol = 2, nrow = 1, 
                         widths = c(3,1), labels = c("A", "B"))
print(grid_venn_IBS)
print(grid_venn_C)
