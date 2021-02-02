

# C-index -----------------------------------------------------------------

if(learn_new_models){
  if(pred_measure == "C"){
    
    for(cancer in cancers_all){
      
      print(paste0("*** Start learning for cancer: ", cancer, " ***"))
      
      # load the data
      source(file = "load_data/load_data_final.R")
      
      # load and store the data
      n_genes_slc = C_df_nested = n_genes_slc_flt = C_df_nested_flt <- 
        data.frame(matrix(ncol = length(methods_vect), nrow = n_rep * K_folds))
      colnames(n_genes_slc) = colnames(C_df_nested) = colnames(n_genes_slc_flt) = colnames(C_df_nested_flt) <- 
        methods_vect
      
      for(method in methods_vect){
        
        print(paste0("*** Start learning for method: ", method, " ***"))
        
        load(file = paste0("data_fit/", cancer, "/", method, "/pred_flt.RData"))
        
        higher <- best_thrs(C_ary_final, "max")
        
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
            IQR_vect_train <- apply(vst_train, 2, IQR)
            
            # pre-filter the data and learn the models
            logCPM_train_std_tmp <- logCPM_train_std[, p_val_cox_train_BH <= higher[1] & 
                                                       IQR_vect_train >= higher[2]]
            logCPM_test_std_tmp <- logCPM_test_std[, p_val_cox_train_BH <= higher[1] & 
                                                     IQR_vect_train >= higher[2]]
            
            # learn a model - no flt
            res <- try(fit <- learn_models(logCPM_train_std, clin_train, method))
            
            # compute the C-index, IBS and p-value
            if(!inherits(res, "try-error")){
              
              # number of genes
              n_genes_slc[ind + k-1, method] <- length(genes_selected_func(fit, "lambda.min")) 
              
              # C-index, IBS and p-value
              C_IBS_pVal_tmp <- C_IBS_pVal_func(fit, clin_train, clin_test, logCPM_train_std, 
                                                logCPM_test_std, method)
              C_df_nested[ind + k-1, method] <- C_IBS_pVal_tmp[1]
            }
            
            # learn a model - flt
            res_flt <- try(fit_flt <- learn_models(logCPM_train_std_tmp, clin_train, method))
            
            # compute the C-index, IBS and p-value
            if(!inherits(res_flt, "try-error")){
              
              # number of genes
              n_genes_slc_flt[ind + k-1, method] <- length(genes_selected_func(fit_flt, "lambda.min")) 
              
              # C-index, IBS and p-value
              C_IBS_pVal_tmp <- C_IBS_pVal_func(fit_flt, clin_train, clin_test, logCPM_train_std_tmp, 
                                                logCPM_test_std_tmp, method)
              C_df_nested_flt[ind + k-1, method] <- C_IBS_pVal_tmp[1]
            }
            
          }
          
          ind <- ind + K_folds
        }
        
      }
      save(C_df_nested, n_genes_slc, C_df_nested_flt, n_genes_slc_flt, 
           file = paste0("data_fit/", cancer, "/C_nested.RData"))
      print(paste0("Data saved in: ", "'data_fit/", cancer, "/C_nested.RData'"))
    }
    
  }
}




# IBS --------------------------------------------------------------------

if(learn_new_models){
  if(pred_measure == "IBS"){
    
    for(cancer in cancers_all){
      
      print(paste0("*** Start learning for cancer: ", cancer, " ***"))
      
      source(file = "load_data/load_data_final.R")
      
      # Surv object for the Cox model
      y_cox <- Surv(time = clin$time, event = clin$status)
      
      # load and store the data
      n_genes_slc = IBS_df_nested = n_genes_slc_flt = IBS_df_nested_flt <- 
        data.frame(matrix(ncol = length(methods_vect), nrow = n_rep * K_folds))
      colnames(n_genes_slc) = colnames(IBS_df_nested) = colnames(n_genes_slc_flt) = colnames(IBS_df_nested_flt) <- 
        methods_vect
      
      for(method in methods_vect){
        
        load(file = paste0("data_fit/", cancer, "/", method, "/pred_flt.RData"))
        
        higher <- best_thrs(IBS_ary_final, "min")
        
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
            IQR_vect_train <- apply(vst_train, 2, IQR)
            
            # pre-filter the data and learn the models
            logCPM_train_std_tmp <- logCPM_train_std[, p_val_cox_train_BH <= higher[1] & 
                                                       IQR_vect_train >= higher[2]]
            logCPM_test_std_tmp <- logCPM_test_std[, p_val_cox_train_BH <= higher[1] & 
                                                     IQR_vect_train >= higher[2]]
            
            # learn a model - no flt
            res <- try(fit <- learn_models(logCPM_train_std, clin_train, method))
            
            # compute the C-index, IBS and p-value
            if(!inherits(res, "try-error")){
              
              # number of genes
              n_genes_slc[ind + k-1, method] <- length(genes_selected_func(fit, "lambda.min")) 
              
              # C-index, IBS and p-value
              C_IBS_pVal_tmp <- C_IBS_pVal_func(fit, clin_train, clin_test, logCPM_train_std, 
                                                logCPM_test_std, method)
              IBS_df_nested[ind + k-1, method] <- C_IBS_pVal_tmp[2]
            }
            
            # learn a model - flt
            res_flt <- try(fit_flt <- learn_models(logCPM_train_std_tmp, clin_train, method))
            
            # compute the C-index, IBS and p-value
            if(!inherits(res_flt, "try-error")){
              
              # number of genes
              n_genes_slc_flt[ind + k-1, method] <- length(genes_selected_func(fit_flt, "lambda.min")) 
              
              # C-index, IBS and p-value
              C_IBS_pVal_tmp <- C_IBS_pVal_func(fit_flt, clin_train, clin_test, logCPM_train_std_tmp, 
                                                logCPM_test_std_tmp, method)
              IBS_df_nested_flt[ind + k-1, method] <- C_IBS_pVal_tmp[2]
            }
            
          }
          
          ind <- ind + K_folds
        }
        
      }
      
      save(IBS_df_nested, n_genes_slc, IBS_df_nested_flt, n_genes_slc_flt, 
           file = paste0("data_fit/", cancer, "/IBS_nested.RData"))
      print(paste0("Data saved in: ", "'data_fit/", cancer, "/IBS_nested.RData'"))
      
    }
  }
}




