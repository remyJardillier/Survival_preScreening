
# Learn models ------------------------------------------------------------


if(learn_new_models){
  for(cancer in cancers_all){
    
    print(paste("*** start learning for", cancer, "***"))
    
    # load the data ---
    source(file = "load_data/load_data_final.R")
    
    n_patients <- nrow(clin_data_cox)
    n_patients
    
    # compute the C-index, IBS and p-value of the univariate Cox ---
    C_df_final = IBS_df_final = p_val_df_final <- 
      data.frame(matrix(ncol = 3, nrow = K_folds * n_rep))
    
    colnames(C_df_final) = colnames(IBS_df_final) = colnames(p_val_df_final) <- 
      c("clin", "mRNA", "both")
    
    ind <- 1
   
    for(i in 1:n_rep){
      
      print(paste0("*** Start learning for repetition number: ", i, " ***"))
      
      flds <- createFolds(1:nrow(count_mRNA), k = K_folds, list = TRUE, returnTrain = FALSE)
      
      for(k in 1:K_folds){
        
        print(paste0("Start learning for fold: ", k, " ***"))
        
        # build training and testing dataset
        id_test <- flds[[k]]
        
        # id of the training dataset
        id_train <- 1:nrow(count_mRNA)
        id_train <- id_train[-id_test]
        
        # testing and training dataset
        count_train <- count_mRNA[id_train,]
        count_test <- count_mRNA[id_test,]
        
        clin_train <- clin[id_train,]
        clin_test <- clin[id_test,]
        
        clin_cox_train <- clin_data_cox[id_train,]
        clin_cox_test <- clin_data_cox[id_test,]
        
        dim(count_train)
        dim(count_test)
        
        # filtering ang logCPM in the training data
        logCPM_list <- log.cpm.cv(count_train, count_test)
        
        logCPM_train <- logCPM_list$logCPM_train
        sd_train <- apply(logCPM_train, 2, sd)
        logCPM_train <- logCPM_train[, sd_train != 0]
        
        logCPM_test <- logCPM_list$logCPM_test
        logCPM_test <- logCPM_test[, colnames(logCPM_train)]
        
        # standardization
        logCPM_train_std <- std_train(logCPM_train)
        logCPM_test_std <- std_test(logCPM_train, logCPM_test)
        
        m_age_train <- mean(clin_cox_train$age)
        sd_age_train <- sd(clin_cox_train$age)
        clin_cox_train$age <- (clin_cox_train$age - m_age_train) / sd_age_train  
        clin_cox_test$age <- (clin_cox_test$age - m_age_train) / sd_age_train 
        
        # mRNA only ---
        res_mRNA <- try(fit <- learn_models(logCPM_train_std, clin_train, "ridge"))
        
        # compute the C-index, IBS, and p-value of the univariate Cox
        if(!inherits(res_mRNA, "try-error")){
          
          # prognostic indices
          beta <- as.numeric(coef(fit, "lambda.min"))
          names(beta) <- colnames(logCPM_train_std) 
          PI_test_mRNA <-  logCPM_test_std %*% beta
          PI_train_mRNA <- logCPM_train_std %*% beta
          
          # choose if the data have to be used standardized or not
          C_IBS_pVal_tmp <- C_IBS_pVal_func(fit, clin_train, clin_test, logCPM_train_std, 
                                            logCPM_test_std, "ridge")
          
          # lambda.min
          C_df_final[ind + k - 1, "mRNA"] <- C_IBS_pVal_tmp[1]
          IBS_df_final[ind + k - 1, "mRNA"] <- C_IBS_pVal_tmp[2]
          p_val_df_final[ind + k - 1, "mRNA"] <- C_IBS_pVal_tmp[3]
          
        }
        
        # clinical only ---
        # learn a model 
        res <- try(fit <- learn_models(clin_cox_train, clin_train, "coxph"))
        
        # compute the C-index, IBS, and p-value of the univariate Cox
        if(!inherits(res, "try-error")){
          
          C_IBS_pVal_tmp <- C_IBS_pVal_func(fit, clin_train, clin_test, clin_cox_train, 
                                            clin_cox_test, "coxph")
          
          C_df_final[ind + k - 1, "clin"] <- C_IBS_pVal_tmp[1]
          IBS_df_final[ind + k - 1, "clin"] <- C_IBS_pVal_tmp[2]
          p_val_df_final[ind + k - 1, "clin"] <- C_IBS_pVal_tmp[3]
        }
        
        # clin + mRNA ---
        # learn a model
        if(!inherits(res_mRNA, "try-error")){
          
          mRNA_clin_train <- cbind(PI = PI_train_mRNA, clin_cox_train)
          mRNA_clin_test <- cbind(PI = PI_test_mRNA, clin_cox_test)
          
          res <- try(fit <- learn_models(mRNA_clin_train, clin_train, "coxph"))
          
          # compute the C-index, IBS, and p-value of the univariate Cox
          if(!inherits(res, "try-error")){
            
            C_IBS_pVal_tmp <- C_IBS_pVal_func(fit, clin_train, clin_test, mRNA_clin_train, 
                                              mRNA_clin_test, "coxph")
            
            C_df_final[ind + k - 1, "both"] <- C_IBS_pVal_tmp[1]
            IBS_df_final[ind + k - 1, "both"] <- C_IBS_pVal_tmp[2]
            p_val_df_final[ind + k - 1, "both"] <- C_IBS_pVal_tmp[3]
          }
        }
        
        
      }
      
      ind <- ind + K_folds
    }
    
    save(C_df_final, IBS_df_final, p_val_df_final, n_patients,
         file = paste0("data_fit/", cancer, "/ridge/pred_clin_mRNA.RData"))
    print(paste0("Data saved in: ", "'data_fit/", cancer, "/ridge/pred_clin_mRNA.RData'"))
  }
}






