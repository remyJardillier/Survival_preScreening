# |beta| = f(-log10(p_val)) ---------------------------------------------------


# learn models ---
for(cancer in cancers_all){
  
  print(paste0("Start learning for: ", cancer))
  
  # load the data ---
  source(file = "load_data/load_data_final.R")
  
  # Surv object for the Cox model
  y_cox <- Surv(time = clin$time, event = clin$status)
  
  # logCPM and normalization
  logCPM <- log.cpm(count_mRNA)
  logCPM_std <- scale(logCPM)
  
  # p-values and beta coefficients
  p_val_vect <- rep(NA, ncol(logCPM_std))
  beta_vect <- rep(NA, ncol(logCPM_std))
  
  for(j in 1:ncol(logCPM_std)){
    
    res <- try(coxph_tmp <- coxph(y_cox ~ logCPM_std[, j]))
    
    if(!inherits(res, "try-error")){
      
      test <- summary(coxph_tmp)
      p_val_vect[j] <- test$logtest[3]
      beta_vect[j] <- test$coefficients[1]
      
    }else{
      p_val_vect[j] <- NA
      beta_vect[j] <- NA
    }
    
  }
  p_val_vect_BH <- p.adjust(p_val_vect, method = "BH")
  
  save(p_val_vect_BH, beta_vect, file = paste0("data_fit/", cancer, "/beta_p_val.RData"))
}

# pdf for all cancers + table of correlations
pdf(file = "pdf/beta_p_val.pdf", paper="a4")

cor_beta_p_val_df <- data.frame(matrix(ncol = length(cancers_all), nrow = 2))
colnames(cor_beta_p_val_df) <- cancers_all
row.names(cor_beta_p_val_df) <- c("cor", "stars")

for(cancer in cancers_all){
  
  print(paste0("Start learning for: ", cancer))
  
  # load the data ---
  load(file = paste0("data_fit/", cancer, "/beta_p_val.RData"))
  
  # correlations
  cor_tmp <- cor.test(-log10(p_val_vect_BH), abs(beta_vect), method = "spearman")
  cor_tmp$p.value
  cor_tmp$estimate
  
  cor_beta_p_val_df["cor", cancer] <- signif(cor_tmp$estimate, 3)
  cor_beta_p_val_df["stars", cancer] <- my_stars.pval(cor_tmp$p.value)
  
  # fit
  fit <- loess(abs(beta_vect)~-log10(p_val_vect_BH))
  
  # plot
  plot(-log10(p_val_vect_BH), abs(beta_vect), 
       main = cancer, xlab = "-log10(p)", ylab = "|beta|")
  mtext(paste(signif(cor_tmp$estimate, 3), "/", my_stars.pval(cor_tmp$p.value)))
  j <- order(-log10(p_val_vect_BH))
  lines(-log10(p_val_vect_BH)[j], fit$fitted[j], col="red", lwd = 3)
  grid()
}
par(mfrow = c(3, 1))
dev.off()
source(file = "functions/graphics.R")

print("Beta - p-value correlation plot saved in 'pdf/beta_p_val.pdf'")

write.xlsx(cor_beta_p_val_df, file = "tables/cor_beta_p_val_df.xlsx", row.names = T)


# abs(beta) = f(-log10(p_val)) - one ---

cancer <- "ACC"

load(file = paste0("data_fit/", cancer, "/beta_p_val.RData"))

# correlations
cor_tmp <- cor.test(-log10(p_val_vect_BH), abs(beta_vect), method = "spearman")
cor_tmp$p.value
cor_tmp$estimate

# fit
fit <- loess(abs(beta_vect)~-log10(p_val_vect_BH))

df_cor <- data.frame(p = -log10(p_val_vect_BH), beta = abs(beta_vect))

cor_beta_p <- ggplot(df_cor, aes(x = p, y = beta)) +
  geom_point() +theme_Publication() +
  xlab("-log10(p)") + ylab("|beta|") +
  ggtitle(cancer) +
  geom_smooth(method=loess,  col = "red")
print(cor_beta_p)

ggsave(cor_beta_p, filename = paste0("pdf/beta_p_val_cor_", cancer, ".pdf"))



# Learn models ------------------------------------------------------------

for(cancer in cancers_all){
  
  print(paste0("*** Start learning for cancer: ", cancer, " ***"))
  
  source(file = "load_data/load_data_final.R")
  
  # Surv object for the Cox model
  y_cox <- Surv(time = clin$time, event = clin$status)
  
  # compute the C-index, IBS and p-value of the univariate Cox
  # for all pre-filtering thresholds 
  C_SIS = IBS_SIS = n_genes_SIS <- rep(NA, K_folds * n_rep)
  
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
      
      # Sure Independance Screening (SIS) ---
      # number of genes to keep in the model
      nsis <- floor(nrow(logCPM_train_std) / (log(nrow(logCPM_train_std))))
      
      # Sure indepedance screening
      res_SIS <- try(fit_SIS <- SIS(x = logCPM_train_std, y = y_cox_train, family = "cox",
                                    penalty = "lasso", tune = "cv", standardize = F, nsis = nsis))
      
      if(!inherits(res_SIS, "try-error") & length(fit_SIS$sis.ix0) > 1){
        
        # ridge regression with the remaining genes
        logCPM_train_std_sis <- logCPM_train_std[, fit_SIS$sis.ix0]
        logCPM_test_std_sis <- logCPM_test_std[, fit_SIS$sis.ix0]
        fit_SIS_ridge <- learn_models(logCPM_train_std_sis, clin_train, method)
        
        n_genes_SIS[ind + k - 1] <- length(fit_SIS$sis.ix0)
        
        # compute prognostic indicator
        C_IBS_pVal_SIS <- C_IBS_pVal_func(fit_SIS_ridge, clin_train, clin_test, logCPM_train_std_sis,
                                          logCPM_test_std_sis, method)
        C_SIS[ind + k - 1] <- C_IBS_pVal_SIS[1]
        IBS_SIS[ind + k - 1] <- C_IBS_pVal_SIS[2]
        
      }else{
        
        # number of genes to keep in the model
        nsis <- floor(nrow(logCPM_train_std) / (2*log(nrow(logCPM_train_std))))
        
        # Sure indepedance screening
        res_SIS <- try(fit_SIS <- SIS(x = logCPM_train_std, y = y_cox_train, family = "cox",
                                      penalty = "lasso", tune = "cv", standardize = F, nsis = nsis))
        
        if(!inherits(res_SIS, "try-error") & length(fit_SIS$sis.ix0) > 1){
          
          # ridge regression with the remaining genes
          logCPM_train_std_sis <- logCPM_train_std[, fit_SIS$sis.ix0]
          logCPM_test_std_sis <- logCPM_test_std[, fit_SIS$sis.ix0]
          fit_SIS_ridge <- learn_models(logCPM_train_std_sis, clin_train, method)
          
          n_genes_SIS[ind + k - 1] <- length(fit_SIS$sis.ix0)
          
          # compute prognostic indicator
          C_IBS_pVal_SIS <- C_IBS_pVal_func(fit_SIS_ridge, clin_train, clin_test, logCPM_train_std_sis,
                                            logCPM_test_std_sis, method)
          C_SIS[ind + k - 1] <- C_IBS_pVal_SIS[1]
          IBS_SIS[ind + k - 1] <- C_IBS_pVal_SIS[2]
          
          
        }else{
          
          # number of genes to keep in the model
          nsis <- floor(nrow(logCPM_train_std) / (4*log(nrow(logCPM_train_std))))
          
          # Sure indepedance screening
          res_SIS <- try(fit_SIS <- SIS(x = logCPM_train_std, y = y_cox_train, family = "cox",
                                        penalty = "lasso", tune = "cv", standardize = F, nsis = nsis))
          
          if(!inherits(res_SIS, "try-error") & length(fit_SIS$sis.ix0) > 1){
            
            # ridge regression with the remaining genes
            logCPM_train_std_sis <- logCPM_train_std[, fit_SIS$sis.ix0]
            logCPM_test_std_sis <- logCPM_test_std[, fit_SIS$sis.ix0]
            fit_SIS_ridge <- learn_models(logCPM_train_std_sis, clin_train, method)
            
            n_genes_SIS[ind + k - 1] <- length(fit_SIS$sis.ix0)
            
            # compute prognostic indicator
            C_IBS_pVal_SIS <- C_IBS_pVal_func(fit_SIS_ridge, clin_train, clin_test, logCPM_train_std_sis,
                                              logCPM_test_std_sis, method)
            C_SIS[ind + k - 1] <- C_IBS_pVal_SIS[1]
            IBS_SIS[ind + k - 1] <- C_IBS_pVal_SIS[2]
            
          }
        }
      }
    }
    
    ind <- ind + K_folds
  }
  
  save(C_SIS, IBS_SIS, n_genes_SIS,
       file = paste0("data_fit/", cancer, "/ridge/SIS.RData"))
}


# boxplot -----------------------------------------------------------------

method <- "EN"

# main - C
cancers_vect <- cancers_all

comp_SIS_C_main <- boxplot_SIS(cancers_vect, cancers_all, "C", method, c(0.35, 1))
print(comp_SIS_C_main)

ggsave(comp_SIS_C_main, filename = "Figures/main_SIS_C.pdf")

# suppl. - IBS
comp_SIS_IBS_suppl <- boxplot_SIS(cancers_vect, cancers_all, "IBS", method, c(-0.05, 0.45))
comp_SIS_IBS_suppl

ggsave(comp_SIS_IBS_suppl, filename = "pdf/suppl_SIS_IBS.pdf")

