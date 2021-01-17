

# Multivariate Cox models -------------------------------------------------

# @description: learn models for elastic net, lasso, ridge and adaptive elastic net
#
# @parameters:
#   - gene_data: matrix of RNA-seq data (column: gene, row: patients)
#   - clin: data frame containing follow-up times and status
#   - method: "EN" for elastic net, "lasso" for lasso, "ridge" for ridge and "AEN"
#             for adaptive elastic net
#   
# @return:
#   - 'cv.glmnet' object
learn_models <- function(gene_data, clin, method){
  
  y_cox <- Surv(clin$time, clin$status)
  
  # Elastic Net
  if(method == "EN"){
    
    fit <- cv.glmnet(gene_data, y_cox, family = "cox",
                                alpha = 0.3, nfolds = 5, grouped = T, standardize = F)
    
    # lasso
  }else if(method == "lasso"){
    
    fit <- cv.glmnet(gene_data, y_cox, family = "cox",
                                alpha = 1, nfolds = 5, grouped = T, standardize = F)
      
    # ridge
  }else if(method == "ridge"){
    
    fit <- cv.glmnet(gene_data, y_cox, family = "cox",
                                alpha = 0, nfolds = 5, grouped = T, standardize = F)
    
  }else if(method == "coxph"){
    
    clin_gene_data <- cbind(y_cox, gene_data)
    fit <- coxph(y_cox~., data = clin_gene_data)
    
  # Adaptive Elastic Net
  }else if(method == "AEN"){
    
      # penalty terms
      fit_ridge <- cv.glmnet(gene_data, y_cox, family = "cox",
                             alpha = 0, nfolds = 5, grouped = T, standardize = F)
      betas_ridge <- as.numeric(coef(fit_ridge, "lambda.min"))
      
      w_AEN <- abs(1 / betas_ridge)
      
      # apply the method to the new data set
      fit <- cv.glmnet(gene_data, y_cox, family = "cox", 
                                nfolds = 5, standardize = F, alpha = 0.3, 
                                grouped = T, penalty.factor = w_AEN) 
  # sPLS
  }else if(method == "sPLS"){
    
    fit <- coxplsDR(Xplan = gene_data, time = clin$time, 
                              time2 = clin$status, scaleX = TRUE, scaleY=FALSE,
                              ncomp=20, allres=TRUE)
  }else if(method == "RF"){
    
    data_RF <- cbind(clin, gene_data)
    
    tune_fit <- tune.rfsrc(Surv(time, status) ~ ., data=data_RF,
                           mtryStart = floor(sqrt(ncol(data_RF) -1)),
                           ntreeTry = 50,
                           nsplit = 10, stepFactor = 2, improve =  0.02, strikeout = 3, maxIter = 25,
                           trace = FALSE, doBest = F)
    fit <- rfsrc(Surv(time, status) ~ ., data=data_RF, mtry = tune_fit$optimal["mtry"], 
                 nodesize = tune_fit$optimal["nodesize"], )
    
    
  }else if(method == "CoxBoost"){
    
    cv_boost <- cv.CoxBoost(time=clin$time,status=clin$status,
                             x=gene_data, 
                             penalty=9*sum(clin$status==1),
                             maxstepno=100, K=10,type="verweij")
    
    fit <- CoxBoost(time=clin$time,status=clin$status,
                      x=gene_data,stepno=cv_boost$optimal.step)
  }
  
  print("Model learned")
  return(fit)
}

# @description: return a vector containing the names of the genes selected
#
# @parameters:
#   - fit: 'cv.glmnet' object
#   - lambda: weight of the penalty
#
# @return:
#   - vector containing the names of the genes selected
genes_selected_func <- function(fit, lambda){
  
  if(class(fit) == "cv.glmnet"){
    # coefficients and active coefficients
    coefs <- coef(fit, s = lambda)
    names_active_coefs <- row.names(coefs)[coefs@i+1]
    return(names_active_coefs)
  }else{
    return(NULL)
  }
  
  
  
}


# Univariate Cox selection ------------------------------------------------

# @description: univariate Cox for one feature
#
# @parameters:
#   - x: feature (one gene)
#   - y_cox: Surv object containing follow-up times and status
#
# @return:
#   - the p-value of the univariate Cox model
p_val_cox_func <- function(x, y_cox){
  
  res <- try(fit <- coxph(y_cox ~ x))
  
  if(!inherits(res, "try-error")){
    
    test <- summary(fit)
    return(test$logtest[3])
    
  }else{
    return(NA)
  }
  
}


# @description: univariate Cox for an entire matrix
#
# @parameters:
#   - gene_data: matrix of RNA-seq data (column: gene, row: patients)
#   - y_cox: Surv object containing follow-up times and status
#
# @return:
#   - the p-value of the univariate Cox model
p_val_univCox_func <- function(gene_data, y_cox){
  
  # compute the p-values
  p_val_univCox_vect <- apply(gene_data, 2, function(x) p_val_cox_func(x, y_cox)) 
  
  print("Model learned")
  
  return(p_val_univCox_vect)
}



# Elastic net -------------------------------------------------------------

# @description: learn models for elastic net, lasso, ridge and adaptive elastic net
#
# @parameters:
#   - gene_data: matrix of RNA-seq data (column: gene, row: patients)
#   - clin: data frame containing follow-up times and status
#   - alpha: weight in the elastic net penalizatinn between l1 and l2 norms
#   
# @return:
#   - 'cv.glmnet' object
learn_model_EN <- function(gene_data, clin, alpha){
  
  y_cox <- Surv(clin$time, clin$status)
  
  fit <- cv.glmnet(gene_data, y_cox, family = "cox",
                   alpha = alpha, nfolds = 5, grouped = T, standardize = F)
  
  print("Model learned")
  return(fit)
}