# @description: compute prediction metrics
#
# @parameters:
#   - fit: cv.glmnet object
#   - clin_train: clinical data of the training dataset
#   - clin_test: clinical data of the testing dataset
#   - data_train: gene epxression data of the training dataset
#   - data_test: gene epxression data of the testing dataset
#   - method: "EN" for elastic net, "lasso" for lasso, "ridge" for ridge and "AEN"
#             for adaptive elastic net
#   - lambda: lambda value of the penalized Cox model
#   
# @return:
#   - vector containing IBS, C-index, and p-value of a univariate Cox model.

C_IBS_pVal_func <- function(fit, clin_train, clin_test, data_train, 
                            data_test, method, lambda = "lambda.min"){
  
  y_cox_test <- Surv(clin_test$time, clin_test$status)
  y_cox_train <- Surv(clin_train$time, clin_train$status)
  
  if(method == "ridge" | method == "EN" | method == "lasso" | method == "AEN"){
    
    if(length(genes_selected_func(fit, lambda)) >= 1){
      
      # prognostic index
      beta <- as.numeric(coef(fit, lambda))
      names(beta) <- colnames(data_train) 
      PI_test <-  data_test %*% beta
      PI_train <- data_train %*% beta
      
    }else{
      return(c(NA, NA, NA))
      #return(c(NA, NA))
    }
    
  }else if(method == "sPLS"){
    
    object <- fit$plsDR_mod
    ncomp <- object$ncomp
    am <- object$loadings$X
    dim(am)
    bm <- object$loadings$Y
    cm <- object$mat.c
    sigma.X = fit$XplanScal
    means.X = fit$XplanCent
    #colnames(object$X) <- gsub("X", "", colnames(object$X))
    
    # testing data
    newdata <- as.matrix(data_test[, match(colnames(object$X), colnames(data_test))])
    ones <- matrix(rep(1, nrow(newdata)), ncol = 1)
    t.pred <- array(0, dim = c(nrow(newdata), ncomp))
    for (h in 1:ncomp) {
      W <- am[, 1:h] %*% solve(t(cm[, 1:h]) %*% am[, 1:h])
      t.pred[, h] <- scale(newdata, center = means.X[match(colnames(object$X), colnames(data_test))], 
                           scale = sigma.X[match(colnames(object$X), colnames(data_test))]) %*% W[, h]
    }
    
    PI_test <- t.pred %*% fit$cox_plsDR$coefficients
    PI_train <- as.matrix(fit$tt_plsDR) %*% fit$cox_plsDR$coefficients
    
    
  }else if(method == "CoxBoost"){
    
    PI_test <- as.vector(predict(fit, newdata = data_test, type = "lp"))
    PI_train <- as.vector(predict(fit, newdata = data_train, type = "lp"))
    
  }else if(method == "RF"){
    
    data_test <- cbind(clin_test, data_test)
    data_train <- cbind(clin_train, data_train)
    
    PI_test <- predict(fit, newdata = data_test)$predicted
    PI_train <- predict(fit, newdata = data_train)$predicted
    
  }else if(method == "coxph"){
    
    PI_test <- as.vector(predict(fit, newdata = data_test, type = "lp"))
    PI_train <- as.vector(predict(fit, newdata = data_train, type = "lp"))
  }
  
  # re-calibrate the model and test for PH asumption
  res <- try(fit_coxph <- coxph(y_cox_train ~ as.vector(PI_train)))
  # PH_test <- cox.zph(fit_coxph)
  # p_valPH <- PH_test$table[2, "p"]
  p_valPH <- 1
  
  if(p_valPH >= 0.00 & !inherits(res, "try-error")){
    
    # test if there is at least 1 event
    if(sum(clin_test$status == 1) != 0){
      PI_train <- PI_train * fit_coxph$coefficients
      PI_test <- PI_test * fit_coxph$coefficients
      
      # concordance
      C <- concordance.index(PI_test, surv.time = clin_test$time, 
                             surv.event = clin_test$status)$c.index 
      
      # # IBS
      # IBS <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25),
      #                type = "brier")$ierror
      res <- try(IBS <- sbrier.score2proba(data.frame(time = clin_train$time, event = clin_train$status, score = PI_train), 
                                           data.frame(time = clin_test$time, event = clin_test$status, score = PI_test), 
                                           method = "cox")$bsc.integrated)
      
      # p-value of the univariate Cox
      p_val <- p_val_cox_func(PI_test, y_cox_test)
      
      if(!inherits(res, "try-error")){
        return(c(C=C, IBS=IBS, p_val))
      }else{
        return(c(C=C, IBS=NA, p_val))
      }
      
    }else{
      return(c(C=NA, IBS=NA, pvalue = NA))
    }
    
    
  }else{
    return(c(C=NA, IBS=NA, pvalue = NA))
  }
  
  
}
