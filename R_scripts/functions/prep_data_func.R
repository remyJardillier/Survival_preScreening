

# logCPM normalization ----------------------------------------------------

# @description: logCPM normalization of a data frame
#
# @parameters: 
#   - valid.count: count data of RNA-seq values
#
# @return:
#   - logCPM data

log.cpm <- function(valid.count){

  # keep the genes that are expressed in at least 0.1% of samples
  valid.count <- t(valid.count)
  vc.dge <- DGEList(counts = valid.count)
  vc.dge.isexpr <- rowSums(cpm(vc.dge) > 1) >= round(dim(vc.dge)[2]*0.01)
  vc.dge <- vc.dge[vc.dge.isexpr,]
  
  # apply CPM normalization with normalization factor
  vc.dge <- calcNormFactors(vc.dge)
  vc.voom <- voom(vc.dge)
  vlc <- t(vc.voom$E)
  vlc <- vlc[complete.cases(vlc),]
  return(vlc)
  # return(log2(t(cpm(vc.dge, prior.count = 1)+1)))
}


# @description: logCPM normalization of a training and testing datasets
#
# @parameters: 
#   - count_train: count data of RNA-seq values (training dataset)
#   - count_test: count data of RNA-seq values (testing dataset)
#
# @return:
#   - list of logCPM data (training and testing)

log.cpm.cv <- function(count_train, count_test){
  
  # keep the genes that are expressed in at least 0.1% of samples
  vc.dge_train <- DGEList(counts = t(count_train))
  vc.dge.isexpr <- rowSums(cpm(vc.dge_train) > 1) >= round(dim(vc.dge_train)[2]*0.01)
  
  # reference column in the training dataset
  q75_libsize <- quantile(vc.dge_train$samples$lib.size, probs = 0.75)
  id_ref <- which.min(abs(vc.dge_train$samples$lib.size - q75_libsize))
  
  # compute the normalization factor for all patients
  count_all <- rbind(count_train, count_test)
  vc.dge_all <- DGEList(counts = t(count_all))
  vc.dge_all <- calcNormFactors(vc.dge_all, refColumn = id_ref)
  
  # remove genes that are not expressed in at least 1% of 
  # the patients of training dataset
  vc.dge_all <- vc.dge_all[vc.dge.isexpr,]
  
  # apply CPM normalization with normalization factor
  vc.voom_all <- voom(vc.dge_all)
  vlc_all <- t(vc.voom_all$E)

  return(list(logCPM_train = vlc_all[1:nrow(count_train),], 
              logCPM_test = vlc_all[-(1:nrow(count_train)),]))
}



# Standardization ---------------------------------------------------------

# @description: standardization of the training dataset
#
# @parameters: 
#   - data_train: training data to standardize
#
# @return:
#   - standardize data

std_train <- function(data_train){
  
  mean_train <- apply(data_train, 2, mean)
  sd_train <- apply(data_train, 2, sd)
  data_train_std <- sweep(data_train, 2, mean_train, FUN="-")
  data_train_std <- sweep(data_train_std,2, sd_train, FUN="/")
  
  return(data_train_std)
}

# @description: standardization of the training and testing datasets
#
# @parameters: 
#   - data_train: training data to standardize
#   - data_test: training data to standardize
#
# @return:
#   - standardize data of testing dataset
std_test <- function(data_train, data_test){
  
  mean_train <- apply(data_train, 2, mean)
  sd_train <- apply(data_train, 2, sd)
  
  data_test_std <- sweep(data_test, 2, mean_train, FUN="-")
  data_test_std <- sweep(data_test_std,2, sd_train, FUN="/")
  
  return(data_test_std)
}
