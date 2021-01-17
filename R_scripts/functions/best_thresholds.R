
# Optimal thresholds ------------------------------------------------------


# @description: compute the optimal thresholds
#
# @parameters: 
#   - ary: array of the C-indices computed with the function 'optimize_pre_filtering'
#   
#
# @return: 
#   - thrs_p_val_opt: optimal supervised threshold
#   - thrs_IQR_opt: optimal unsupervised threshold
#   - obj: highest median C-index

best_thrs <- function(ary, max_min){
  
  # remove data with more than 10 errors
  id_no_NA <- apply(ary, c(1,2), function(x) sum(is.na(x)) < 15)
  
  # compute the medians
  med_ary <- matrix(ncol = ncol(id_no_NA), nrow = nrow(id_no_NA))
  colnames(med_ary) <- colnames(id_no_NA)
  row.names(med_ary) <- row.names(id_no_NA)
  for(i in 1:nrow(id_no_NA)){
    for(j in 1:ncol(id_no_NA)){
      if(id_no_NA[i,j])
        med_ary[i,j] <- median(ary[i,j,], na.rm = T)
    }
  }
   
  
  # melted data frame
  med_ary <- melt(med_ary)
  if(sum(is.na(med_ary$value)) == nrow(med_ary)){
    return(c(thrs_p_val_opt = 1,
             thrs_IQR_opt = 0,
             obj = NA))
  }else if(max_min == "max"){
    id_opt <- which.max(med_ary$value)
  }else if(max_min == "min"){
    id_opt <- which.min(med_ary$value)
  }
  
  
  return(c(thrs_p_val_opt = med_ary[id_opt,1],
              thrs_IQR_opt = med_ary[id_opt,2],
              obj = med_ary[id_opt,3]))
  
}


# Other functions ---------------------------------------------------------

# @description: compute the inices of the n lowest number in a vector
#               used to compute the indices of the n lowest p-values
#
# @parameters: 
#   - x: a numeric vector
#   - n: the number of indices with lowest values in x to return
#
# @return:
#   - indices of the n lowest number in the vector x
which.minn <- function(x,n=1){
  if (n==1)
    which.min(x)
  else
  {
    if (n>1){
      ii <- order(x,decreasing=FALSE)[1:min(n,length(x))]
      ii[!is.na(x[ii])]
    }
    else {
      stop("n must be >=1")
    }
  }
}

# @description: compute the inices of the n highest number in a vector
#               used to compute the indices of the n lowest p-values
#
# @parameters: 
#   - x: a numeric vector
#   - n: the number of indices with highest values in x to return
#
# @return:
#   - indices of the n lowest number in the vector x
which.maxn <- function(x,n=1){
  if (n==1)
    which.max(x)
  else
  {
    if (n>1){
      ii <- order(x,decreasing=T)[1:min(n,length(x))]
      ii[!is.na(x[ii])]
    }
    else {
      stop("n must be >=1")
    }
  }
}





