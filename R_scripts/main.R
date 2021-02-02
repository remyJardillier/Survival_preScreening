# Load the packages needed ------------------------------------------------

source(file = "load_data/load_packages.R")

# load functions
invisible(lapply(list.files(path = "functions/", pattern = "[.]R$", 
                            recursive = TRUE, full.names = T), source))

# list of cancers
cancers_all <- list.files("../RData_cancer_screen",  )
cancers_all <- tools::file_path_sans_ext(cancers_all)

# create folders to save the object in .RData files for each cancer
methods_vect <- c("ridge", "EN")
for(cancer in cancers_all){
  if (!dir.exists(paste0("data_fit/", cancer)))
    dir.create(paste0("data_fit/", cancer), recursive = T)
  
  for(method in methods_vect){
    if (!dir.exists(paste0("data_fit/", cancer, "/", method)))
      dir.create(paste0("data_fit/", cancer, "/", method), recursive = T)
  }
}

# create folders for the figures, pdf, and tables
for(cancer in cancers_all){
  if (!dir.exists("Figures"))
    dir.create("Figures", recursive = T)
  
  if (!dir.exists("tables"))
    dir.create("tables", recursive = T)
  
  if (!dir.exists("pdf"))
    dir.create("pdf", recursive = T)
}

# C-index for all cancers -------------------------------------------------

# parameters for the fit
K_folds <- 3 # 5 in the paper
n_rep <- 2 # 10 in the paper
perc_min_clin <- 90

# compute the C-indices with mRNA data (ridge), clinical data, and both
learn_new_models <- T # choose if new models have to be learned
source(file = "1.1_pred_clin_mRNA.R")

# sort the cancers according to the median C-index computed with mRNA-seq data
med_C_all <- rep(NA, length(cancers_all))
names(med_C_all) <- cancers_all

for(cancer in cancers_all){
  
  # load the data
  load(file = paste0("data_fit/", cancer, "/ridge/pred_clin_mRNA.RData"))
  
  # median C-index for ridge
  med_C_all[cancer] <- median(C_df_final$mRNA, na.rm = T)  
}

med_C_sorted <- sort(med_C_all, decreasing = T)
med_C_sorted

cancers_all <- names(med_C_sorted)

# Choose cancer -----------------------------------------------------------

p_val_test <- rep(NA, length(cancers_all))
names(p_val_test) <- cancers_all

for(cancer in cancers_all){
  
  # load the data
  load(file = paste0("data_fit/", cancer, "/ridge/pred_clin_mRNA.RData"))
  
  # median C-index for ridge
  p_val_test[cancer] <- wilcox.test(C_df_final$mRNA, mu = 0.6, 
                                    alternative = "greater")$p.value
}

# BH correction 
s_level <- 0.2 # 0.01 in the paper
p_val_test_BH <- p.adjust(p_val_test, method = "BH")
sum(p_val_test_BH <= s_level)
length(p_val_test_BH)
p_val_test_BH[p_val_test_BH <= s_level]
hist(p_val_test_BH)

# boxplot to choose cancers
source(file = "1.2_choose_cancers_plot.R")

# plot boxplots of C-indices with mRNA, clinical and both for all cancers choosed
source(file = "1.3_pred_clin_mRNA_plot.R")

# characteristics of the cancers choosed
cancers_all <- names(p_val_test_BH[p_val_test_BH <= s_level])
source(file = "2_ch_cancers.R")


# Learn models with different filtering thresholds ------------------------

# mean-variance trend for the IQR and Variance Stabilizatio Transformation (VST)
cancer <- "ACC"
source(file = "3.1_VST.R")

# parameter for the pre-filtering and the fit
thrs_p_val <- c(0.2,1)
thrs_IQR <- c(0, 1.5)

# !!! this step can take some times !!!
method <- "EN"
learn_new_models <- T # choose if new models have to be learned
source(file = "3.2_flt_learn_models.R")

method <- "ridge"
learn_new_models <- T # choose if new models have to be learned
source(file = "3.2_flt_learn_models.R")

# nested cross-validataion ---
# learn new models with filtering thresholds computed in the previous step

# learn models for nested CV
pred_measure <- "C" 
learn_new_models <- T # choose if new models have to be learned
source(file = "3.3_flt_nested_CV.R")

pred_measure <- "IBS" 
learn_new_models <- T # choose if new models have to be learned
source(file = "3.3_flt_nested_CV.R")

# produce the main figures ---
# cancer and prediction measure for the main figure
cancer <- "ACC"
type <- "C"

source(file = "3.4_raw_vs_flt_plot.R")


# Irrelevant genes selected without pre-screening -------------------------

# choose cancer and method
cancer <- "ACC"
method <- "EN"

source(file = "4_slc_irrelevant_genes.R")

# Benchmark of penalization methods after pre-screening -------------------

source(file = "5_benchmark.R")

# Comparison with ISIS ----------------------------------------------------

learn_new_models <- T # choose if new models have to be learned
source(file = "6_SIS.R")
