This file contains instructions on how to use the code associated to the paper "Benchmark of lasso-like penalized Cox models reveals large differences between cancers and improved performance with two-dimensional pre-screening".


How to download the TCGA data (optional, a sample with 100 patients is provided)
==================================

To overcome runing time and size issues, the data available in 'Rdata_cancer_screen/' contain only 100 patients for 3 cancers (ACC, BRCA and KIRC) in .RData files.

Clinical and mRNA-seq datasets were obtained using the Broad GDAC FIREHOSE utility.

1. Go on the website https://gdac.broadinstitute.org
2. Choose the cancer and click on 'Browse' in the 'Data' column
3. Clinical data:
	3.a. On the 'Clinical' panel, click on 'Merge_Clinical'. A .zip file is downloaded.
	3.b. Open this .zip file and copy the file 'cancer-name'.clin.merged.txt' 
	(e.g. 'KIRC.clin.merged.txt') into the folder '../data_cancer/'cancer-name'/' (e.g. '../data_cancer/KIRC/').
4. mRNA data:
	4.a. On the 'mRNASeq' panel, click on 'mRNAseq_Preprocess'. A zip file is downloaded.
	4.b. Open this .zip file and copy the file ''cancer-name'.uncv2.mRNAseq_RSEM_all.txt' 
	(e.g. 'KIRC.uncv2.mRNAseq_RSEM_all.txt') into the folder '../data_cancer/'cancer-name'/' (e.g. '../data_cancer/KIRC/').	

In this script, we use RData files that we created to load RNA-seq and clinical data. To use TCGA data downloaded as mentionned above, please use the commented part of the file 'R_scripts/load_data/load_data_final.R'.

How to run the code 
==================================

1. Open the file 'R_scripts/Survival_preScreening.Rproj' to open RStudio environment.
2. Open the 'R_scripts/main.R' file and follow the instruction in the script. 

As the running time is high, results will be saved all along the code in 'R_scripts/data_fit' into .RData files. Please set 'learn_new_models' variable to 'TRUE' if you want to learn new models. Otherwise, previous models (if they exist) will be loaded.


Running times 
==================================

The article value for the number of models to be learned is 50 (10 repetitions of 5-fold cross-validation). It can take a lot of times with 6 IQR thresholds, 6 p-value thresholds, and 50 models learned (typically 15 hours), and we advice the user to first set K_folds = 3 (3 folds) and n_rep = 2 (2 repetitions)) for a first check.
