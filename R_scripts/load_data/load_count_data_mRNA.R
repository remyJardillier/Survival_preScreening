#  --- Working directories ---
wd_cancer <- paste0("../data_cancer/", cancer, "/")

count_mRNA <- t(read.csv(file = paste0(wd_cancer, cancer, ".uncv2.mRNAseq_raw_counts.txt"),
                    stringsAsFactors = F, sep = "\t", header =F))
count_mRNA[1:5, 1:5]

# colnames, rownames, NA values
colnames(count_mRNA) <- count_mRNA[1,]
count_mRNA <- count_mRNA[-1,]
row.names(count_mRNA) <-  count_mRNA[,1]
count_mRNA <- count_mRNA[,-1]
count_mRNA <- as.matrix(count_mRNA)
class(count_mRNA) <- "numeric"

count_mRNA[1:5, 1:5]

# keep only primary solid tumor samples 
# (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)
freq_sample_types <- as.data.frame(table(substr(row.names(count_mRNA), 14, 15)))
freq_sample_types[,2] <- paste0("(", freq_sample_types[,2], ")") 
print(paste("Sample types mRNA:", paste(paste(freq_sample_types[,1], freq_sample_types[,2]), collapse = ", ")))

if(cancer == "LAML"){
  id_prim_tumor <- as.numeric(substr(row.names(count_mRNA), 14, 15)) == 3
}else{
  id_prim_tumor <- as.numeric(substr(row.names(count_mRNA), 14, 15)) == 1
}

count_mRNA <- count_mRNA[id_prim_tumor, ]
row.names(count_mRNA) <- substr(row.names(count_mRNA),1,12)
dim(count_mRNA)
count_mRNA[1:5, 1:5]

