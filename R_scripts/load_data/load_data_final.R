load(file = paste0("../RData_cancer_screen/", cancer, ".RData"))

# keep patients that are both in clinical and count data
com_pat <- intersect(row.names(clin), row.names(count_mRNA))
com_pat <- intersect(com_pat, row.names(clin_data_cox))
clin <- clin[com_pat, ]
count_mRNA <- count_mRNA[com_pat,]
clin_data_cox <- clin_data_cox[com_pat,]
length(com_pat)

# remove genes with NA values 
id_NA_gene <- apply(count_mRNA, 2, function(x) sum(which(is.na(x))))
id_NA_gene <- id_NA_gene > 0
sum(id_NA_gene)
count_mRNA <- count_mRNA[, !id_NA_gene]
print(paste0(sum(id_NA_gene), " gene(s) removed due to NA values"))

# Surv object for the Cox model
y_cox <- Surv(time = clin$time, event = clin$status)