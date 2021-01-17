
print(paste("*** start learning for", cancer, "***"))
wd_cancer <- paste0("../data_cancer/", cancer, "/")

# --- import and prepare clinical data ---
clin_data <- t(read.csv(file = paste0(wd_cancer, cancer, ".clin.merged.txt"),
                        stringsAsFactors = F, sep = "\t", header = F))

colnames(clin_data) <- clin_data[1,]
clin_data <- clin_data[-1,]
row.names(clin_data) <- toupper(clin_data[, "patient.bcr_patient_barcode"])
clin_data <- as.data.frame(clin_data)
n_patients <- nrow(clin_data)

# age
if("patient.primary_pathology.age_at_initial_pathologic_diagnosis" %in% colnames(clin_data)){
  
  age <- as.numeric(as.character(clin_data$patient.primary_pathology.age_at_initial_pathologic_diagnosis))
  
}else if("patient.age_at_initial_pathologic_diagnosis" %in% colnames(clin_data)){
  
  age <- as.numeric(as.character(clin_data$patient.age_at_initial_pathologic_diagnosis))
  
}

sum(is.na(age))
length(age)

# gender
gender <- factor(clin_data$patient.gender)
levels(gender)

if(length(levels(gender)) > 1){
  clin_data_cox <- data.frame(age = age, gender = gender)
}else{
  clin_data_cox <- data.frame(age = age)
}

# grade
if("patient.neoplasm_histologic_grade" %in% colnames(clin_data)){
  
  grade <- as.character(clin_data$patient.neoplasm_histologic_grade)
  id_NA_grade <- is.na(grade)
  
  if(sum(!id_NA_grade) / n_patients >= perc_min_clin/100){
    
    grade <- factor(grade)
    levels(grade)
    
    clin_data_cox <- cbind(clin_data_cox, grade = grade)
  }
  
}


# T
if("patient.stage_event.tnm_categories.pathologic_categories.pathologic_t" %in% colnames(clin_data)){
  
  id_NA_T <- is.na(clin_data$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t)
  
  if(sum(!id_NA_T) / n_patients >= perc_min_clin/100){
    
    T_stage <- clin_data$patient.stage_event.tnm_categories.pathologic_categories.pathologic_t
    T_stage <- factor(substr(T_stage, start = 1, stop = 2))
    levels(T_stage)
    
    clin_data_cox <- cbind(clin_data_cox, T_stage = T_stage)
  }
}

# N
if("patient.stage_event.tnm_categories.pathologic_categories.pathologic_n" %in% colnames(clin_data)){
  
  id_NA_N <- is.na(clin_data$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n)
  
  if(sum(!id_NA_N) / n_patients >= perc_min_clin/100){
    
    N_stage <- clin_data$patient.stage_event.tnm_categories.pathologic_categories.pathologic_n
    N_stage <- factor(substr(N_stage, start = 1, stop = 2))
    levels(N_stage)
    
    clin_data_cox <- cbind(clin_data_cox, N_stage = N_stage)
  }
}

# M
if("patient.stage_event.tnm_categories.pathologic_categories.pathologic_m" %in% colnames(clin_data)){
  
  id_NA_M <- is.na(clin_data$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m)
  
  if(sum(!id_NA_M) / n_patients >= perc_min_clin/100){
    
    M_stage <- clin_data$patient.stage_event.tnm_categories.pathologic_categories.pathologic_m
    M_stage <- factor(substr(M_stage, start = 1, stop = 2))
    levels(M_stage)
    
    clin_data_cox <- cbind(clin_data_cox, M_stage = M_stage)
  }
}


row.names(clin_data_cox) <- row.names(clin_data)
str(clin_data_cox)

# remove patients with NA data
id_NA_pat <- apply(clin_data_cox, 1, function(x) sum(is.na(x)) > 0)
sum(id_NA_pat)
clin_data_cox <- clin_data_cox[!id_NA_pat, ]
str(clin_data_cox)

# common patients with mRNA and miRNA RNA-seq data
com_patient <- intersect(row.names(count_mRNA), row.names(clin_data_cox))
length(com_patient)

count_mRNA <- count_mRNA[com_patient, ]
clin_data_cox <- clin_data_cox[com_patient, ]
clin <- clin[com_patient, ]
