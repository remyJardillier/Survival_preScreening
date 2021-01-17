#download outcome data from Liu et al. Cell 2018
#freely accessible on PubMed Central: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/
upd.Surv <- read.xlsx("../data_cancer/surv_time/NIHMS978596-supplement-1.xlsx",
                      sheet = "TCGA-CDR") 
#clean up data
upd.Surv <- upd.Surv[,-1]
upd.Surv$OS <- as.character(upd.Surv$OS) %>%
  as.numeric()
upd.Surv$OS.time <- as.character(upd.Surv$OS.time) %>%
  as.numeric()
upd.Surv$PFI <- as.character(upd.Surv$PFI) %>%
  as.numeric()
upd.Surv$PFI.time <- as.character(upd.Surv$PFI.time) %>%
  as.numeric()

#keep only data for cancer type analyzed here
clin <- upd.Surv[upd.Surv$type == cancer
                 ,c("bcr_patient_barcode",
                    "OS", "OS.time",
                    "PFI", "PFI.time")] %>%
  droplevels(.)
rownames(clin) <- clin$bcr_patient_barcode
clin <- clin[, -1]
rm(upd.Surv)

clin[1:5,]

# choose PFI or OS
if(type == "PFI") {
  clin <- clin[, c("PFI", "PFI.time")]
  colnames(clin) <- c("status", "time")
} else {
  clin <- clin[, c("OS", "OS.time")]
  colnames(clin) <- c("status", "time")
}
clin <- clin[!is.na(clin$time),]
clin <- clin[clin$time > 0,]
clin$time <- clin$time / 365.25
