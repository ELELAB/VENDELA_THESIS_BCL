library(SummarizedExperiment)
library(TCGAbiolinks)
setwd("1.TCGA_download_prepare")
# Get functions from "TCGA_download_functions.r"
source("TCGA_download_functions.r")

#-------------------------------Define sample types--------------------------------------------

# Choose tumor sample from project (long name). For example: "Primary solid Tumor"
tumor_sample <- "Primary solid Tumor"

# Short name for tumor sample (must match chosen tumor sample). 
# For example (if tumor sample is "Primary solid Tumor"): "TP" 
tumor_short <- "TP"

# Choose normal sample from project (long name). For example: "Solid Tissue Normal
normal_sample <- "Solid Tissue Normal"

# Short name for normal sample (must match chosen normal sample). 
# For example (if normal sample is "Solid Tissue Normal"): "NT" 
normal_short <- "NT"

# List of cancer types which can filter for tumor purity. For example: c("BRCA", "GBM", "LUAD", ...)
cancer_types <- c("BRCA",
                  "GBM",
                  "LUAD",
                  "UCEC",
                  "KIRC",
                  "HNSC",
                  "THCA",
                  "LUSC",
                  "PRAD",
                  "COAD",
                  "STAD",
                  "BLCA",
                  "LIHC",
                  "KIRP",
                  "ESCA",
                  "READ",
                  "KICH",
                  "CHOL")


#----------------For-loop to download and aggregate TCGA data from all cancer types in cancer_types-----------------------


for(cancer in cancer_types){

  project <- paste("TCGA-", cancer, sep = "")
  
  # Download TCGA data
  tumor_rda <- paste0(tolower(cancer), ".exp.rda")
  query.exp <- get_samples(project, tumor_sample, normal_sample)
  
  # Aggregate TCGA data. If cancer type = BRCA, keep the females only
  if(cancer == "BRCA"){
    gdc_brca <- get_GDC(query.exp, tumor_rda)
    gdc <- female(gdc_brca, tumor_rda)
  } else{
    gdc <- get_GDC(query.exp, tumor_rda)
  }
}  
