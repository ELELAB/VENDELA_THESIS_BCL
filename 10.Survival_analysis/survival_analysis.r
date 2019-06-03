#### survival analysis
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(plyr)
library(ggfortify)
library(ggplot2)
library(survminer)


setwd("10.Survival_analysis")
source("survival_analysis_functions.r")
#source("/data/user/vendela/case_studies/code/data_curation/TCGAbiolinks_functions.R")


# TCGA cancer types
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
                  "CHOL"
)

# HUGO gene names
geneNames <- c("BCL2",
               "BCL2L1", #Bcl-xL
               "BCL2L2", #Bcl-w
               "MCL1",
               "BCL2A1", #Bfl-1
               "BCL2L10",
               "BAX",
               "BAK1",
               "BOK"
)


for(cancer in cancer_types){

  #tumor_rda <- get(load(paste("/data/user/vendela/case_studies/data/",tolower(cancer),".exp.rda", sep = "")))
  dataFilt <- get(load(paste0("../1.TCGA_download_prepare/filtering/",cancer, "_dataFilt.rda")))
  dataFilt_mat <- as.matrix(dataFilt)
  for (gene in geneNames){
    gene_pos <- which(rownames(dataFilt) == gene)
    if(length(gene_pos)>0){

      clinical <- get(load(paste0("../1.TCGA_download_prepare/clinical/",cancer,"_clinical.rda")))
      clinical_survival <- get_survival_table(dataFilt_mat,clinical,gene, 0.25,0.75,"KM")
      survival_plot(clinical_survival,gene,cancer,25)
      p.val <- p_func(clinical_survival)
    }
  }
}

