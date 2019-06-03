library(SummarizedExperiment)
library(TCGAbiolinks)

#------------------------ Functions to download and aggregate TCGA gene expression data -----------------------------------------------

# Download samples
get_samples <- function(project, tumor_sample, normal_sample){
  dir.create("data", showWarnings = FALSE)
  setwd("data")
  cat("Downloading query. Cancer type =", cancer, "\n")
  query.exp <- GDCquery(project = project,
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - Counts",
                        sample.type = c(tumor_sample, normal_sample),
                        legacy = FALSE)
  GDCdownload(query.exp)
  setwd("..")
  return(query.exp)
  
}

# Download .rda files from GDC
get_GDC <- function(query.exp, tumor_rda){
  cat("Downloading .rda files. Cancer type =", cancer, "\n")
  setwd("data")
  tumor_rd <- GDCprepare(query = query.exp, save = TRUE, save.filename = tumor_rda)
  setwd("..")
  return(tumor_rd)
  
}

# Sort out females from BRCA
female <- function(gdc_brca, tumor_rda){
  cat("Sorting out Females from ",cancer, "\n")
  setwd("data")
  female <- subset(x = gdc_brca, select = (gender == "female"))
  save(female , file = tumor_rda)
  setwd("..")
  return(female)
  
}