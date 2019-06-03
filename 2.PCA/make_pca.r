library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)


setwd("2.PCA")
source("pca_functions.r")


# List of cancer types which can NOT filter for tumor purity. Foe example: c("STAD", "ESCA", ...)
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

for(cancer in cancer_types){
  dataFilt <- get(load(paste0("../1.TCGA_download_prepare/filtering/",cancer,"_dataFilt.rda")))
  
  # Do PCA
  # filename_pca = file name on .rda outpur from pca
  filename_pca <- paste0(cancer, "_pca.rda")
  dataFilt_pca <- get_pca(dataFilt, filename_pca)
  sample_types <- get_sample_types(dataFilt)
  
  
  
  # Plot variations
  #summary_file = file name on summary file
  #cumulative_plotfile = file name on plotted cumulative summary .pdf file
  #variance_plotfile = file name on plotted variance .pdf file
  #pca_plotfile = file name on PCA .pdf plot
  

  summary_file <- paste0(cancer, "_variations_summary.rda")
  cumulative_plotfile <- paste0(cancer, "_cumulative_summary.pdf")
  variance_plotfile <- paste0(cancer, "_PC_variance.pdf")
  pca_plotfile <- paste0(cancer, "_PCA_plot.pdf")
  plot_variations(dataFilt_pca, sample_types, summary_file, variance_plotfile, cumulative_plotfile, pca_plotfile)
}



