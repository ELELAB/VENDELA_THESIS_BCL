library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("3.TCGA_DEA")

# Get functions to perform a differential expression analysis
source("TCGA_DEA_functions.r")

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

# HUGO gene names to focus on 
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

# Colours to use in barplot by cancer types. Must be same length as 'geneNames'
colours_genes <- c("red", "orange", "yellow", "#33FF33","#009999","#0099FF", "#333399", "purple", "pink")

# Colours to use in barplot by gene name. Must be same length as 'cancer_types'
colours_cancer <- c("#FF0000FF", "#FF5500FF", "#FFAA00FF", "#FFFF00FF", "#AAFF00FF", "#55FF00FF", "#00FF00FF", "#00FF55FF", "#00FFAAFF",
                    "#00FFFFFF", "#00AAFFFF", "#0055FFFF", "#0000FFFF", "#5500FFFF", "#AA00FFFF", "#FF00FFFF", "#FF00AAFF", "#FF0055FF")

for(cancer in cancer_types){
  
  dataFilt <- get(load(paste0("../1.TCGA_download_prepare/filtering/",cancer, "_dataFilt.rda")))
  
  # Divide dataFilt depending on sample type (normal or primary tumor samples)
  dataN <- sample_types(colnames(dataFilt), dataFilt,"NT")
  dataT <- sample_types(colnames(dataFilt), dataFilt, "TP")
  
  # Get a summary on both normal and primary tumor samples on 
  filename_NT_summary <- paste0(cancer, "_NT_summary.txt")
  filename_TP_summary <- paste0(cancer, "_TP_summary.txt")
  NT_summary <- do_summary(dataN, "NT", filename_NT_summary)
  TP_summary <- do_summary(dataT, "TP", filename_TP_summary)

  # Perform a differential expression analysis with limma+voom 
  # filename = name on output DEA file
  filename_DEA <- paste0(cancer, "_dataDEGs.rda")
  dataDEGs <- get_DEGs(dataN, dataT, filename_DEA)
  # Sort out the genes in geneNames
  dataDEGs_genes <- dataDEGs[which(rownames(dataDEGs) %in% geneNames),]
  save(dataDEGs_genes, file = paste0("DEA/",cancer, "_dataDEGs_genes.rda"))
  
  # Do a barplot
  # filename_barplot = name on barplot .pdf file
  filename_barplot <- paste0(cancer, "_DEA_barplot.pdf")
  DEA_barplot(dataFilt, dataDEGs_genes, colours_genes, geneNames, filename_barplot)
  
  # Do volcano plot
  # filename_volcano = name on volcano plot .pdf name
  filename_volcano <- paste0(cancer, "_volcano_plot.pdf")
  DEA_volcanoplot(dataDEGs, dataDEGs_genes, filename_volcano)
}

# Make barplots by gene name in 'geneNames' list
# showing differential expression in all cancer types in 'cancer_types' -list
# Note that all DEA must be prepared for all cancer types in 'cancer_types'before running
plot_by_genes(cancer_types, geneNames, colours_cancer)
