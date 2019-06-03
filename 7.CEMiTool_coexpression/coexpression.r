#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("CEMiTool")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(CEMiTool)
library(ggplot2)

setwd("7.CEMiTool_coexpression/")

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


for(cancer in cancer_types){
  
  dataFilt <- get(load(paste0("../1.TCGA_download_prepare/filtering/",cancer, "_dataFilt.rda")))
  tum_pos <- c()
  tumor_barcodes <- TCGAquery_SampleTypes(colnames(dataFilt), "TP")
  dataFilt_tp <- dataFilt[, which(colnames(dataFilt) %in% tumor_barcodes)]
  
  
  # Do a limma + voom transformation and get results as log2-counts per million (logCPM)
  dataFilt_voom = voom(dataFilt_tp)$E
  dataFilt_ok <- as.data.frame(dataFilt_voom)
  
  # Find coexpression in that cancer type and divide into hubs and modules
  cem <- cemitool(dataFilt_ok)
  cem #to have basics information on the modules and selected genes
  filter_expr(cem)
  nmodules(cem)
  head(module_genes(cem))
  modules<- module_genes(cem)
  hubs <- get_hubs(cem,10)
  hubs_M1 <- as.data.frame(hubs[["M1"]])
  hubs_M2 <- as.data.frame(hubs[["M2"]])
  hubs_M3 <- as.data.frame(hubs[["M3"]])
  hubs_M4 <- as.data.frame(hubs[["M4"]])
  summary <- mod_summary(cem)
  
  # Plot results
  cem <- plot_profile(cem)
  plots <- show_plot(cem, "profile")
  plots[1]
  gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
  gmt_in <- read_gmt(gmt_fname)
  
  # perform over representation analysis
  cem <- mod_ora(cem, gmt_in)
  # plot ora results - in the example for module 1
  cem <- plot_ora(cem)
  plots <- show_plot(cem, "ora")
  plots[1]
  
  # read interactions with the CEMiTool default 
  int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
  int_df <- read.delim(int_fname)
  head(int_df)
  
  # plot interactions - in the example for module 1
  interactions_data(cem) <- int_df # add interactions
  cem <- plot_interactions(cem) # generate plot
  plots <- show_plot(cem, "interaction") # view the plot for the first module
  plots[1]
  
  # read interactions with i2d 
  i2d <- read.table("i2d_geneName.txt")
  i2d_df <- as.data.frame(i2d)
  
  # plot interactions with i2d - in the example for module 1
  
  interactions_data(cem) <- i2d_df # add interactions
  cem <- plot_interactions(cem) # generate plot
  plots <- show_plot(cem, "interaction") # view the plot for the first module
  plots[1]
  
  cat("saving files")
  generate_report(cem, output_format= "pdf_document", force=TRUE)
  
  # write analysis results into files
  write_files(cem, force=TRUE)
  
  # save all plots
  save_plots(cem, "all")
  
}

