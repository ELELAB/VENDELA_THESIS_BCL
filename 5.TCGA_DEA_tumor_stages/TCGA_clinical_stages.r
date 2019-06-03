library(SummarizedExperiment)
library(TCGAbiolinks)

source("3.TCGA_DEA/TCGA_DEA_functions.r")
source("4.TCGA_DEA_subtypes/TCGA_subtypes_functions.r")

setwd("5.TCGA_DEA_tumor_stages/")
#-------------------------------Define sample types--------------------------------------------

# List of cancer types. For example: c("BRCA","GBM","LUAD", ...)
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

# List of HUGO gene names to sort out for plots
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

# Colours to use in barplot. Must be same length as 'geneNames'
colours_genes <- c("red", "orange", "yellow", "#33FF33","#009999","#0099FF", "#333399", "purple", "pink")

for(cancer in cancer_types){
  project <- paste("TCGA-", cancer, sep = "")
  
  #Get clinical information from '1.TCGA_download_prepare/clinical'
  dataClin <- get(load(paste0("../1.TCGA_download_prepare/clinical/",cancer,"_clinical.rda")))

  # Get tumor_stages  
  tumor_stages <- unique(dataClin$tumor_stage)
  
  # Remove "not reported"
  tumor_stages <- tumor_stages[- which(tumor_stages == "not reported")]
  # For every tumor stage in the cancer type, do...
  for(stage in tumor_stages){
    
    # Sort out samples that belong to that stage from the clinical information
    dataClin_sub <- dataClin[which(dataClin$tumor_stage == stage),]
    
    # Don't continue if there is less than 10 samples and sort out samples that belong to that stage
    if(nrow(dataClin_sub) > 9){
      barcodes <- dataClin_sub$submitter_id
      barcodes_long <- get_full_barcodes(cancer, barcodes, tumor_sample, normal_sample)
      NT <- TCGAquery_SampleTypes(barcodes_long, "NT")
      TP <- TCGAquery_SampleTypes(barcodes_long, "TP")
      setwd("../5.TCGA_DEA_tumor_stages")
      
      # Do not do DEA if there is less than 5 samples for each sample type 
      if(length(NT) > 4 && length(TP) > 4){
        
        # Get dataFilt from '1.TCGA_download_prepare/filtering'
        dataFilt <- get(load(paste0("../1.TCGA_download_prepare/filtering/",cancer,"_dataFilt.rda")))
        dataFilt_stage <- dataFilt[,which(colnames(dataFilt) %in% barcodes_long)]
        dataNT_stage <- sample_types(colnames(dataFilt_stage), dataFilt_stage, "NT")
        dataTP_stage <- sample_types(colnames(dataFilt_stage), dataFilt_stage, "TP")
        
        # Remove " " in stage name to make it more suitable for filenames
        stage_name <- gsub(stage, pattern = " ", replacement = "_")
        # Get a summary on both normal and primary tumor samples on 
        filename_NT_sum <- paste0(cancer, "_", stage_name, "_NT_summary.txt")
        filename_TP_sum <- paste0(cancer, "_", stage_name, "_TP_summary.txt")
        NT_summary <- do_summary(dataNT_stage, normal_short, filename_NT_sum)
        TP_summary <- do_summary(dataTP_stage,tumor_short, filename_TP_sum)
        
        # Perform a differential expression analysis with limma+voom
        # filename_DEA = name on output file from DEA
        filename_DEA <- paste0(cancer, "_", stage_name, "_dataDEGs.rda")
        dataDEGs_stage <- get_DEGs(dataNT_stage, dataTP_stage, filename_DEA)
        
        # Sort out the genes in geneNames
        dataDEGs_stage_genes <- dataDEGs_stage[which(rownames(dataDEGs_stage) %in% geneNames),]
        save(dataDEGs_stage_genes, file = paste0("DEA/",cancer,"_",stage_name, "_dataDEGs_genes.rda"))
        
        # Do a barplot on the results from the DEA 
        # filename_barplot = name on .pdf barplot file
        filename_barplot <- paste0(cancer,"_",stage_name,"_DEGs_barplot.pdf")
        DEA_subtype_barplot(dataFilt_stage, stage_name, dataDEGs_stage_genes, colours_genes, geneNames, filename_barplot)
       
      }
    }
  }
}  
