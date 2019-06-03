library(SummarizedExperiment)
library(TCGAbiolinks)
library(GenomicDataCommons)
library(biomaRt)

setwd("4.TCGA_DEA_subtypes")
source("TCGA_subtypes_functions.r")

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

for(cancer in cancer_types){
  
  #Get subtype information 
  tabPanCancer <- as.data.frame(PanCancerAtlas_subtypes())
  cancer_subs <- tabPanCancer[which(tabPanCancer$cancer.type == cancer),]
  subtypes <- unique(cancer_subs$Subtype_Selected)
  
  # Subtypes to not include in the analysis
  donts <- c("NotAssigned", "Normal", paste0(cancer,"_LGG.NA"), paste0(cancer, ".NA"))
  donts_len <- which(subtypes %in% donts)
  if(length(donts_len) > 0){
    subtypes <- subtypes[- which(subtypes %in% donts)]
  }
  for(subtype in subtypes){
    barcodes <- cancer_subs[which(cancer_subs$Subtype_Selected == subtype),]$pan.samplesID
    if(length(barcodes) > 4){
      barcodes_long <- get_full_barcodes(cancer, barcodes, tumor_sample, normal_sample)
      NT <- TCGAquery_SampleTypes(barcodes_long, "NT")
      TP <- TCGAquery_SampleTypes(barcodes_long, "TP")
      if(length(NT) > 4 && length(TP) > 4){
        # Get dataFilt from '1.TCGA_download_prepare/filtering'
        dataFilt <- get(load(paste0("../1.TCGA_download_prepare/filtering/",cancer,"_dataFilt.rda")))
        dataFilt_sub <- dataFilt[,which(colnames(dataFilt) %in% barcodes_long)]
        dataNT_sub <- sample_types(colnames(dataFilt_sub), dataFilt, "NT")
        dataTP_sub <- sample_types(colnames(dataFilt_sub), dataFilt, "TP")
        
        # Get a summary on both normal and primary tumor samples on 
        filename_NT_sum <- paste0(cancer, "_", subtype, "_NT_summary.txt")
        filename_TP_sum <- paste0(cancer, "_", subtype, "_TP_summary.txt")
        NT_summary <- do_summary(dataNT_sub, "NT", filename_NT_sum)
        TP_summary <- do_summary(dataTP_sub, "TP", filename_TP_sum)
        
        # Perform a differential expression analysis with limma+voom
        # filename_DEA = name on output file from DEA
        filename_DEA <- paste0(cancer, "_", subtype, "_dataDEGs.rda")
        dataDEGs_sub <- get_DEGs(dataNT_sub, dataTP_sub, filename_DEA)
        
        # Sort out the genes in geneNames
        dataDEGs_sub_genes <- dataDEGs_sub[which(rownames(dataDEGs_sub) %in% geneNames),]
        save(dataDEGs_sub_genes, file = paste0("DEA/",cancer,"_",subtype, "_dataDEGs_genes.rda"))
        
        # Do a barplot on the results from DEA
        # filename_barplot = name on barplot .pdf file
        filename_barplot <- paste0(cancer,"_",subtype,"_DEGs_barplot.pdf")
        DEA_subtype_barplot(dataFilt_sub, subtype, dataDEGs_sub_genes, colours_genes, geneNames, filename_barplot)
      } 
    }
  }
}  
