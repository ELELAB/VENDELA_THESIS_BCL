library(TCGAbiolinks)
library(ggplot2)
library(reshape2)
setwd("1.TCGA_download_prepare")
dir.create("plots", showWarnings = FALSE)

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

# Which sample types to find in expression data (short letter code)
sample_types_short <- c("TP", "NT")

for(cancer in cancer_types){
  
  # Get dataFilt from '1.TCGA_download_prepare/filtering'
  dataFilt <- get(load(paste0("filtering/",cancer,"_dataFilt.rda")))
  
  # Do a limma+voom transformation. 
  dataFilt_voom <- limma::voom(dataFilt)
  dataFilt_voom <- dataFilt_voom[[1]]
  
  # Get expression values for the genes in 'geneNames'
  dataFilt_voom_sub <- dataFilt_voom[which(rownames(dataFilt_voom) %in% geneNames),] 
  
  # Do a boxplot with sample types as variables
  df_t <- as.data.frame(t(dataFilt_voom_sub))
  df_t$sample_type <- rep("na", length(rownames(df_t)))
  for(sample in sample_types_short){
    dataFilt_sample <- dataFilt_voom_sub[,TCGAquery_SampleTypes(colnames(dataFilt_voom_sub), sample)]
    if(ncol(dataFilt_sample) > 1){
      df_t[colnames(dataFilt_sample),]$sample_type <- sample
    }
  }
  col_list<- c("blue", "gold")
  cols <- col_list[seq(length(unique(df_t$sample_type)))]
  dfm <- melt(df_t)
  ggplot(data=dfm) + 
    geom_boxplot(aes(x=sample_type, y = value),
                 colour = rep(cols,length(unique(dfm$variable)))) + 
    facet_wrap(~variable) +
    ggtitle(cancer) +
    xlab("Sample types") + 
    ylab("log2 gene expression")
  ggsave(paste0("plots/",cancer, "_boxplot.pdf"))
}
