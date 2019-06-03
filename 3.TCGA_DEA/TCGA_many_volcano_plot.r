
library(SummarizedExperiment)
library(TCGAbiolinks)

#---------------------------------Define cancer types----------------------------------


# List of cancer types. For example: c("BRCA", "GBM", "LUAD", ...)
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
#-------------------------------------Do volcano plot----------------------------------------------------

for(cancer in cancer_types){
  dir.create("plots")
  #Get dataDEGs files 
  dataDEGs <- get(load(paste(cancer,"_dataDEGs.rda", sep = "")))
  dataDEGs_genes <- dataDEGs[which(rownames(dataDEGs) %in% geneNames),]
  
  # Do a volcano plot
  TCGAVisualize_volcano(
    x = dataDEGs$logFC,
    y = dataDEGs$adj.P.Val,
    x.cut = 1,
    y.cut = -log10(0.05),
    names = rownames(dataDEGs_genes),
    color = c("black","red","darkgreen"),
    names.size = 3,
    xlab = " Gene expression fold change (Log2)",
    legend = "State",
    title = paste(cancer," Volcano plot"),
    show.names = "both",
    highlight = rownames(dataDEGs_genes),
    filename = paste("plots/"cancer,"_volcano.pdf", sep = ""))
}
