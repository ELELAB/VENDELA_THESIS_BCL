
#----------- Functions to perform a differential expression analysis on filtered TCGA files ------------------

# Get sample types from filtered data
sample_types <- function(barcodes, dataFilt, sample_type_short){
  cat(paste("Finding ",sample_type_short, " sample types", "\n"))
  dataSm <- TCGAquery_SampleTypes(barcodes,sample_type_short)
  st <- which(dataSm %in% colnames(dataFilt))
  dataSm <- dataSm[st]
  dataST <- dataFilt[,dataSm]
  return(dataST)
}

# Get a summary for the tumor and normal samples
do_summary <- function(dataN, normal, filename){
  dir.create("summary", showWarnings = FALSE)
  dataN_sum <- t(dataN)
  ens_normal_sum <- summary(dataN_sum)
  ens_normal_sum <- t(ens_normal_sum)
  write.table(ens_normal_sum, paste0("summary/",filename), sep = "\t", quote = FALSE)
}

# Do a differential expression analysis

get_DEGs <- function(dataN, dataT, filename){
  cat(paste("Differential expression analysis. Cancer type = ", cancer, "\n"))
  dir.create("DEA", showWarnings = FALSE)
  dataDEGs <- TCGAanalyze_DEA(mat1 = dataN,
                              mat2 = dataT,
                              pipeline = "limma", 
                              batch.factors = "TSS",
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 0.05 ,
                              logFC.cut = 0,
                              voom = TRUE)
  save(dataDEGs, file= paste0("DEA/",filename))
  return(dataDEGs)
}

DEA_barplot <- function(dataFilt, dataDEGs_genes, colours_genes, geneNames, filename_barplot){
  dir.create("plots/genes", showWarnings = FALSE)
  
  # Make dataframe for barplot
  genes_order <- as.data.frame(geneNames)
  genes_order$geneNames <- as.character(genes_order$geneNames) 
  genes_order$logFC <- rep(0,nrow(genes_order))
  genes_order$colours <- as.character(colours_genes)
  for(row in rownames(dataDEGs_genes)){
    genes_order[which(genes_order$geneNames == row),]$logFC <- dataDEGs_genes[which(rownames(dataDEGs_genes) == row),]$logFC
  }
  
  # Add information about missing values
  not_deg <- as.character(genes_order[which(genes_order$logFC == 0),]$geneNames)
  not_filt <- as.character(genes_order[- which(genes_order$geneNames %in% rownames(dataFilt)),]$geneNames)
  nofilt_text <- ""
  nodeg_text <- ""
  if(length(not_filt) > 0){
    not_deg <- not_deg[-which(not_deg %in% not_filt)]
    nofilt_text <- "** Gene not present after filtering"
    for(nofilt in not_filt){
      genes_order[which(genes_order$geneNames == nofilt),]$geneNames <- paste0(nofilt, "**")
    }
    if(length(not_deg) > 0){
      for(nodeg in not_deg){
        genes_order[which(genes_order$geneNames == nodeg),]$geneNames <- paste0(nodeg, "*")
        nodeg_text <- "* Gene not differentially expressed"
      }
    }
  }else{
    if(length(not_deg) > 0){
      for(nodeg in not_deg){
        genes_order[which(genes_order$geneNames == nodeg),]$geneNames <- paste0(nodeg, "*")
        nodeg_text <- "* Gene not differentially expressed"
      }
    }
  }
  
  # Make barplot
  pdf(paste0("plots/", filename_barplot))
  bp <- barplot(genes_order$logFC,
                names.arg = genes_order$geneNames,
                ylim = c(-4,6),
                cex.axis=0.9,
                cex.names=0.7,
                main = paste(cancer, "differential expression"),
                ylab = "logFC",
                col = genes_order$colours,
                axes = TRUE)
  mtext(nodeg_text, side = 3, adj = 0)
  mtext(nofilt_text, side = 3, adj = 1) 
  dev.off()
  return(bp)
}

DEA_volcanoplot <- function(dataDEGs, dataDEGs_genes, filename_volcano){
  dir.create("plots", showWarnings = FALSE)

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
    filename = paste0("plots/", filename_volcano))
}

plot_by_genes <- function(cancer_types, geneNames, colours_cancer){
  dir.create("plots", showWarnings = FALSE)
  all_DEGs <- data.frame()
  
  # Make dataframe to use in barplot with information about all cancer types in 'cancer_types'-list
  for(cancer in cancer_types){
    
    # Which colour should that cancer type be represented by
    col <- colours_cancer[which(cancer_types == cancer)]
    
    # Get DEG file
    dataDEGs_genes <- get(load(paste0("DEA/",cancer, "_dataDEGs_genes.rda")))
    
    # Get filtered file from the pre-processing steps
    dataFilt <- get(load(paste0("../1.TCGA_download_prepare/filtering/",cancer, "_dataFilt.rda")))
    
    # Make dataframe
    genes_order <- as.data.frame(geneNames)
    genes_order$geneNames <- as.character(genes_order$geneNames)
    genes_order$cancer <- rep(cancer, nrow(genes_order))
    genes_order$logFC <- rep(0,nrow(genes_order))
    genes_order$colours <- as.character(rep(col, nrow(genes_order)))
    for(row in rownames(dataDEGs_genes)){
      genes_order[which(genes_order$geneNames == row),]$logFC <- dataDEGs_genes[which(rownames(dataDEGs_genes) == row),]$logFC
    }
    
    # Get information about missing values
    not_deg <- as.character(genes_order[which(genes_order$logFC == 0),]$geneNames)
    not_filt <- as.character(genes_order[- which(genes_order$geneNames %in% rownames(dataFilt)),]$geneNames)
    
    if(length(not_filt) > 0){
      not_deg <- not_deg[-which(not_deg %in% not_filt)]
      
      for(nofilt in not_filt){
        genes_order[which(genes_order$geneNames == nofilt),]$cancer <- paste0(cancer, "**")
      }
    }  
    if(length(not_deg) > 0){
      for(nodeg in not_deg){
        genes_order[which(genes_order$geneNames == nodeg),]$cancer <- paste0(cancer, "*")
        
      }
    }
    
    all_DEGs <- rbind(all_DEGs, genes_order)
    save(all_DEGs, file = "DEA/all_DEA_combined.rda")
  }
  
  # Add information about missing values
  for(gene in geneNames){
    nofilt_text <- ""
    nodeg_text <- ""
    genestar <- paste0(gene, "*")
    genestarstar <- paste(gene, "**")
    genes_var <- c(gene, genestar, genestarstar)
    genes_plot <- all_DEGs[which(all_DEGs$geneNames %in% genes_var),]

    for(genrow in nrow(genes_plot)){
      if(length(grep("*", genes_plot[genrow,]$cancer))>0){
        nodeg_text <- "* Not differentially expressed"
      }
      if(length(grep("*", genes_plot[genrow,]$cancer))>0){
        nofilt_text <- "** Not present after filtering"
      }
    }
    
    # Prepare y-axis proportions
    max_yaxis <- sort(genes_plot$logFC)[nrow(genes_plot)]
    max_yaxis <- max_yaxis + 2
    min_yaxis <- sort(genes_plot$logFC)[1]
    min_yaxis <- min_yaxis - 2
    
    # Make barplot
    pdf(paste0("plots/genes/", gene, "_DEA_barplot.pdf"))
    bp <- barplot(genes_plot$logFC,
                  names.arg = genes_plot$cancer,
                  ylim = c(min_yaxis,max_yaxis),
                  cex.axis=1.3,
                  cex.names=1,
                  main = paste(gene, "differential expression"),
                  ylab = "logFC",
                  col = genes_plot$colours,
                  las = 2,
                  axes = TRUE)
    mtext(nodeg_text, side = 3, adj = 0)
    mtext(nofilt_text, side = 3, adj = 1) 
    dev.off()
    
  }
  
}