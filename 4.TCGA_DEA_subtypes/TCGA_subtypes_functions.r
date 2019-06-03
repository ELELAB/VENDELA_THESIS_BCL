#----------------------------------- Functions to perform DEA on cancer subtypes -------------------------------------

get_full_barcodes <- function(cancer, barcodes, tumor_sample, normal_sample){
  
  setwd("../1.TCGA_download_prepare/data")
  project <- paste0("TCGA-", cancer)
  stage.query.exp <- GDCquery(project = project,
                              data.category = "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification",
                              workflow.type = "HTSeq - Counts",
                              barcode = barcodes,
                              sample.type = c(tumor_sample, normal_sample),
                              legacy = FALSE)
  
  #GDCdownload(query.exp)
  setwd("../../4.TCGA_DEA_subtypes")
  barcodes_long <- getResults(stage.query.exp,cols="cases")
  return(barcodes_long)
}



DEA_subtype_barplot <- function(dataFilt_sub, subtype, dataDEGs_sub_genes, colours_genes, geneNames, filename_barplot){
  dir.create("plots", showWarnings = FALSE)
  
  # Make dataframe for barplot
  genes_order <- as.data.frame(geneNames)
  genes_order$geneNames <- as.character(genes_order$geneNames) 
  genes_order$logFC <- rep(0,nrow(genes_order))
  genes_order$colours <- as.character(colours_genes)
  for(row in rownames(dataDEGs_sub_genes)){
    genes_order[which(genes_order$geneNames == row),]$logFC <- dataDEGs_sub_genes[which(rownames(dataDEGs_sub_genes) == row),]$logFC
  }
  
  # Add information about missing values
  not_deg <- as.character(genes_order[which(genes_order$logFC == 0),]$geneNames)
  not_filt <- as.character(genes_order[- which(genes_order$geneNames %in% rownames(dataFilt_sub)),]$geneNames)
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
          cex.names=0.6,
          main = paste(cancer, subtype, "differential expression"),
          ylab = "logFC",
          col = genes_order$colours,
          axes = TRUE)
  mtext(nodeg_text, side = 3, adj = 0)
  mtext(nofilt_text, side = 3, adj = 1) 
  dev.off()
  return(bp)
}

