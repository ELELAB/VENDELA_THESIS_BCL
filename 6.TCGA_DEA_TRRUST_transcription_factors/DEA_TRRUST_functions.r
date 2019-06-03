make_tf_df <- function(dataDEGs_tf_genes, tflist){
  for(gene in rownames(dataDEGs_tf_genes)){
    tf_gene <- c()
    genpos <- which(rownames(dataDEGs_tf_genes) == gene)
    for(tf in 1:length(tflist)){
      for(genpos1 in 1:length(tflist[[tf]])){
        if(strsplit(tflist[[tf]][genpos1], "_")[[1]][1] == gene){
          gene1 <- paste0(".", tf)
          gene11 <- ""
          if(strsplit(tflist[[tf]][genpos1], "_")[[1]][2] == "Activation"){
            gene11 <- paste0(gene1, "+")
          }else{
            gene11 <- paste0(gene1, "-")
          }
          
        }
      }
      tf_gene <- paste0(gene, gene11)
      rownames(dataDEGs_tf_genes)[genpos] <- tf_gene
    }
  }
  
  return(dataDEGs_tf_genes)
}

barplot_tf <- function(dataDEGs_tf_genes, filename_barplot, geneNames, cols){
  
  dir.create("plots", showWarnings = FALSE)
  # Do barplot
  
  
  pdf(paste0("plots/", filename_barplot), onefile = TRUE)
  par(mar=c(8.5, 4.1, 4.1, 2.1))
  barplot(dataDEGs_tf_genes$logFC,
              names.arg = rownames(dataDEGs_tf_genes),
              ylim = c(-4,6),
              cex.axis=1,
              cex.names=1,
              main = paste(cancer, " transcription factor differential expression"),
              ylab = "logFC",
              col = cols,
              las = 2,
              axes = TRUE,
              beside = TRUE)
  mtext(paste0("1 = ",geneNames[1], ", 2 = ", geneNames[2], ", 3 = ", geneNames[3], ", 4 = ", geneNames[4], ", 5 = ", geneNames[5], "\n +/- = Activator/Repressor"),
        side = 1, line = 6.5)
  dev.off()
}