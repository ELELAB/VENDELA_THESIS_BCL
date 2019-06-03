# Perform a PCA and return an prcomp object
get_pca <- function(dataFilt, filename_pca){
  dir.create("PCA", showWarnings = FALSE)
  data_t <- t(dataFilt)
  data_pca <- prcomp(data_t, scale. = TRUE, center = TRUE)
  save(data_pca, file = paste0("PCA/",filename_pca))
  return(data_pca)
}

# Get information on which sample are from tumor and which are from normal
get_sample_types <- function(dataFilt){
  obj_names <- colnames(dataFilt)
  tcgaIDs <- get_IDs(dataFilt)
  obj_groups <- as.factor(tcgaIDs$condition)
  obj_groups <- as.data.frame(obj_groups)
  sample_types <- cbind(obj_names, obj_groups)
  return(sample_types)
}

# Get PC
plot_variations <- function(dataFilt_pca, sample_types, summary_file, variance_plotfile, cumulative_plotfile, pca_plotfile){
  dir.create("summary", showWarnings = FALSE)
  dir.create("plots", showWarnings = FALSE)
  
  # Do a summary on the prcomp file
  dataFilt_pcasum <- summary(dataFilt_pca)
  sum_importance <- dataFilt_pcasum$importance
  save(sum_importance, file = paste0("summary/",summary_file))
  
  # Do a cumulative summary on variances
  summary_percent_cumsum <- sum_importance[3, ] * 100
  summary_truncate_cumsum <- summary_percent_cumsum[1:64, drop = FALSE]
  summary_truncate_cumsum<- as.data.frame(summary_truncate_cumsum)
  colnames(summary_truncate_cumsum) <- c("Var")
  
  # Plot cumuliative summary
  ggplot(summary_truncate_cumsum, aes(x = rownames(summary_truncate_cumsum) , y = Var)) +
    coord_cartesian(xlim = c(0, 10))+
    geom_bar(stat = "identity", color = "darkorchid1", fill = "purple") + 
    scale_x_discrete(limits=unique(as.character(rownames(summary_truncate_cumsum))), breaks = (rownames(summary_truncate_cumsum)[c(1,50,100,150,200,250,300,350,400,450,500,550,600)])) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +xlab("PCs") +
    ylab("% of Cumulative Variance Explained by PC components") +
    ggtitle("Cumulative Variance explained by different PCs") + theme(plot.title = element_text(hjust = 0.5)) +
    ggsave(paste0("plots/", cumulative_plotfile))
  
  # Get variances
  prop_of_variance <- sum_importance[2,]
  variance_percent <- prop_of_variance * 10
  variance_percent <- as.data.frame(variance_percent)
  variance_truncate <- variance_percent[1:10, , drop = FALSE]
  
  # Do a PC1-PC2 plot
  PC1_var <- variance_percent["PC1",]
  PC2_var <- variance_percent["PC2",]
  xlab <- paste0("PC1"," (",PC1_var,"%)")
  ylab <- paste0("PC2"," (",PC2_var,"%)")
  
  ggplot(as.data.frame(dataFilt_pca$x), aes(x=PC1, y=PC2, color = factor(sample_types$obj_groups))) +
    geom_point(size = 3, alpha = 0.3) + xlab(xlab) + ylab(ylab)  + 
    labs(title = paste0(cancer," PCA"), color = "Groups\n") +
    scale_color_manual(labels = c("Tumor", "Normal"),values = c("blue", "red")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggsave(paste0("plots/",pca_plotfile))
  
  # Plot variance by PC
  ggplot(variance_truncate, aes(x = rownames(variance_truncate), y = variance_truncate$variance_percent)) +
    geom_bar(stat = "identity", color = "darkorchid1", fill = "purple") + scale_x_discrete(limits=unique(as.character(rownames(variance_truncate)))) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +xlab("PCs") + ylab("% of Variance Explained by PC components") +
    ggtitle("A View of % of Variance explained by different PCs") + theme(plot.title = element_text(hjust = 0.5)) +
    ggsave(paste0("plots/", variance_plotfile))

}




