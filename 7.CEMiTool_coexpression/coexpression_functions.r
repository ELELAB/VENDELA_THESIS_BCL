

get_samples <- function(project, tumor_sample){
  cat(paste("Downloading query", "\n"))
  query.exp <- GDCquery(project = project,
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - Counts",
                        sample.type = c(tumor_sample),
                        legacy = FALSE)
  GDCdownload(query.exp)
  return(query.exp)
}

# Download .rda files from GDC
get_GDC <- function(query.exp){
  cat("Downloading .rda files." , "\n")
  tumor_rd <- GDCprepare(query = query.exp)
  return(tumor_rd)
}
