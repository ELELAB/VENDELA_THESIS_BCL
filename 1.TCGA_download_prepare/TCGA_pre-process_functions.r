library(SummarizedExperiment)
library(TCGAbiolinks)

#------------------------ Functions for pre-processing TCGA data ---------------------------------------

# Get clinical information
get_clinical <- function(project){
  cat("Downloading clinical information. Cancer type =", cancer, "\n")
  dir.create("clinical", showWarnings = FALSE)
  dataClin <- GDCquery_clinic(project = project,"clinical")
  save(dataClin, file = paste0("clinical/",cancer,"_clinical.rda", sep = ""))
  return(dataClin)
}

# Prepare samples
prepare_data <- function(gdc){
  cat("Preparing data. Cancer type =", cancer, "\n")
  dir.create("pre-processing", showWarnings = FALSE)
  dataPrep <- TCGAanalyze_Preprocessing(object = gdc,
                                        cor.cut = 0.6,
                                        filename = paste0("pre-processing/dataPrep_",cancer, "_array_array.png"))
  rownames(dataPrep) <- rowData(gdc)[match(rownames(dataPrep), rowData(gdc)[,"ensembl_gene_id"]),"external_gene_name"]
  save(dataPrep, file = paste0("pre-processing/",cancer,"_dataPrep.rda"))
  return(dataPrep)
} 

# Normalize samples
normalize_data <- function(dataPrep){
  cat("Normalizing data. Cancer type = ",cancer,"\n" )
  dir.create("normalization", showWarnings = FALSE)
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                        geneInfo = geneInfo,
                                        method = "gcContent")
  save(dataNorm,file=paste0("normalization/",cancer,"_dataNorm.rda", sep = ""))
  return(dataNorm)
}

# Filter samples
filter_data <- function(dataNorm){
  cat("Filter data. Cancer type =", cancer,"\n")
  dir.create("filtering", showWarnings = FALSE)
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                    method = "quantile",
                                    qnt.cut =  0.25)
  save(dataFilt, file=paste0("filtering/",cancer,"_dataFilt.rda", sep = ""))
  return(dataFilt)
}
