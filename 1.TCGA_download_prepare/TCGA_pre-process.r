setwd("1.TCGA_download_prepare")

# Get functions for pre-processing TCGA data
source("TCGA_pre-process_functions.r")

# List of cancer types to pre-process TCGA data. There must be downloaded and aggregated GDC files for every cancer_type
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
                  "CHOL")


#------------------------ Pre-process TCGA data for every cancer type in cancer_types---------------------------------------------

for(cancer in cancer_types){ 

  # Get TCGA project that corresponds to the cancer type
  project <- paste0("TCGA-", cancer)

  # get the GDC .rda file from TCGA_download
  gdc <- get(load(paste0("data/",tolower(cancer), ".exp.rda")))

  # Get clinical information
  dataClin <- get_clinical(project)
  
  # Pre-process gdc
  dataPrep <- prepare_data(gdc)
  
  # Normalize the pre-processed data with gc-count
  dataNorm <- normalize_data(dataPrep)
  
  # Filter normalized samples, with mean filtering threshold is 0.25
  dataFilt <- filter_data(dataNorm) 
}
