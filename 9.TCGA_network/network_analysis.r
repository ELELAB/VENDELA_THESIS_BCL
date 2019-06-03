library(igraph)
library(biomaRt)
library(dplyr)
library(SummarizedExperiment)
setwd("9.TCGA_network")
source("network_analysis_functions.r")


#read iid annotation file
annot_file <- read.delim(file = "human_annotated_PPIs.txt", sep = "\t")

#list of uniprot names for proteins
proteins_list <- c("P10415", #bcl-2
              "Q07817", #bcl-xl
              "Q92843", #bcl-w
              "Q16548", #bfl-1
              "Q07820", #mcl-1
              "Q9HD36", #bcl2l10
              "Q07812", #bax
              "Q16611", #bak1
              "Q9UMX3" #bok
)

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

# Remove duplicated interactions
iid_dup <- iid_remove_duplicates(annot_file, protein_list)

# Remove interactions with less than 2 pmids avaiable
iid_pmid <- iid_remove_with_few_pmid(iid_dup)

# Get DEG values for every interacting protein from the differential expression analysis file
iid_ens <-iid_get_ens_DEG(iid_pmid, cancer_types)

#iid_ens <- get(load("protein_interactions_DEA.rda"))
iid_multiple_ens <- iid_get_multiple_ens(iid_ens)

# From iid_multiple_ens, find which ensembl IDs to remove and remove them
ensembl_remove <- c("ENSG00000258643", "ENSG00000205571")
iid_single <- iid_get_single_ens(iid_ens, ensembl_remove) 

# Prepare file for networkplot
iid_network <- prepare_for_networkplot(iid_single, iid_pmid, cancer_types)

# Do a network plot
do_network_plot(iid_network)
