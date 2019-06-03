setwd("./6.TCGA_DEA_TRRUST_transcription_factors")
source("DEA_TRRUST_functions.r")


# This information has been collected from TRRUST database
# transcription factors with less than 2 pmid or unknown regulatory  effect have been removed 
# Add name on transcription factor and 
#"_Activation" if it is activating BCL2 gene expression
#"_Repression" if it is repressing BCL2 gene expression
bcl2_regulators <- c("E2F1_Activation", "GLI1_Activation", "GLI2 _Activation" , "HOXA1_Activation" , "ING4_Repression",  "NFKB1_Both", "PAWR_Repression", "RELA_Activation", "RELA_Repression", "STAT3_Activation" , "TP53_Repression", "ZNRD1_Activation")

bcl2l1_regulators <- c("NFKB1_Activation", "RELA_Activation", "STAT3_Activation")

# No information about BCL2L2 regulators in TRRUST data base
#bcl2l2_regulators <- c()

# No information about BCL2L10 in TRRUST data base
#bcl2l10_regulators <- c()

bcl2a1_regulators <- c("NFKB1_Activation", "RELA_Activation")

mcl1_regulators <- c("STAT3_ Activation")

bax_regulators <- c("ING1_Activation","TP53_Activation")

# No information about BAK1
#bak1_regulators <- c()

# No information about BOK
#bok_regulators <- c()

# Get all transcription factor names
all_tf_func <- c(bcl2_regulators, bcl2l1_regulators, bcl2a1_regulators, mcl1_regulators, bax_regulators)
all_tf <- c()
for(tf in 1:length(all_tf_func)){
  all_tf <- c(all_tf, strsplit(all_tf_func[tf], "_")[[1]][1])
  
}
all_tf <- unique(all_tf)

# List with the different regulators. Remember their positions (the barlpot names will use those positions!)
tflist <- list(bcl2_regulators, bcl2l1_regulators, bcl2a1_regulators, mcl1_regulators, bax_regulators)

# List of HUGO gene names for which transcription factor information was found. Make this same order as 'tflist'
geneNames <- c("BCL2", "BCL2L1", "BCL2A1", "MCL1", "BAX")

#Colours to use in barplot. Must have same length as 'all_tf'
colours_tf <- c("#FF0000FF","#FF8000FF","#FFFF00FF","#80FF00FF","#00FF00FF","#00FF80FF","#00FFFFFF","#0080FFFF","#0000FFFF","#8000FFFF","#FF00FFFF","#FF0080FF")


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


for(cancer in cancer_types){
  # Get dataDEGs for the transcription factors
  dir.create("DEA", showWarnings = FALSE)
  dataDEGs <- get(load(paste0("../3.TCGA_DEA/DEA/",cancer, "_dataDEGs.rda")))
  dataDEGs_tf_genes <- dataDEGs[which(rownames(dataDEGs) %in% all_tf),]
  save(dataDEGs_tf_genes, file = paste0("DEA/",cancer, "_transcription_factors_DEA.rda"))
  dataDEGs_tf_genes <- as.data.frame(dataDEGs_tf_genes)
  tf_df <- make_tf_df(dataDEGs_tf_genes, tflist)
  cols <- colours_tf[which(all_tf %in% rownames(dataDEGs_tf_genes))]
  # Do barplot
  #filename_barplot = name of pdf barplot file 
  filename_barplot <- paste0(cancer, "_TRRUST_tf_barplot.pdf")
  barplot_tf(tf_df, filename_barplot, geneNames, cols)

}  
 
