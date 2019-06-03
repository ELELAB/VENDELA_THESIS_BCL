Vendela Rissler
Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark

This repository contains curated RNASEQ data from tumor and normal solid tissue samples obtained from TCGA, and the recount2 initiative along with their analyses. The repository was made with the intent of openly sharing both the raw input data used at the time of the analyses and the R-scripts employed to carry out the study.

We here present a TCGA Pan-Cancer analysis. Data from 18 different cancer types from The Cancer Genome Atlas (TCGA) have been used to analyse genetic variations in the BCL-2 family.


The repository contains the following folders which need to be used sequentially. Each folders contain a different README file with further information:

(1) 1.TCGA_download_prepare: Downloading, aggregating and preparing TCGA data
(2) 2.PCA: A principal component analysis
(3) 3.TCGA_DEA: A differential expression analysis
(4) 4.TCGA_DEA_subtypes: A differential expression analysis on cancer subtypes
(5) 5.TCGA_DEA_tumor_stages: A differential expression analysis on tumor stages
(6) 6.TCGA_DEA_TRRUST_transcription_factors: Differential expression analysis results on transcription factors that regulate BCL-2 gene expression. The transcription factors were selected from TRRUST data base
(7) 7.CEMiTOOL_coexpression: A co-expression analysis using CEMiTool
(8) 8.TCGA_somatic_variations: Finding mutations (SNPs) and their impact
(9) 9.TCGA_network: A network analysis using information from the IID2 annotation file
(10) 10.Survival_analysis: A Kaplan-Meier survival analysis, analysing gene expression in correlation to survival


Used packages:

TCGAbiolinks version 2.11.3 
CEMiTool version 1.0.3
SummarizedExperiment version 1.8.1
ggplot2 version 3.1.0
stats version 3.4.4
dplyr version 0.7.8
biomaRt version 2.34.2
igraph version 1.1.2
GenomicDataCommons version 1.2.0
limma version 3.34.9
survival version 2.43.3
survminer version 0.4.3
plyr version 1.8.4
ggfortify version 0.4.6


We suggest to use Rstudio version 3.4.4 to run the scripts of interest so that you can follow the analyses one line at the time and digest the results.

