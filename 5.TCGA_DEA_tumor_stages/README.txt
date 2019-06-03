To perform a differential expression analysis on tumour stages run:
a) TCGA_clinical_stages.r. Functions will be called from 3.TCGA_DEA/TCGA_DEA_functions.r and 4.TCGA_DEA_subtypes/TCGA_subtypes_functions.r. Run TCGA_clinical_stages.r to:
	1) Find every tumour stage for every cancer type, and sort out samples that belong to those. We will not continue with stages containing less than 5 samples per sample type (tumor and normal)
	2) Do a summary on every sample type for every tumour stage. Files will be saved in a directory called 'summary', which is created while running the script
	2) Do a differerential expression analysis on the remaining stages. Results will be saved in a directory called 'DEA', which is created while running the script
	3) Produce barplots of the differentially expressed genes in every stage. Results will be saved as .pdf files in a directory called 'plots', which is created while running the script
 

