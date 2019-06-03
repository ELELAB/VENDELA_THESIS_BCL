To perform a differential expression analysis on subtypes run:
a) TCGA_subtypes.r with its corresponding functions in b) TCGA_subtypes_functions.r to do:
	1) Finding every subtype for every cancer type, and sorting out samples that belong to those. We will not continue with subtypes containing less than 5 samples per sample type (tumor and normal)
	2) Do a summary on every sample type for every sub type. Files will be saved in a directory called 'summary', which is created while running the script
	2) Do a differerential expression analysis on the remaining subtypes. Results will be saved in a directory called 'DEA', which is created while running the script
	3) Produce barplots of the differentially expressed genes. Results will be saved as .pdf files in a directory called 'plots', which is created while running the script
 

