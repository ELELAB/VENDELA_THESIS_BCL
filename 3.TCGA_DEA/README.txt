Do a differential expression analysis by comparing gene expression in tumor samples to normal samples from the filtered data produced in 1.TCGA_pre-process.
Run:
a) TCGA_DEA.r toghether with its corresponding functions b) TCGA_DEA_functions to produce:
	1) Summary files on the two different sample types, tumour and normal. The files will be saved in a directory called 'summary', which is created when running the script
	2) Differential expression analysis results from both all genes and only the genes of interest. The results will be saved as .rda files in a directory called 'DEA*, which is created when running the script
	3) Bar plots of the differentially expressed genes. Both by cancer types and genes. All plots will be saved as .pdf files in a directory called 'plots', which is created when running the script
 	4) Volcano plots on the differentially expressed genes. Plots will be saved in a the same directory as the barplots ('plots')

