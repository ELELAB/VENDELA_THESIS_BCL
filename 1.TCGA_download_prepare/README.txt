Following scripts are provided for downloading, aggregating, pre-processing and plotting 

(1) a) TCGA_download.r, together with corresponding functions in b) TCGA_download_functions.r
Here TCGA data from GDC will be first downloaded, and then aggregared to a .rda file. A all data and .rda files will be saved in a directory called 'data', which is created when running the script

(2) a) TCGA_pre-process.r together with corresponding functions in b) TCGA_pre-process_functions.r
Here the aggregated .rda file from 'TCGA_download.r' (1.a) will be pre-processed in three steps:
	1) pre-processing, where samples with a pearson correlation lower than 0.6 will be removed. The results will be saved as .rda files in a drectory called 'pre-processing', which is created when running the script
	2) normalization, where pre-processed file will be normalized to GC-count and sequencing depth. The results will be saved as .rda files in a directory called 'normalization', which is created when running the script
	3) filtering, where samples after the normalization step with low expression across samples will be removed. Results will be saved as .rda files in a directory called 'filtering', which is created when running the script
Clinical information will also be collected and saved as .rda files in a directory called 'clinical'

(3) TCGA_boxplot.r, to produce boxlpots on the filtered files. The boxplots will picture gene expression. Output will be .pdf files saved in a directory called 'plots', which is created when running the script  
