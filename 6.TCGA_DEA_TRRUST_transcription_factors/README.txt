To perform this analysis, you have to:
(1) go to https://www.grnpedia.org/trrust/ and collect the transcription factors from the TRRUST data base that regulate expression of the genes of interest. Keep only those with 2 or more pmids, and a known regulatory effect.
(2) Make a vector with tha name of the transcription factors and their regulatory effect, for example:
	tf_effect <- c("tf1_Repressor", "tf2_Activator") (see example in DEA_TRRUST_tf.r)
(3) Run DEA_TRRUST_tf.r with its corresponding functions in DEA_TRRUST_functions.r to:
	a)Sort out the transcription factors from the file from the differential expression analysis, and add information about their regulatory effect. The results will be saved in a directory called 'DEA', which is created while running the script
	b) Plot DEA results as barplots. They will be saved as .pdf files in a directory called 'plots', which is created while running the script 
