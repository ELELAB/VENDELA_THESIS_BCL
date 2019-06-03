Before doing a network analysis, you need to have an IID2 annotatino file in protein-protein interactions (see 'human_annotated_PPIs.txt')

Then, run a) network_analysis.r with its corresponding functions in b) network_analysis_functions.r to:
	1) Convert uniprot ids to ensembl ids
	2) get DEG values for every interacting protein. Results is saved in a directory called 'protein_interactions', which is created while running the script
	3) plot interactions as a network plot. Results will be saved as a pdf. file in a directory called 'plots', which is created while running the script
