# gene-gini-matlab
Codes for "What are housekeeping genes?"
Chintan J. Joshi, Eyleen J. O’Rourke, Nathan E. Lewis --biorXiv link--

# NOTES


# Loading transcriptomics data (rnaData)
The transcriptomics data is required as a MATLAB structure containing following fields:
a. gene: list of gene names
b. value: a matrix whose each column contains the expression of genes across a tissue/context and each row contains expression of a gene across all tissues/contexts
c. genesymbols (optional): any alternative gene names that the user may want to keep track of
d. Tissue: names of conditions
