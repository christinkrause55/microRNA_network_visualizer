## microRNA network visualizer (miRNA NVis)
The idea behind this tool was to identify interaction networks between microRNAs and mRNAs of their potential target genes beyond web-based database query. Therefore, first an adjacency matrix was generated based to interaction databases, which is followed by seed prediction based on a simple seed match. The generation of the adjacency matrix is based on MATLAB, along with high throughput correlation and linear regression analysis between the microRNA expression and (metabolic) traits of the underlying study population. Targeted network visualisation, heatmap and cluster analysis was generated with follow-up matlab scripts. Alternatively, the seed match based prediction and visualisation functions were also generated using R to make this tool available for a borader spectrum of scientists.  

### Implementation in R

These functions, containing a predicition and visualisation method, work for mouse and human datasets. An adjacency matrix is included for both species to link information from public databases (validated miRNA-mRNA-interactions) to the seed prediction hits. Own expression data should be used to visualise edge weights in the networks. After publication of the data, a comprehensive repository regarding obese human liver expression in metabolic diseases will be available.

### Implementation in MATLAB
This tool framework implemented in MATLAB (R2020a, TheMathWorks) was used for the high throughput analysis of microRNA microarray data sets.
A database connection and database which manages the underlying data will be published upon publication to ensure reproducible data usage. Currently only mouse based data processing is uploaded to this repo. A similar framework for human data will be uploaded after publication of the respective manuscript. 

The getSeedPrediction.m script can be used as a standalone function. An implementation of the seed prediction tool which is independent of input sequences for the respective potential target gene UTR is also located in the R directory. 
