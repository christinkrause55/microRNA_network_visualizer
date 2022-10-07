## microRNA network visualizer (miRNA NVis)
This tool framework implemented in MATLAB (R2020a, TheMathWorks) can be used for high throughput analysis of microRNA microarray data sets.
A database connection and database which manages the underlying data will be published upon publication to ensure reproducible data usage. 

Currently only mouse based data processing is uploaded to this repo. A similar framework for human data will be uploaded after publication of the respective manuscript. 

The getSeedPrediction.m script can be used as a standalone function. An implementation of the seed prediction tool which is independent of input sequences for the respective potential target gene UTR is also located in the R directory. This script works for mouse and human datasets. An adjacency matrix is included for both species to link information from public databases (validated miRNA-mRNA-interactions) to the seed prediction hits.  
