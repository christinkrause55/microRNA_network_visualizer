## microRNA network visualizer (miRNA NVis)
The idea behind this tool was to identify interaction networks between microRNAs and mRNAs of their potential target genes beyond web-based database query. Therefore, first an adjacency matrix was generated based to interaction databases, which is followed by seed prediction based on a simple seed match. The generation of the adjacency matrix is based on MATLAB, along with high throughput correlation and linear regression analysis between the microRNA expression and (metabolic) traits of the underlying study population. Targeted network visualisation, heatmap and cluster analysis was generated with follow-up matlab scripts. Alternatively, the seed match based prediction and visualisation functions were also generated using R to make this tool available for a broader spectrum of scientists.  

### Implementation in R

These functions, containing a prediction and visualisation method, work for mouse and human datasets. An adjacency matrix is included for both species to link information from public databases (validated miRNA-mRNA-interactions) to the seed prediction hits. Own expression data should be used to visualise edge weights in the networks. After publication of the data, a comprehensive repository regarding obese human liver expression for microRNAs and targets in metabolic diseases will be available. 

Nevertheless, also without own data, potential links can be visualised by creating an empty data frame with genes and microRNAs of interest. An example file for common metabolic genes is supplied as an example script. 

![plot](https://github.com/christinkrause55/microRNA_network_visualizer/blob/main/predictionplot.png)

Red edges indicate predicted (as 8mer or 7mer) target-microRNA interactions, whereby blue edges emerge from database entries. Key and title can be added in a customised fashion after execution of the drawPlot function. 

### Implementation in MATLAB
This tool framework implemented in MATLAB (R2020a, TheMathWorks) was used for the high throughput analysis of microRNA microarray data sets using linear regression and correlation models. Besides evaluation, MATLAB was mainly used to generate adjacency matrices from database entries. Therefore, a comprehensive SQLite database was generated from publuc data repositories, which was used to query target genes from every mature microRNA. The resulting matrices for human and mouse can be found in the R directory. The connection and database which manages the underlying tables will be published upon publication to ensure reproducible data usage. Currently only the frame for mouse based data processing is uploaded to this repo. A similar framework for human data will be uploaded after publication of the respective manuscript. 

The getSeedPrediction.m script can be used as a standalone function. An implementation of the seed prediction tool which is independent of input sequences for the respective potential target gene UTR is also located in the R directory. 

### Dependencies
Python script for processing fasta transcript output from ensembl biomart needs numpy (v1.19.5) and pandas (v1.3.1) with python 3.9. Resulting files are needed to use the seedPredicition functions in R (complete) or MATLAB (targeted) since these files contain sequence information. This script is included in case of working with other species than mouse or human. 

The pipeline in R (4.1.0) needs further attached packages dplyr (1.0.10), igraph (1.3.5), stringr (1.4.1), stringi (1.7.6) and biomaRt (2.50.3) to the base packages stats, graphics, grDevices utils, datasets, methods and base.   
