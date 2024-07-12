This directory shows the minimal set to use the prediction tool in MATLAB. The workspace variable contains examples for miRNA and mRNA input variables, microRNA expression values and IDs. 

Start this tool by loading the workspace and opening of the script 'drawGraph.m'. The variable 'traitTable' needs to be replaced by microRNA expression values. Moreover, a SQLite database 'GeneExpression.db' 
with target gene expression tables (form of 'GeneName_expression') and the columns 'SampleName' and 'dCT' is needed. Edit microRNA and mRNA names and execute the script.
