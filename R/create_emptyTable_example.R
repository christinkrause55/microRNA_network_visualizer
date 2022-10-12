### Generate empty matrix as input

# Function to get genes based on pathway  identifier
getGenelist_mouse <- function(go_ID){
  library(biomaRt)
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
  gene.data <- getBM(attributes=c('external_gene_name'), filters = 'go', values =  go_ID, mart = ensembl)
  return(gene.data)
}

# Choice of different metabolic pathways
pathway_list <- list()
pathway_list <- append(pathway_list, list("Insulin receptor signaling pathway"="GO:0008286"))
pathway_list <- append(pathway_list, list("Gluconeogenesis"="GO:0006094"))
pathway_list <- append(pathway_list, list("Lipid biosynthetic process"="GO:0008610"))
pathway_list <- append(pathway_list, list("Canonical glycolysis"="GO:0061621"))
pathway_list <- append(pathway_list, list("Sterol biosynthetic process"="GO:0016126"))
pathway_list <- append(pathway_list, list("lipid metabolism"="GO:0006629"))
pathway_list <- append(pathway_list, list("Glycogen synthesis"="GO:0005978"))
pathway_list <- append(pathway_list, list("Glycogen breakdown"="GO:0005980"))
pathway_list <- append(pathway_list, list("Insulin response"="GO:0032868"))
pathway_list <- append(pathway_list, list("Lipolysis"="GO:0016042"))
pathway_list <- append(pathway_list, list("Lipogenesis"="GO:0008610"))
pathway_list <- append(pathway_list, list("Cholesterol metabolism"="GO:0008203"))

pt <- 12
genes <- getGenelist_mouse(pathway_list[[pt]])
microRNA_list <- c("mmu-let-7e-5p","mmu-miR-16-5p","mmu-miR-34a-5p","mmu-miR-149-5p","mmu-miR-182-5p")

ids <- c(genes$external_gene_name, microRNA_list)

# create empty frame
empty_table <- data.frame(matrix(1,4,length(ids)))
colnames(empty_table) <- ids

# Call functions
results <- miRNA_targetome_prediction(species = "mmu", 
                                      addToAdjacency = TRUE, geneNames = genes$external_gene_name, 
                                      complete_mR = TRUE, microRNANames = microRNA_list, complete_miR = TRUE)

drawPlot(microRNA_list, empty_table, results$adjacency_matrix, strict=FALSE, layout = "circle", emptymode = TRUE)
