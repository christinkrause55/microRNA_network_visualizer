#' Plot networks between microRNA and target genes based on created adjacency matrix
#' 
#' @description 
#' The drawPlot function uses the information from the adjacency matrix to create network plots.
#' Moreover, it calculates correlation estimates and p values from RNAseq or qPCR data to create
#' weights for edges.
#' 
#' @param microRNA_list A list of microRNAs matching to the used adjacency matrix to use as central regulators
#' @param counts A complete dataframe with microRNA and target gene values. Column names should contain microRNA 
#' or gene names matching to the adjacency matrix without the numeric identifier. If no correlation should be performed,
#' just fill the dataframe with ones. Data can be prepared by the preparation function.
#' @param adjacency An adjacency matrix, which could be created by the miRNA_targetome_predicition function. Row should
#' contain microRNA names, columns should contain target gene names.
#' @param strict Boolean to apply a strict filter of negative estimate and significant p value
#' @param threshold Adjusted p value filter threshold for a strict plotting
#' @param layout Indicates, which layout the plot should have. Use ["default","circle","tree"] 
#' @param emptymode Boolean to indicate whether real count data was used
#' 
#' @return A plot with the regulatory network
#' 
#' @author Christin Krause

drawPlot <-
  function(microRNA_list = c("hsa-miR-182-5p"),
           counts,
           adjacency, strict=FALSE, threshold = 0.05,
           layout = "default", emptymode = FALSE) {
    
    suppressMessages(library(dplyr))
    suppressMessages(library(igraph))
    
    if(strict & emptymode){
      message("This won't work, strict is reversed")
      strict <- FALSE
    }
    
    # Setup correlation dataframes
    final_results_complete <- data.frame()
    for (k in c("X1", "padj", "estimate")) {
      final_results_complete[[k]] <- as.numeric()
    }
    for (k in c("target", "regulator", "gene")) {
      final_results_complete[[k]] <- as.character()
    }
    
    RNAseq_correlate_p <- data.frame()
    RNAseq_correlate_r <- data.frame()
    
    for (k in colnames(counts)) {
      RNAseq_correlate_p[[k]] <- as.numeric()
      RNAseq_correlate_r[[k]] <- as.numeric()
    }
    
    # Fill dataframes
    for (microRNA_name in microRNA_list) {
      message(microRNA_name)
      
      # Perform correlation analysis
      miR_position <- which(colnames(counts) == microRNA_name)
      for (k in 1:ncol(counts)) {
        if(!emptymode){
          test <-
            cor.test(as.numeric(counts[, k]),
                     counts[, miR_position])
          RNAseq_correlate_p[1, k] <- test$p.value
          RNAseq_correlate_r[1, k] <- test$estimate
        }else{
          RNAseq_correlate_p[1, k] <- 1
          RNAseq_correlate_r[1, k] <- -0.01
        }
      }
      
      RNAseq_correlate_pt <- data.frame(t(RNAseq_correlate_p))
      RNAseq_correlate_rt <- data.frame(t(RNAseq_correlate_r))
      
      RNAseq_correlate_pt$padj <-
        p.adjust(RNAseq_correlate_pt$X1, method = 'BH')
      RNAseq_correlate_pt$estimate <- RNAseq_correlate_rt$X1
      
      # Filter and extract relevant regulator - target interactions
      if(strict){
      relevant <-
        RNAseq_correlate_pt[RNAseq_correlate_pt$padj < threshold &
                             RNAseq_correlate_pt$estimate < 0 , ]
      }else{
        relevant <- RNAseq_correlate_pt
      }
      
      c_adj <- strsplit(colnames(adjacency), split = "_")
      new_c_adj <- character()
      for (i in 1:length(c_adj)) {
        new_c_adj <- append(new_c_adj, c_adj[[i]][1])
      }
      
      for (gene in rownames(relevant)) {
        if (gene %in% new_c_adj) {
          relevant$target[rownames(relevant) == gene] <-
            adjacency[which(rownames(adjacency) == microRNA_name), which(new_c_adj == gene)][1]

          relevant$regulator[rownames(relevant) == gene] <-
            microRNA_name
        }
      }
      
      final_results <- na.omit(relevant[relevant$target == 1 | relevant$target == 2,])
      final_results$gene <- rownames(final_results)
      rownames(final_results) <- NULL
      final_results_complete <-
        rbind(final_results_complete, final_results)
    }
    
    # Plot network
    g <- character()
    for (i in 1:nrow(final_results_complete)) {
      g <- append(g, final_results_complete$regulator[i])
      g <- append(g, final_results_complete$gene[i])
    }
    
    if(!is.na(g[1])){
      net <- graph(edges = g, directed = FALSE)
      
      v_names <- rep("target", length(V(net)$name))
      v_size <- rep(10, length(V(net)$name))
      V_label_size <- rep(1.25, length(V(net)$name))
      for (microRNA_name in microRNA_list) {
        v_names[which(V(net)$name == microRNA_name)] <- 'miRNA'
        v_size[which(V(net)$name == microRNA_name)] <- 7.5
        V_label_size[which(V(net)$name == microRNA_name)] <- 1.5
      }
      V(net)$type <- v_names
      
      e_size <- numeric()
      e_color <- numeric()
      for (i in 1:nrow(final_results_complete)) {
        if (final_results_complete$estimate[i] > 0) {
          e_size[i] <- 1
          e_color[i] <- 1
        }else{
          if(!emptymode){
            if(final_results_complete$estimate[i] > 0){
              e_size[i] <- (final_results_complete$estimate[i] * 10) + 1
            }else{
              e_size[i] <- (final_results_complete$estimate[i] * -10) + 1
            }
            
          }else{
            e_size[i] <- 1
          }
          if(final_results_complete$target[i] == 1){
            e_color[i] <- 2
          }else{
            e_color[i] <- 3
          }
        }
      }
      e_size[e_size>5] <- 5
      
      if(layout == "default"){
        plot(
          net,
          vertex.size = v_size,
          vertex.color = c("goldenrod2", "moccasin")[1 + as.numeric(!V(net)$type ==
                                                                        "miRNA")],
          vertex.label.color = "black",
          vertex.label.degree = 45,
          vertex.frame.color = "azure",
          vertex.label.dist = 0.5,
          vertex.label.cex = V_label_size,
          edge.color = c("darkgrey","steelblue1","red1")[e_color],
          edge.width = e_size
        )
      }else if(layout == "circle"){
        plot(
          net,
          vertex.size = v_size,
          vertex.color = c("coral", "lightskyblue")[1 + as.numeric(!V(net)$type ==
                                                                        "miRNA")],
          layout=layout_in_circle,
          vertex.label.color = "black",
          vertex.label.degree = 45,
          vertex.frame.color = "azure",
          vertex.label.dist = 1,
          vertex.label.cex = V_label_size,
          edge.color = c("darkgrey","steelblue1","red1")[e_color],
          edge.width = e_size
        )
      } else{
        plot(
          net,
          vertex.size = v_size,
          vertex.color = c("coral", "lightskyblue")[1 + as.numeric(!V(net)$type ==
                                                                     "miRNA")],
          layout=layout_as_tree,
          vertex.label.color = "black",
          vertex.label.degree = 45,
          vertex.frame.color = "azure",
          vertex.label.dist = 0.5,
          vertex.label.cex = V_label_size,
          edge.color = c("darkgrey","steelblue1","red1")[e_color],
          edge.width = e_size
        )
      }
    }else{
      plot(make_ring(1))
    }
    return(final_results_complete)
  }
