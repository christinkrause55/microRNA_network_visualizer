#' Seed prediction of miRNAs for the regulation of target genes by base pairing with the 3' UTR sequence
#' 
#' @description 
#' seedPredictionTargetome returns a table with all information for the miRNA targetome
#' 
#' @details 
#' This function predicts miRNA binding mode (8mer, 7mer or 6mer) and analysis parameters of pairing 
#' conditions, such as AU content around the predicted seed, proximity to the canonical Pumilio binding motif,
#' relative position of the seed in the 3'UTR and additional backbone binding.
#' 
#' @param miRNA_ID name of a specific miRNA, such as hsa-miR-182-5p or mmu-miR-182-5p
#' @param mRNA_ID name of a specific mRNA, such as "LRP6" in human or "Lrp6" in mouse
#' @param species which species to analyse [hsa, mmu]
#' @param ownAdj boolean whether tool should use own matrix instead of predefined one
#' @param adjacency path to or dataframe of an adjacency matrix NxM whitch N rows as miRNA and M columns for mRNA. 
#' @param geneNames matrix with gene names of interest for a reduced set of mRNAs; set 'complete_mR' as TRUE for usage
#' @param complete_miR boolean if complete miRNA transcriptome should be checked
#' @param complete_mRNA boolean of complete mRNA transcriptome should be checked as potential targetome
#' @param nameOfTable name of the output table if 'writeTo is set TRUE
#' @param writeTo boolean whether table should be exported as CSV file
#' @param addToAdjacency boolean which indicates whether the prediction results should be added to an adjacency matrix
#' 
#' @return list with targetome and adjacency matrix
#' 
#' @author Christin Krause
#'
seedPredictionTargetome <- function(miRNA_ID = "hsa-miR-182-5p" ,
                                    mRNA_ID = "LRP6" ,
                                    nameOfTable = "output.csv",
                                    writeTo = FALSE,
                                    species = "hsa",
                                    adjacency = readRDS("miRNA_targetome_adjacency.rds"),
                                    ownAdj = FALSE,
                                    geneNames = "",
                                    complete_miR = FALSE,
                                    complete_mR = FALSE,
                                    addToAdjacency = FALSE)
{
  ### imports
  suppressMessages(library(stringi))
  suppressMessages(library(stringr))
  suppressMessages(library(dplyr))
  
  if(species == "hsa"){
    miRNA = readRDS(file = "human_miRNA.rds")
    mRNA = readRDS(file = "human_mRNA.rds")
    if(!ownAdj){
      adjacency = readRDS("miRNA_targetome_adjacency.rds")
    }
  }else if (species == "mmu"){
    miRNA = readRDS(file = "murine_miRNA.rds")
    mRNA = readRDS(file = "murine_mRNA.rds")
    if(!ownAdj){
      adjacency = readRDS("miRNA_targetome_adjacency_mouse.rds")
    }
  }else{
    message("Species not supported")
  }
  
  if(length(geneNames) > 1){
    mRNA_reduced <- mRNA[0,]
    for (gene in geneNames){
      t <- which(mRNA$Gene == gene)
      mRNA_reduced <- rbind(mRNA_reduced, mRNA[t,])
    }
    mRNA <- mRNA_reduced
  }
  
  ### Basic function for sequence processing
  # Re-Write Seed as reverse-complement sequence
  process_seed <- function(miR_transcript) {
    seed <- substring(miR_transcript, 2, 8)
    seed <- gsub('U', 'T', seed)
    # Reverse
    seed <- stri_reverse(seed)
    # Complement
    single <- unlist(strsplit(seed, ''))
    new_seed <- as.character()
    for (i in 1:length(single)) {
      if (single[i] == 'A')
        new_seed <- paste(new_seed, 'T')
      if (single[i] == 'C')
        new_seed <- paste(new_seed, 'G')
      if (single[i] == 'T')
        new_seed <- paste(new_seed, 'A')
      if (single[i] == 'G')
        new_seed <- paste(new_seed, 'C')
    }
    return(gsub(' ', '', new_seed))
  }
  
  # Get backbone pairing
  upstream_seed <- function(miR_transcript) {
    seed <- substring(miR_transcript, 11, 17)
    return(process_seed(seed))
  }
  
  # Consider longest transcript UTR
  process_mRNA <- function(geneName, mRNA) {
    seq <- tryCatch({
      mRNA <- mRNA[mRNA$Gene == geneName,]
      s <- 0
      seq <- ""
      for (item in mRNA$Sequence) {
        if (stri_length(item) > s) {
          s <- stri_length(item)
          seq <- item
        } else
          s <- s
        seq <- seq
      }
    },
    error = function(cond) {
      message("Could not find gene")
      return(NA)
    },
    finally = {
      message(paste("Longes transcript for", geneName, "is", nchar(seq), sep = " "))
      return(seq)
    })
    return(seq)
  }
  
  ### Functions for seed match prediction
  ### returns positions: Fist half represents start, second half respective end positions.
  contains_8mer <- function(miRNA_transcript, mRNA_transcript) {
    positions <-
      unlist(str_locate_all(mRNA_transcript, paste(miRNA_transcript, 'A', sep = "")))
    return(positions)
    
  }
  
  contains_7mer_A1 <- function(miRNA_transcript, mRNA_transcript) {
    positions <-
      unlist(str_locate_all(mRNA_transcript, paste(
        substring(miRNA_transcript, 2, 8), 'A', sep = ""
      )))
    return(positions)
  }
  
  contains_7mer_m8 <- function(miRNA_transcript, mRNA_transcript) {
    positions <-
      unlist(str_locate_all(mRNA_transcript, miRNA_transcript))
    return(positions)
  }
  
  contains_6mer <- function(miRNA_transcript, mRNA_transcript) {
    positions <-
      unlist(str_locate_all(mRNA_transcript, substring(miRNA_transcript, 2, 8)))
    return(positions)
  }
  
  contains_offset_6mer <-
    function(miRNA_transcript, mRNA_transcript) {
      positions <-
        unlist(str_locate_all(mRNA_transcript, substring(miRNA_transcript, 1, 6)))
      return(positions)
    }
  
  ### Functions for seed position quality estimation
  contains_AU <- function(pos, seq) {
    limit_seed_dw <- pos - 30
    limit_seed_up <- pos + 7 + 30
    if (limit_seed_up > nchar(seq)) {
      limit_seed_up <- nchar(seq)
    }
    
    if (limit_seed_dw >= 1 & limit_seed_up <= nchar(seq)) {
      roi <- substring(seq, limit_seed_dw, limit_seed_up)
    } else if (limit_seed_dw <= 0 & limit_seed_up <= nchar(seq)) {
      roi <- substring(seq, 1, limit_seed_up)
    } else if (limit_seed_dw >= 1 & limit_seed_up > nchar(seq)) {
      roi <- substring(seq, limit_seed_dw, nchar(seq) + limit_seed_up)
    } else{
      roi <- substring(seq, 1, nchar(seq))
    }
    counter_A_T <- 0
    for (i in 1:nchar(roi)) {
      if (substring(roi, i, i) == 'T' | substring(roi, i, i) == 'A') {
        counter_A_T <- counter_A_T + 1
      }
    }
    return(counter_A_T / nchar(roi))
  }
  
  contains_relative <- function(pos, seq) {
    return(pos / nchar(seq))
    # ToADD: Estimate whether the position is considered 'effective' - fist or last quarter of UTR
    # This estimation is currently in the responsibility of the scientist who executes seed matching
  }
  
  contains_Pumilio_count <- function(pos, seq) {
    distance <- 0
    pumi <- unlist(str_locate_all(seq, regex("TGTA[A,C,T,G]ATA")))
    if (length(pumi > 0))
      return(length(pumi / 2))
    else
      return(NA)
  }
  
  contains_Pumilio <- function(pos, seq) {
    distance <- 0
    pumi <- unlist(str_locate_all(seq, regex("TGTA[A,C,T,G]ATA")))
    if (length(pumi > 0))
      return(min(pumi - pos))
    else
      return(NA)
  }
  
  contains_backbone_binding <-
    function(pos, miRNA_transcript, seq) {
      contains <- FALSE
      seed <- upstream_seed(miRNA_transcript)
      limit_upstream_seed <- pos - 5 - 6
      
      if (limit_upstream_seed >= 1) {
        roi <- substring(seq, limit_upstream_seed, pos - 1)
      } else{
        roi <- substring(seq, 1, pos - 1)
      }
      
      pos_second_seed <- list()
      size_second_seed <- list()
      
      for (i in 1:2) {
        positions <-
          unlist(str_locate_all(roi, substring(seed, i, nchar(seed))))
        if (length(positions) > 0) {
          contains <- TRUE
          pos_second_seed <-
            append(pos_second_seed, positions + limit_upstream_seed - 1)
          message(i)
          size_second_seed <-
            append(size_second_seed, nchar(substring(seed, i, nchar(seed))))
        }
        positions <-
          unlist(str_locate_all(roi, substring(seed, i, nchar(seed) - 1)))
        if (length(positions) > 0) {
          contains <- TRUE
          pos_second_seed <-
            append(pos_second_seed, positions + limit_upstream_seed - 1)
          size_second_seed <-
            append(size_second_seed, nchar(substring(seed, i, nchar(seed))))
        }
      }
      
      tmp <- unlist(pos_second_seed)
      tmp2 <- unlist(size_second_seed)
      outlist <- list()
      if (length(tmp) > 0) {
        for (i in 1:length(tmp) / 2) {
          outlist <- append(outlist, paste(tmp[i], tmp2[i], sep = "-"))
        }
      }
      if (length(outlist) > 0)
        return(paste(unlist(outlist), sep = ";"))
      else
        return("-")
    }
  
  contains_prox_same <- function(targetome) {
    
  }
  
  contains_prox_other <- function(targ) {
    matches <- which(!is.na(targ$Pos_8mer) | !is.na(targ$Pos_7mer_m8) | !is.na(targ$Pos_7mer_A1) | !is.na(targ$Pos_6mer))
    for(m in matches){
      mOI <- which(!is.na(targ[m,3:6])) + 2
      value <- targ[m,mOI]
      distance <- 100
      
      if(value+distance <= nchar(process_mRNA(targ$mRNA[m], mRNA))){
        upper <- value + distance
      } else{
        upper <- nchar(process_mRNA(targ$mRNA[m], mRNA))
      }
      
      if(value > distance){
        lower <- value - distance
      } else{
        lower <- 0
      }
      
      others <- which(targ[,3:6] >= lower & targ[,3:6] <= upper)
    }
  }
  
  ### Functions for handling prediction output
  # Add information from targetome to adjacency matrix
  addPredictionToAdjacency <- function(targetome, adj){
    targetome <- targetome[!(is.na(targetome$Pos_8mer) & is.na(targetome$Pos_7mer_m8) & is.na(targetome$Pos_7mer_A1) & is.na(targetome$Pos_6mer)),]
    r_adj <- rownames(adj)
    c_adj <- strsplit(colnames(adj), split ="_")
    # To find matches betweem targetome and adj
    new_c_adj <- character()
    for(i in 1:length(c_adj)){
      new_c_adj <-append(new_c_adj, c_adj[[i]][1])
    }
    
    r <- unique(targetome$miRNA)
    for (miR in r){
      r_pos <- which(r_adj == miR)
      c <- unique(targetome$mRNA[targetome$miRNA == miR])
      for(mR in c){
        c_pos <- which(new_c_adj == mR)
        adj[r_pos, c_pos] <- 2
      }
    }
    return(adj)
  }
  
  # Initialize output table
  targetome <-
    data.frame(
      "miRNA" = as.character(),
      "mRNA" = as.character(),
      "Pos_8mer" = as.numeric(),
      "Pos_7mer_m8" = as.numeric(),
      "Pos_7mer_A1" = as.numeric(),
      "Pos_6mer" = as.numeric(),
      "Pos_offset_6mer" = as.numeric(),
      "AU_content_30nt" = as.numeric(),
      "In_prox_same" = as.numeric(),
      "In_prox_other" = as.numeric(),
      "miR_binding_backbone" = as.character(),
      "relative_within_UTR" = as.numeric(),
      "Pumilio_min_dist" = as.numeric(),
      "Pumilio_count" = as.numeric()
    )
  
  ### Start seed match
  ### START OF LOOPs
  if (complete_miR & !complete_mR) {
    # Check all miRNAs for special mRNA
    message("CASE1: Check all miRNAs for a specific mRNA")
    table_idx = 1
    gene_sequence <- process_mRNA(mRNA_ID, mRNA)
    for (miR in 1:nrow(miRNA)) {
      message(paste(
        "check binding for",
        miRNA$TranscriptID[miR],
        "on",
        mRNA_ID,
        sep = " "
      ))
      seed <- process_seed(miRNA$Sequence[miR])
      
      positions <-
        contains_8mer(miRNA_transcript = seed,
                      mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_8mer[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
      
      positions <-
        contains_7mer_A1(miRNA_transcript = seed,
                         mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_7mer_A1[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
      
      positions <-
        contains_7mer_m8(miRNA_transcript = seed,
                         mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_7mer_m8[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
      
      positions <-
        contains_6mer(miRNA_transcript = seed,
                      mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_6mer[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
      
      positions <-
        contains_offset_6mer(miRNA_transcript = seed,
                             mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_offset_6mer[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
    }
  }
  else if (!complete_miR & complete_mR) {
    # Check all mRNAs for special miRNA
    message("CASE2: Check all mRNAs for a specific miRNA")
    table_idx = 1
    
    for (mR in 1:nrow(mRNA)) {
      gene_sequence <- process_mRNA(mRNA$Gene[mR], mRNA)
      if(mRNA$Gene[mR] %in% targetome$mRNA){
        next
      }
      
      message(paste(
        "check binding for",
        miRNA_ID,
        "on",
        mRNA$Gene[mR],
        sep = " "
      ))
      
      seed <- process_seed(miRNA$Sequence[miRNA$TranscriptID == miRNA_ID])
      
      miR <- which(miRNA$TranscriptID == miRNA_ID)
      mRNA_ID <- mRNA$Gene[mR]
      
      positions <-
        contains_8mer(miRNA_transcript = seed,
                      mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_8mer[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
      
      positions <-
        contains_7mer_A1(miRNA_transcript = seed,
                         mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_7mer_A1[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
      
      positions <-
        contains_7mer_m8(miRNA_transcript = seed,
                         mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_7mer_m8[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
      
      positions <-
        contains_6mer(miRNA_transcript = seed,
                      mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_6mer[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
      
      positions <-
        contains_offset_6mer(miRNA_transcript = seed,
                             mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_offset_6mer[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
    }
  }
  else if (!complete_miR & !complete_mR) {
    # Check specific mRNA with specific miRNA
    message("CASE3: Check specific mRNAs for a specific miRNA")
    table_idx = 1
    gene_sequence <- process_mRNA(mRNA_ID, mRNA)
    message(paste("check binding for", miRNA_ID, "on", mRNA_ID, sep = " "))
    seed <-
      process_seed(miRNA$Sequence[miRNA$TranscriptID == miRNA_ID])
    
    positions <-
      contains_8mer(miRNA_transcript = seed, mRNA_transcript = gene_sequence)
    if (length(positions > 0)) {
      for (p in 1:(length(positions) / 2)) {
        targetome[nrow(targetome) + 1,] <- NA
        targetome$Pos_8mer[table_idx] <- positions[p]
        targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
        targetome$mRNA[table_idx] <- mRNA_ID
        targetome$AU_content_30nt[table_idx] <-
          contains_AU(positions[p], gene_sequence)
        targetome$relative_within_UTR[table_idx] <-
          contains_relative(positions[p], gene_sequence)
        targetome$Pumilio_min_dist[table_idx] <-
          contains_Pumilio(positions[p], gene_sequence)
        targetome$Pumilio_count[table_idx] <-
          contains_Pumilio_count(positions[p], gene_sequence)
        targetome$miR_binding_backbone[table_idx] <-
          contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
        table_idx <- table_idx + 1
      }
    }
    positions <-
      contains_7mer_A1(miRNA_transcript = seed, mRNA_transcript = gene_sequence)
    if (length(positions > 0)) {
      for (p in 1:(length(positions) / 2)) {
        targetome[nrow(targetome) + 1,] <- NA
        targetome$Pos_7mer_A1[table_idx] <- positions[p]
        targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
        targetome$mRNA[table_idx] <- mRNA_ID
        targetome$AU_content_30nt[table_idx] <-
          contains_AU(positions[p], gene_sequence)
        targetome$relative_within_UTR[table_idx] <-
          contains_relative(positions[p], gene_sequence)
        targetome$Pumilio_min_dist[table_idx] <-
          contains_Pumilio(positions[p], gene_sequence)
        targetome$Pumilio_count[table_idx] <-
          contains_Pumilio_count(positions[p], gene_sequence)
        targetome$miR_binding_backbone[table_idx] <-
          contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
        table_idx <- table_idx + 1
      }
    }
    
    positions <-
      contains_7mer_m8(miRNA_transcript = seed, mRNA_transcript = gene_sequence)
    if (length(positions > 0)) {
      for (p in 1:(length(positions) / 2)) {
        targetome[nrow(targetome) + 1,] <- NA
        targetome$Pos_7mer_m8[table_idx] <- positions[p]
        targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
        targetome$mRNA[table_idx] <- mRNA_ID
        targetome$AU_content_30nt[table_idx] <-
          contains_AU(positions[p], gene_sequence)
        targetome$relative_within_UTR[table_idx] <-
          contains_relative(positions[p], gene_sequence)
        targetome$Pumilio_min_dist[table_idx] <-
          contains_Pumilio(positions[p], gene_sequence)
        targetome$Pumilio_count[table_idx] <-
          contains_Pumilio_count(positions[p], gene_sequence)
        targetome$miR_binding_backbone[table_idx] <-
          contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
        table_idx <- table_idx + 1
      }
    }
    
    positions <-
      contains_6mer(miRNA_transcript = seed, mRNA_transcript = gene_sequence)
    if (length(positions > 0)) {
      for (p in 1:(length(positions) / 2)) {
        targetome[nrow(targetome) + 1,] <- NA
        targetome$Pos_6mer[table_idx] <- positions[p]
        targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
        targetome$mRNA[table_idx] <- mRNA_ID
        targetome$AU_content_30nt[table_idx] <-
          contains_AU(positions[p], gene_sequence)
        targetome$relative_within_UTR[table_idx] <-
          contains_relative(positions[p], gene_sequence)
        targetome$Pumilio_min_dist[table_idx] <-
          contains_Pumilio(positions[p], gene_sequence)
        targetome$Pumilio_count[table_idx] <-
          contains_Pumilio_count(positions[p], gene_sequence)
        targetome$miR_binding_backbone[table_idx] <-
          contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
        table_idx <- table_idx + 1
      }
    }
    
    positions <-
      contains_offset_6mer(miRNA_transcript = seed, mRNA_transcript = gene_sequence)
    if (length(positions > 0)) {
      for (p in 1:(length(positions) / 2)) {
        targetome[nrow(targetome) + 1,] <- NA
        targetome$Pos_offset_6mer[table_idx] <- positions[p]
        targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
        targetome$mRNA[table_idx] <- mRNA_ID
        targetome$AU_content_30nt[table_idx] <-
          contains_AU(positions[p], gene_sequence)
        targetome$relative_within_UTR[table_idx] <-
          contains_relative(positions[p], gene_sequence)
        targetome$Pumilio_min_dist[table_idx] <-
          contains_Pumilio(positions[p], gene_sequence)
        targetome$Pumilio_count[table_idx] <-
          contains_Pumilio_count(positions[p], gene_sequence)
        targetome$miR_binding_backbone[table_idx] <-
          contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
        table_idx <- table_idx + 1
      }
    }
    
    
  }
  else {
    # Check all with all
    table_idx = 1
    message("CASE4: Check all mRNAs for all miRNA")
    for (mR in 1:nrow(mRNA)) {
      gene_sequence <- process_mRNA(mRNA$Gene[mR], mRNA)
      mRNA_ID <- mRNA$Gene[mR]
      message(mRNA_ID)
      for (miR in 1:nrow(miRNA)) {
        message(paste(
          "check binding for",
          miRNA$TranscriptID[miR],
          "on",
          mRNA_ID,
          sep = " "
        ))
        seed <- process_seed(miRNA$Sequence[miR])
        
        positions <-
          contains_8mer(miRNA_transcript = seed,
                        mRNA_transcript = gene_sequence)
        if (length(positions > 0)) {
          for (p in 1:(length(positions) / 2)) {
            targetome[nrow(targetome) + 1,] <- NA
            targetome$Pos_8mer[table_idx] <- positions[p]
            targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
            targetome$mRNA[table_idx] <- mRNA_ID
            targetome$AU_content_30nt[table_idx] <-
              contains_AU(positions[p], gene_sequence)
            targetome$relative_within_UTR[table_idx] <-
              contains_relative(positions[p], gene_sequence)
            targetome$Pumilio_min_dist[table_idx] <-
              contains_Pumilio(positions[p], gene_sequence)
            targetome$Pumilio_count[table_idx] <-
              contains_Pumilio_count(positions[p], gene_sequence)
            targetome$miR_binding_backbone[table_idx] <-
              contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
            table_idx <- table_idx + 1
          }
        }
        
        positions <-
          contains_7mer_A1(miRNA_transcript = seed,
                           mRNA_transcript = gene_sequence)
        if (length(positions > 0)) {
          for (p in 1:(length(positions) / 2)) {
            targetome[nrow(targetome) + 1,] <- NA
            targetome$Pos_7mer_A1[table_idx] <- positions[p]
            targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
            targetome$mRNA[table_idx] <- mRNA_ID
            targetome$AU_content_30nt[table_idx] <-
              contains_AU(positions[p], gene_sequence)
            targetome$relative_within_UTR[table_idx] <-
              contains_relative(positions[p], gene_sequence)
            targetome$Pumilio_min_dist[table_idx] <-
              contains_Pumilio(positions[p], gene_sequence)
            targetome$Pumilio_count[table_idx] <-
              contains_Pumilio_count(positions[p], gene_sequence)
            targetome$miR_binding_backbone[table_idx] <-
              contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
            table_idx <- table_idx + 1
          }
        }
        
        positions <-
          contains_7mer_m8(miRNA_transcript = seed,
                           mRNA_transcript = gene_sequence)
        if (length(positions > 0)) {
          for (p in 1:(length(positions) / 2)) {
            targetome[nrow(targetome) + 1,] <- NA
            targetome$Pos_7mer_m8[table_idx] <- positions[p]
            targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
            targetome$mRNA[table_idx] <- mRNA_ID
            targetome$AU_content_30nt[table_idx] <-
              contains_AU(positions[p], gene_sequence)
            targetome$relative_within_UTR[table_idx] <-
              contains_relative(positions[p], gene_sequence)
            targetome$Pumilio_min_dist[table_idx] <-
              contains_Pumilio(positions[p], gene_sequence)
            targetome$Pumilio_count[table_idx] <-
              contains_Pumilio_count(positions[p], gene_sequence)
            targetome$miR_binding_backbone[table_idx] <-
              contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
            table_idx <- table_idx + 1
          }
        }
        
        positions <-
          contains_6mer(miRNA_transcript = seed,
                        mRNA_transcript = gene_sequence)
        if (length(positions > 0)) {
          for (p in 1:(length(positions) / 2)) {
            targetome[nrow(targetome) + 1,] <- NA
            targetome$Pos_6mer[table_idx] <- positions[p]
            targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
            targetome$mRNA[table_idx] <- mRNA_ID
            targetome$AU_content_30nt[table_idx] <-
              contains_AU(positions[p], gene_sequence)
            targetome$relative_within_UTR[table_idx] <-
              contains_relative(positions[p], gene_sequence)
            targetome$Pumilio_min_dist[table_idx] <-
              contains_Pumilio(positions[p], gene_sequence)
            targetome$Pumilio_count[table_idx] <-
              contains_Pumilio_count(positions[p], gene_sequence)
            targetome$miR_binding_backbone[table_idx] <-
              contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
            table_idx <- table_idx + 1
          }
        }
      }
      
      positions <-
        contains_offset_6mer(miRNA_transcript = seed,
                             mRNA_transcript = gene_sequence)
      if (length(positions > 0)) {
        for (p in 1:(length(positions) / 2)) {
          targetome[nrow(targetome) + 1,] <- NA
          targetome$Pos_offset_6mer[table_idx] <- positions[p]
          targetome$miRNA[table_idx] <- miRNA$TranscriptID[miR]
          targetome$mRNA[table_idx] <- mRNA_ID
          targetome$AU_content_30nt[table_idx] <-
            contains_AU(positions[p], gene_sequence)
          targetome$relative_within_UTR[table_idx] <-
            contains_relative(positions[p], gene_sequence)
          targetome$Pumilio_min_dist[table_idx] <-
            contains_Pumilio(positions[p], gene_sequence)
          targetome$Pumilio_count[table_idx] <-
            contains_Pumilio_count(positions[p], gene_sequence)
          targetome$miR_binding_backbone[table_idx] <-
            contains_backbone_binding(positions[p], miRNA$Sequence[miR], gene_sequence)
          table_idx <- table_idx + 1
        }
      }
    }
  }
  
  if(addToAdjacency){
    message("Adapt adjacency matrix")
    adjacency <- addPredictionToAdjacency(targetome, adjacency)
  }
  
  if(writeTo){
    message(paste("Write to file", nameOfTable, sep = " "))
    write.csv(targetome, file = nameOfTable, row.names = FALSE)
  }
  
  message("Generate output")
  output <- list("targetome_table" = targetome, "adjacency_matrix" = adjacency)
  return(output)
}