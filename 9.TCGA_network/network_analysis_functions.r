# Functions for protein-protein interactions 

iid_remove_duplicates <- function(annot_file, protein_list){
  uniprot1 <- annot_file[which(as.character(annot_file$uniprot1) %in% proteins_list),]
  uniprot2 <- annot_file[which(as.character(annot_file$uniprot2) %in% proteins_list),]
  iid <- merge(uniprot1[,1:6], uniprot2[,1:6], all = TRUE)
  onetwo <- paste(iid$uniprot1, iid$uniprot2)
  twoone <- paste(iid$uniprot2, iid$uniprot1)
  iid_dup <- iid[-which(onetwo %in% twoone),]
  return(iid_dup)
  
}

iid_remove_with_few_pmid <- function(iid_dup){
  iid_pmid <- data.frame()
  for(i in seq(1:length(rownames(iid_dup)))){
    pmid <- length(strsplit(as.character(iid_dup[i,]$pmids), ";")[[1]])
    if(pmid > 1){
      d_f <- data.frame(iid_dup[i,])
      iid_pmid <- rbind(iid_pmid, d_f)
    }
  }
  return(iid_pmid)
}

iid_remove_same_gene1_gene2 <- function(iid_pmid){
  remove_row_list <- c()
  for(row_numb in seq(1,length(rownames(iid_pmid)))){
    if(as.character(iid_pmid[row_numb,]$symbol1) == as.character(iid_pmid[row_numb,]$symbol2)){
      remove_row_list <- c(remove_row_list, row_numb)
    }
  }
  iid_same <- iid_pmid[-remove_row_list,]
  write.csv(iid_same, "protein_interactions_filt.csv")
  return(iid_same)
}


iid_get_ens_DEG <- function(iid_pmid, cancer_types){
  iid_pmid <- iid_pmid[,c(1,2,3,4)]
  iid_pmid$uniprot1 <- as.character(iid_pmid$uniprot1)
  iid_pmid$uniprot2 <- as.character(iid_pmid$uniprot2)
  uni_all <- data.frame()
  
  for(cancer in cancer_types){
    dir.create("protein_interactions", showWarnings = FALSE)
    cat("Working on ", cancer, "\n")
    
    # Get gdc data from TCGA_download 
    gdc <- get(load(paste0("../1.TCGA_download_prepare/data/",tolower(cancer),".exp.rda")))
    
    # Get differential expression analysis data
    dataDEGs <- get(load(paste0("../3.TCGA_DEA/DEA/",cancer, "_dataDEGs.rda")))
    
    # Change rownames from HUGO gene names to ensembl ids
    rownames(dataDEGs) <- rowData(gdc)[match(rownames(dataDEGs), rowData(gdc)[,"external_gene_name"]),"ensembl_gene_id"]
    dataDEGs <- as.data.frame(dataDEGs)
    
    # Get matching uniprot id for every ensembl id
    mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    bm_file <- getBM(filters= "ensembl_gene_id", 
                     attributes= c("ensembl_gene_id","uniprotswissprot"),
                     values= rownames(dataDEGs),
                     mart= mart)
    
    ens_pro <- bm_file[-which(bm_file$uniprotswissprot == ""),]

    uniprot_df <- as.data.frame(ens_pro$uniprotswissprot)
    colnames(uniprot_df) <- "uniprot_id"
    uniprot_df$ensembl_id <- rep("NA", length(rownames(uniprot_df)))
    uniprot_df$DEA <- rep("NA", length(rownames(uniprot_df)))
    uniprot_df$interaction_uniprot <- rep("NA", length(rownames(uniprot_df)))
    iid_pmid_uniprots <- c(iid_pmid$uniprot1, iid_pmid$uniprot2)
    iid_pmid_uniprots <- unique(iid_pmid_uniprots)
    for(uniprots in seq(1, length(iid_pmid_uniprots))){
      int_1_2 <- ""
      int_1_2 <- ""
      int <- iid_pmid_uniprots[uniprots]
      int_1 <- which(iid_pmid$uniprot1 == int)
      int_2 <- which(iid_pmid$uniprot2 == int)
      int_1_2 <- iid_pmid[int_1,]$uniprot2
      int_2_1 <- iid_pmid[int_2,]$uniprot1
      all_ints <- c(int_1_2, int_2_1)
      all_ints <- unique(all_ints)
      all_interactions <- paste(all_ints, collapse = ";")
      uni_df_pos <- which(uniprot_df$uniprot_id == int)
      if(as.numeric(length(uni_df_pos)) > 0){
        uniprot_df[which(uniprot_df$uniprot_id == int),]$interaction_uniprot <- all_interactions
      }
    }
    
    uni_prot <- uniprot_df[-which(uniprot_df$interaction_uniprot == "NA"),]
    for(unis in uni_prot$uniprot_id){
      lengen <- ens_pro[which(ens_pro$uniprotswissprot == unis),]$ensembl_gene_id
      if(length(lengen) > 1){
        ensembl_DE <- c()
        lengen <- sort(lengen)
        lengen_l <- paste(lengen, collapse = ";")
        uni_prot[which(uni_prot$uniprot_id == unis),]$ensembl_id <- lengen_l
        for(ens_long in lengen){
          ens_DE <- dataDEGs[which(rownames(dataDEGs) == ens_long),]$logFC
          ensembl_DE <- c(ensembl_DE, ens_DE)
        }
        ens_DEG <- paste(ensembl_DE, collapse = ";")
        uni_prot[which(uni_prot$uniprot_id == unis),]$DEA <- ens_DEG
      }
      else if(length(lengen) == 1){
        uni_prot[which(uni_prot$uniprot_id == unis),]$ensembl_id <- lengen
        ens_DEG <- dataDEGs[which(rownames(dataDEGs) == lengen),]$logFC
        uni_prot[which(uni_prot$uniprot_id == unis),]$DEA <- ens_DEG
      }
    }
    
    uni_prot$cancer <- rep(cancer, length(rownames(uni_prot)))
    uni_all <- rbind(uni_all, uni_prot)
  }
  save(uni_all, file = "protein_interactions/protein_interactions_DEA.rda")
  return(uni_all)
}

iid_get_multiple_ens <- function(iid_ens){
  iid_ens_short <- data.frame()
  for(unipr in seq(1,length(rownames(iid_ens)))){
    numbers_ens <- as.numeric(length(strsplit(as.character(iid_ens[unipr,]$ensembl_id), ";")[[1]]))
    if(numbers_ens > 1){
      longlong <- iid_ens[unipr,]
      iid_ens_short <- rbind(iid_ens_short,longlong)
    }
  }
  uniq_prot <- unique(as.character(iid_ens_short$uniprot_id))
  uniq_prot <- as.data.frame(uniq_prot)
  colnames(uniq_prot) <- "uniprot_id"
  uniq_prot$ensembl_id <- rep("NA", length(rownames(uniq_prot))) 
  for(uni in uniq_prot$uniprot_id){
    protein_ids_pos <- which(iid_ens_short$uniprot_id == uni)
    ensembls <- iid_ens_short[protein_ids_pos,]$ensembl_id
    ens_list <- strsplit(ensembls, ";")[[1]]
    uni_ens_list <- unique(ens_list)
    uni_list <- paste(uni_ens_list, collapse = ";")
    uniq_prot[which(uniq_prot$uniprot_id == uni),]$ensembl_id <- uni_list
  }
  
  return(uniq_prot)
}

iid_get_single_ens <- function(iid_ens, ensembl_remove){
  iid_single <- iid_ens
  for(unipr in seq(1,length(rownames(iid_ens)))){
    get_all_ens <- strsplit(as.character(iid_ens[unipr,]$ensembl_id), ";")[[1]]
    find_remove_ens <- which(get_all_ens %in% ensembl_remove)
    if(length(find_remove_ens) > 0){
      get_all_deg <- strsplit(as.character(iid_ens[unipr,]$DEA), ";")[[1]]
      remove_deg <- get_all_deg[-find_remove_ens]
      remove_ens <- get_all_ens[-find_remove_ens]
      if(length(remove_deg) < 1){
        remove_deg <- "NA"
      }
      if(length(remove_ens) < 1){
        remove_ens <- "NA"
      }
      
      iid_single$ensembl_id[unipr] <- remove_ens
      iid_single$DEA[unipr] <- remove_deg
      
    }
  }
  nans <- which(iid_single$ensembl_id == "NA")
  if(length(nans) > 0){
    iid_single <- iid_single[- nans,]
  }

  return(iid_single)
}

prepare_for_networkplot <- function(iid_single, iid_pmid, cancer_types){
  iid_network <- iid_pmid[,c(1,2,3,4)]
  iid_network$uniprot1 <- as.character(iid_network$uniprot1)
  iid_network$uniprot2 <- as.character(iid_network$uniprot2)
  for(cancer in cancer_types){
    deg1 <- as.data.frame(rep(0, length(rownames(iid_network))))
    colnames(deg1) <- c(paste0(cancer, "_DEG1", sep = ""))
    deg2 <- as.data.frame(rep(0, length(rownames(iid_network))))
    colnames(deg2) <- c(paste0(cancer, "_DEG2", sep = ""))
    iid_network <- cbind(iid_network, deg1)
    iid_network <- cbind(iid_network, deg2)
    iid_cancer <- iid_single[which(iid_single$cancer == cancer),]
    for(prot_row in seq(1, length(rownames(iid_cancer)))){
      protein_id <- iid_cancer[prot_row,]$uniprot_id
      diffexp <- iid_cancer[prot_row,]$DEA
      
      protein1 <- which(iid_network$uniprot1 == protein_id)
      iid_network[protein1,length(colnames(iid_network))-1] <- diffexp
      
      protein2 <- which(iid_network$uniprot2 == protein_id)
      iid_network[protein2,length(colnames(iid_network))] <- diffexp
    }
  }
  return(iid_network)
}

do_network_plot <- function(iid_network){
  dir.create("plots", showWarnings = FALSE)
  for(col_number in seq(5,length(colnames(iid_network)),2)){
    pos1 <- col_number
    pos2 <- col_number + 1
    iid_plot <- iid_network[,c(3,4,pos1, pos2)]
    iid_plot[,3] <- as.numeric(iid_plot[,3])
    iid_plot[,4] <- as.numeric(iid_plot[,4])
    del_row <- c()
    for(rownumb in seq(1,length(rownames(iid_plot)))){
      if(iid_plot[rownumb,3] == 0 || iid_plot[rownumb,4] == 0){
        del_row <- c(del_row, rownumb)
      } 
    }
    iid_plot <- iid_plot[-del_row,]
    
    colour_list <- c()
    for(val in seq(1,length(rownames(iid_plot)))){
      if((iid_plot[val,3] >= 0) && (iid_plot[val,4] >= 0)){
        colour_list <- c(colour_list, "purple")
      }
      else if((iid_plot[val,3] < 0) && (iid_plot[val,4] < 0)){
        colour_list <- c(colour_list, "purple")
      }
      else if((iid_plot[val,3] < 0 ) && (iid_plot[val,4] >= 0)){
        colour_list <- c(colour_list, "red3")
      }
      else if((iid_plot[val,3] >= 0) && (iid_plot[val,4] < 0)){
        colour_list <- c(colour_list, "blue4")
      }
    } 
    cancer <- strsplit(as.character(colnames(iid_plot)[3]), "_")[[1]][1]
    corr_graph <- graph_from_data_frame(iid_plot)
    pdf(paste0("plots/",cancer,"_network_plot.pdf", sep = ""))
    plot(corr_graph, edge.color = colour_list, main = paste0(cancer, " network plot"), edge.arrow.size=0.3, vertex.label.cex = 0.3, vertex.size = 10)
    dev.off()
  }
}
