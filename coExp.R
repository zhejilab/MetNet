library("optparse")

option_list = list(
  make_option(c("-d", "--microarray_data"), type="character",
              help="Microarray RMA values", metavar="character"),
  make_option(c("-g", "--genes"), type="character", 
              help="Intersted set of genes", metavar="character"),
  make_option(c("-p", "--patience"), type="character", 
              help="Termination Patience [default = 25]", metavar="character"),
  make_option(c("-s", "--stabilization_threshold"), type="character", 
              help="Stabilization Threshold [default = 0.01]", metavar="character"),
  make_option(c("-f", "--fdr"), type="character",
              help="FDR [default = 0.01]", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

mat2dist <- function(input_matrix){  
  d_mat <- dist(input_matrix, method = "euclidean")
  d_list0 <- unlist(as.data.frame(as.matrix(d_mat)))
  
  zero_num <- table(d_list0 == 0)['TRUE']
  d_list0 <- d_list0[! d_list0 == 0]
  d_list <- c(d_list0, rep(0, zero_num - nrow(input_matrix)))
  
  return (unname(d_list))
}

shuffle_cor_mat <- function(input_matrix){
  temp <- sample(unlist(as.data.frame(input_matrix)))
  return (matrix(temp, ncol = nrow(input_matrix)))
}

cut_dend_tree <- function(data_frame, cutoff){
  d_cor <- dist(data_frame,method = "euclidean") 
  hc_cor <- hclust(d_cor, method = "complete")
  return (cutree(hc_cor, h=cutoff))}

auto_shuffle <- function(input_matrix, shuffle_num = 1000, fdr, break_threshold = 0.01, patience=25){
  
  ini <- c()
  qile <- c()
  
  for (shuffle_i in 1:shuffle_num){
    
    rand_mat <- cor(t(pc_data[sample(rownames(pc_data), nrow(input_matrix)), ]), method = "spearman")
    ini <- c(ini, mat2dist(rand_mat))
    qile <- c(qile, quantile(ini, 1-fdr))
    
    if (shuffle_i > patience+1){
      check_vector <- qile[(length(qile)-patience+1):length(qile)]
      recent_range <- max(check_vector) - min(check_vector)
      if (recent_range < break_threshold){
        break
      }
    }
  }
  return (qile)
}


cuttree_k2 <- function(data_frame){
  d_cor <- dist(data_frame, method = "euclidean")
  hc_cor <- hclust(d_cor, method = "complete")
  return (cutree(hc_cor,k=2))}

cor_matrix2percentile <- function(input_matrix, percentile){  #excluding diagonal values
  d_list0 <- unlist(as.data.frame(as.matrix(input_matrix)))
  
  zero_num <- table(d_list0 == 1)['TRUE']
  d_list0 <- d_list0[! d_list0 == 1]
  d_list <- c(d_list0, rep(1, zero_num - nrow(input_matrix)))
  
  return (quantile(d_list, percentile))
}

recursion_func <- function(mean_cen_pc_data, gene_list, threshold1=0.017, threshold2=0.024){
  
  temp_pc_df <- mean_cen_pc_data[rownames(mean_cen_pc_data) %in% gene_list, ]
  temp_cor <- cor(t(temp_pc_df), method = "spearman")
  temp_cor_percentile1 <- cor_matrix2percentile(temp_cor, 0.1)
  temp_cor_percentile2 <- cor_matrix2percentile(temp_cor, 0.5)
  
  if (temp_cor_percentile1 > threshold1 & temp_cor_percentile2 > threshold2){
    record_list <<- c(record_list, paste(gene_list, collapse = "--"))
  }else{
    new_cut <- cuttree_k2(temp_cor)
    for (k_i in 1:2){
      temp_new_gene_list <<- names(new_cut[new_cut==k_i])
      if (length(temp_new_gene_list) >= 10){
        recursion_func(mean_cen_pc_data, temp_new_gene_list, threshold1, threshold2)
      }else{
        droppped_list <<- c(droppped_list, temp_new_gene_list)
      }
    }
  }
}

data <- read.csv(option_list[[1]], row.names = 1, check.names = F)
mean_cen_data  <- data - rowMeans(data)
mean_cen_data 

cor_mat <- cor(t(pc_data[rownames(pc_data) %in% option_list[[2]], ]), method = "spearman")
dim(cor_mat)

fdrs <- auto_shuffle(cor_mat, fdr = 0.01, shuffle_num = 1000, fdr, break_threshold = option_list[[3]], patience=option_list[[4]])

subtree <- cut_dend_tree(cor_mat, fdrs[length(fdrs)])

gene_clusters <- c()

for (i in unique(subtree)){
  record_list <- c()
  droppped_list <- c()
  recursion_func(mean_cen_pc_data, subgos[[GO_ID]], threshold1, threshold2)

  if (length(droppped_list)!= 0){
    names(droppped_list) <- rep("Removed", length(droppped_list))
  }
  
  num <- 1
  for (i in record_list){
    temp_cluster <- unlist(str_split(i, "--"))
    names(temp_cluster) <- rep(paste(num), length(temp_cluster))
    gene_clusters <- c(gene_clusters, temp_cluster)
    num <- num + 1
  }

  gene_clusters <- c(gene_clusters, droppped_list)
  write.csv(data.frame(id = gene_clusters, cluster = names(gene_clusters)), paste0(option_list[[5]], "/tree_", i, ".csv"))

}
