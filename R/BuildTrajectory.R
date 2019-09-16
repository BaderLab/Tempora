#' Build mutual information-based network
#'
#' Reduce the dimensionality of the pathway enrichment matrix and build the information-based cluster-cluster network
#' @param object A Tempora object containing a gene expression matrix and metadata
#' @param pathwaygmt A database of pathways oirganized as a .gmt file
#' @param method Method used to calculate pathway enrichment profile. Can be "gsva", "ssgsea", "zscore" or "plage". See ?gsva for more information.
#' @param similarity_threshold Percent of similarity for cutoff. Default at 0.01
#' @export
#' @examples tempora_data <- ImportSeuratObject(seurat_object)
#' BuildTrajectory
BuildTrajectory <- function(object, n_pcs, similarity_threshold=0.01){

  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }


  if (n_pcs > ncol(object@cluster.pathways.dr$rotation)){
    stop("Number of PCs selected exceeds number of PCs calculated")
  }

  devtools::use_package("bnlearn")

  significant_pathways_list <- gsva_pca <- list()
  for (i in 1:n_pcs){
    genes_scaled <- scale(object@cluster.pathways.dr$rotation[,i])
    significant_pathways_list[[i]] <- object@cluster.pathways[which(rownames(object@cluster.pathways) %in% names(which(genes_scaled[,1] > 1.5 | genes_scaled[,1] < -1.5))), ]
    gsva_pca[[i]] <- colMeans(significant_pathways_list[[i]])
  }

  gsva_pca <- Reduce(rbind, gsva_pca)
  rownames(gsva_pca) <- paste0("PC", seq(1:nrow(gsva_pca)))

  mi_network <- aracne(as.data.frame(gsva_pca))
  edges_df <- as.data.frame(mi_network$arcs)
  edges_df$to <- as.numeric(as.character(edges_df$to))
  edges_df$from <- as.numeric(as.character(edges_df$from))
  edges_df$from_clusterscore <- unlist(sapply(edges_df$from, function(x) object@cluster.metadata$Cluster_time_score[which(rownames(object@cluster.metadata) == x)]))
  edges_df$to_clusterscore <- unlist(sapply(edges_df$to, function(x) object@cluster.metadata$Cluster_time_score[which(rownames(object@cluster.metadata) == x)]))


  edges_df$direction <- ifelse((abs(edges_df$to_clusterscore - edges_df$from_clusterscore)/(0.5*(edges_df$to_clusterscore + edges_df$from_clusterscore))) < similarity_threshold, "bidirectional", "unidirectional")
  edges_df <- edges_df[-which(edges_df$from_clusterscore > edges_df$to_clusterscore), ]
  edges_df$id <- ifelse(as.numeric(edges_df$from) > as.numeric(edges_df$to), paste0(edges_df$from, edges_df$to), paste0(edges_df$to, edges_df$from))
  edges_df <- edges_df[!duplicated(edges_df$id), ]
  edges_df <- edges_df[, -6]
  edges_df$type <-  ifelse(edges_df$direction == "bidirectional", 3, 1)

  object@trajectory <- edges_df
  object@n.pcs <- n_pcs
  return(object)
}

# testobj <- BuildTrajectory(testobj, n_pcs=2, similarity_threshold = 0.03)


