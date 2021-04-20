#' Build mutual information-based network
#'
#' Build the information-based cluster-cluster network using reduced-dimensionality pathway enrichment profiles of all clusters.
#' This network connects related cell types and states across multiple time points. Taking advantage of the available time information, Tempora also infers the directions of all connections in a trajectory that go from early to late clusters.
#' @param object A Tempora object containing a gene expression matrix and metadata
#' @param n_pcs Number of principal components to be used in building the network.
#' @param difference_threshold Percent of permissible difference between the temporal scores of two clusters to determine the direction of their connections. The temporal scores are calculated based on based on the clusters' composition of cells from each timepoint. The directions of edges connecting pairs of clusters will only be determined for cluster pairs with difference in their time scores higher than the threshold. Other edges will remain undirected. Default at 0.01
#' @param loadings Threshold of PCA loadings for pathways to be used in trajectory construction. The higher the loading, the more the pathway contributes to a principal component. Default at 0.4.
#' @export
#' @importFrom bnlearn aracne
#' @importFrom methods new validObject
#' @importFrom stats p.adjust prcomp screeplot
#' @importFrom reshape2 dcast
#' @examples \dontrun{tempora_data <- BuildTrajectory(tempora_data, n_pcs=10, difference_threshold=0.01, loadings=0.4)}
#' BuildTrajectory
BuildTrajectory <- function(object, n_pcs, difference_threshold=0.01, loadings = 0.4){

  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }


  if (n_pcs > ncol(object@cluster.pathways.dr$rotation)){
    stop("Number of PCs selected exceeds number of PCs calculated")
  }

  significant_pathways_list <- gsva_pca <- list()
  for (i in 1:n_pcs){
    genes_scaled <- scale(object@cluster.pathways.dr$rotation[,i])
    significant_pathways_list[[i]] <- object@cluster.pathways[which(rownames(object@cluster.pathways) %in% names(which(genes_scaled[,1] > loadings | genes_scaled[,1] < (-1*loadings)))), ]
    gsva_pca[[i]] <- colMeans(significant_pathways_list[[i]])
  }

  gsva_pca <- Reduce(rbind, gsva_pca)
  rownames(gsva_pca) <- paste0("PC", seq(1:nrow(gsva_pca)))

  mi_network <- bnlearn::aracne(as.data.frame(gsva_pca))
  edges_df <- as.data.frame(mi_network$arcs)
  edges_df$to <- as.numeric(as.character(edges_df$to))
  edges_df$from <- as.numeric(as.character(edges_df$from))
  edges_df$from_clusterscore <- unlist(sapply(edges_df$from, function(x) object@cluster.metadata$Cluster_time_score[object@cluster.metadata$Id == x]))
  edges_df$to_clusterscore <- unlist(sapply(edges_df$to, function(x) object@cluster.metadata$Cluster_time_score[object@cluster.metadata$Id == x]))


  edges_df$direction <- ifelse((abs(edges_df$to_clusterscore - edges_df$from_clusterscore)/(0.5*(edges_df$to_clusterscore + edges_df$from_clusterscore))) < difference_threshold, "bidirectional", "unidirectional")
  edges_df <- edges_df[-which(edges_df$from_clusterscore > edges_df$to_clusterscore), ]
  edges_df$id <- ifelse(as.numeric(edges_df$from) > as.numeric(edges_df$to), paste0(edges_df$from, edges_df$to), paste0(edges_df$to, edges_df$from))
  edges_df <- edges_df[!duplicated(edges_df$id), ]
  edges_df <- edges_df[, -6]
  edges_df$type <-  ifelse(edges_df$direction == "bidirectional", 3, 1)

  object@trajectory <- edges_df
  object@n.pcs <- n_pcs
  return(object)
}



