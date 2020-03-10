#' Calculate pathway enrichment profile
#'
#' Calculate cluster average gene expression profile and determine the pathway enrichment profile of each cluster
#' @param exprMatrix A gene expression matrix, with genes in rows and cells in columns.
#' @param cell_markers A list of possible cell types found in the dataset and their marker genes.
#' @param cluster_labels A named vector of cluster identifier for each cell in the gene expression matrix
#' @export
#' @importFrom GSVA gsva
#' @importFrom ggplot2 ggplot
#' @importFrom reshape2 melt
#' @return A vector of cell types inferred from the expression of marker genes provided
#' CalculatePWProfiles
IdentifyCellTypes <- function(exprMatrix, cluster_labels, cell_markers){
  exprMatrix_bycluster <- list()
  for (i in sort(as.numeric(unique(cluster_labels)))){
    exprMatrix_bycluster[[i]] <- rowMeans(exprMatrix[, which(colnames(exprMatrix) %in% names(cluster_labels)[which(cluster_labels == i)])])
  }
  exprMatrix_bycluster <- do.call(cbind, exprMatrix_bycluster)
  colnames(exprMatrix_bycluster) <- sort(as.numeric(unique(cluster_labels)))
  rownames(exprMatrix_bycluster) <- rownames(exprMatrix)

  cell_type_classifier <- GSVA::gsva(exprMatrix_bycluster, cell_markers, parallel.sz=1)
  cell_type_classifier_m <- reshape2::melt(rownames_to_column(as.data.frame(cell_type_classifier), var="Cell_type"))
  colnames(cell_type_classifier_m)[2] <- "Cluster"
  cell_type_classifier_m$Cluster <- sub("^V", "", cell_type_classifier_m$Cluster)
  cell_type_classifier_m$Cluster <- factor(cell_type_classifier_m$Cluster, levels=seq(1:length(unique(cell_type_classifier_m$Cluster))))

  celltype_plot <- ggplot2::ggplot(cell_type_classifier_m, aes(x=Cluster, y=Cell_type, size=value, color=value)) + geom_point(alpha = 0.8, stroke=0) + scale_x_discrete() + theme_bw()
  celltype_plot <- celltype_plot+scale_size(range = c(1, 12), limits=c(-0.1, 1)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_gradient(low = "grey100",  high = "black", space = "Lab", limit = c(-0.1, 1)) + ggtitle("Cell type prediction scores")
  plot(celltype_plot)

  cell_types <- apply(cell_type_classifier, 2, function(x) paste(rownames(cell_type_classifier)[which(x>0.5)], collapse="/"))
  if (any(cell_types == "")){
    cell_types[cell_types==""] <- "Unclassified"
  }
#   if (any(apply(cell_type_classifier, 2, max) < 0.5)){
#   message(paste("Cluster(s)", names(which(apply(cell_type_classifier, 2, max) < 0.5)), "cannot be classified with confidence using the provided markers and received the label 'Other'. Please revise marker sets and re-run if possible, or
#                 manually alter cluster labels by setting object@clusterlabel"))
#   cell_types[which(apply(cell_type_classifier, 2, max) < 0.5)] <- "Other"
# }
  return(cell_types)
}








