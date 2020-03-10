#' Visualize the trajectory
#'
#' Reduce the dimensionality of the pathway enrichment matrix and build the information-based cluster-cluster network
#' @param object A Tempora object
#' @param layout Layout method for the trajectory plot. Can be "Sugiyama" (hierarchical graph drawing) or "force_directed" (force directed graph drawing)
#' @param ... Any additional arguments to plot.igraph
#' @examples \dontrun{tempora_data <- PlotTrajectory(tempora_data)}
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph layout_with_sugiyama graph_from_data_frame plot.igraph E layout_with_graphopt
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis legend par points text
#' @importFrom methods new validObject
#' @importFrom stats p.adjust prcomp screeplot
#' @importFrom scales rescale
#' @importFrom reshape2 dcast
#'
#'
PlotTrajectory <- function(object, layout, ...){

  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }

  if (is.null(object@trajectory)){
    stop("BuildTrajectory has not been run. See ?Tempora::BuildTrajectory for details")
  }

  if (!layout %in% c("Sugiyama", "force_directed")){
    stop("Layout method not supported. Layout method can be 'Sugiyama' or 'force_directed'")
  }

  edge_graph <- igraph::graph_from_data_frame(d=object@trajectory, vertices = object@cluster.metadata, directed = T)
  if (layout=="Sugiyama") {
    l <- igraph::layout_with_sugiyama(edge_graph, hgap=5, maxiter = 500)
    l$layout[,2] <- 3-(rescale(object@cluster.metadata$Cluster_time_score, to=c(0,3)))
    if (length(levels(object@meta.data$Timepoints)) > 9){
      colours <- colorRampPalette(RColorBrewer::brewer.pal(7, "YlOrRd"))
      plot.igraph(edge_graph, ylim=c(-1,1), layout = l$layout, ylab = "Inferred time", vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                  vertex.pie.color=list(colours(length(levels(object@meta.data$Timepoints)))), pie.border=list(rep("white", 4)), vertex.frame.color="white", edge.arrow.size = 0.8, edge.width = 1.5, vertex.label.family="Arial",
                  vertex.label.color="black", edge.lty = E(edge_graph)$type, vertex.label.cex=1.5, ...)
      axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
      legend("topleft", legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border="black")
    } else {
      colours <- brewer.pal(length(levels(object@meta.data$Timepoints)), "YlOrRd")
      plot.igraph(edge_graph, ylim=c(-1,1), ylab = "Inferred time", layout = l$layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                  vertex.pie.color=list(colours), pie.border=list(rep("white", length(levels(object@meta.data$Timepoints)))), vertex.frame.color="white",
                  edge.arrow.size = 0.8, edge.width = 1.5, vertex.label.family="Arial", vertex.label.color="black", edge.lty = E(edge_graph)$type, vertex.label.cex=1.5, vertex.size = 21, ...)
      legend("topleft", legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border = "black")
      axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
    }
    object@layouts <- l$layout
  } else if (layout=="force_directed") {
    l <- igraph::layout_with_graphopt(edge_graph, spring.length = 10)
    l[,2] <- 3-(scales::rescale(object@cluster.metadata$Cluster_time_score, to=c(0,3)))
    if (length(levels(object@meta.data$Timepoints)) > 9){
      colours <- colorRampPalette(RColorBrewer::brewer.pal(7, "YlOrRd"))
      plot.igraph(edge_graph, ylim=c(-1,1), layout = l, ylab = "Inferred time", vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                  vertex.pie.color=list(colours(length(levels(object@meta.data$Timepoints)))), pie.border=list(rep("white", 4)), vertex.frame.color="white", edge.arrow.size = 0.8, edge.width = 1.5, vertex.label.family="Arial",
                  vertex.label.color="black", edge.lty = E(edge_graph)$type,  ...)
      axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
      legend("topleft", legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border="black")
    } else {
      colours <- brewer.pal(length(levels(object@meta.data$Timepoints)), "YlOrRd")
      plot.igraph(edge_graph, ylim=c(-1,1), ylab = "Inferred time", layout = l, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                  vertex.pie.color=list(colours), pie.border=list(rep("white", length(levels(object@meta.data$Timepoints)))), vertex.frame.color="white",
                  edge.arrow.size = 0.8, edge.width = 1.5, vertex.label.family="Arial", vertex.label.color="black", edge.lty = E(edge_graph)$type, vertex.size = 21, ...)
      legend("topleft", legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border = "black")
      axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
    }
    object@layouts <- l
  }

  validObject(object)
  return(object)
}



