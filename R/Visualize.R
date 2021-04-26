#' Visualize the trajectory
#'
#' Reduce the dimensionality of the pathway enrichment matrix and build the information-based cluster-cluster network
#' @param object A Tempora object
#' @param layout A 2-column matrix containing the x- and y-coordinates of all vertices in the graph. If NULL, the layout will be obtained using Sugiyama layout
#' @param ... Any additional arguments to plot.igraph
#' @examples \dontrun{tempora_data <- PlotTrajectory(tempora_data)}
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph layout_with_sugiyama graph_from_data_frame plot.igraph E
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis legend par points text
#' @importFrom methods new validObject
#' @importFrom stats p.adjust prcomp screeplot
#' @importFrom scales rescale
#' @importFrom reshape2 dcast
#' @importFrom magrittr '%>%'
#'
#'
PlotTrajectory <- function(object, layout=NULL, ...){

  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }

  if (is.null(object@trajectory)){
    stop("BuildTrajectory has not been run. See ?Tempora::BuildTrajectory for details")
  }

  edge_graph <- igraph::graph_from_data_frame(d=object@trajectory, vertices = object@cluster.metadata, directed = T)

  if (is.null(layout)){
    l <- igraph::layout_with_sugiyama(edge_graph, layers = object@cluster.metadata$Cluster_time_score, maxiter = 1000)
    #l$layout[,2] <- 3-(rescale(object@cluster.metadata$Cluster_time_score, to=c(0,3)))
    if (length(levels(object@meta.data$Timepoints)) > 9){
      colours <- colorRampPalette(RColorBrewer::brewer.pal(7, "YlOrRd"))
      plot.igraph(edge_graph, ylim=c(-1,1), layout = l$layout, ylab = "Inferred time", vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                  vertex.pie.color=list(colours(length(levels(object@meta.data$Timepoints)))), pie.border=list(rep("white", 4)), vertex.frame.color="white", edge.arrow.size = 0.5, edge.width = 1.5, vertex.label.family="Arial",
                  vertex.label.color="black", edge.lty = E(edge_graph)$type, ...)
      axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
      legend("topleft", legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border="black")
    } else {
      colours <- brewer.pal(length(levels(object@meta.data$Timepoints)), "YlOrRd")
      plot.igraph(edge_graph, ylim=c(-1,1), ylab = "Inferred time", layout = l$layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                  vertex.pie.color=list(colours), pie.border=list(rep("white", length(levels(object@meta.data$Timepoints)))), vertex.frame.color="white",
                  vertex.label.family="Arial", vertex.label.color="black", edge.lty = E(edge_graph)$type,...)
      legend("topleft", legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border = "black")
      axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
    }
    object@layouts <- l$layout

  } else {
    if (length(levels(object@meta.data$Timepoints)) > 9){
      colours <- colorRampPalette(RColorBrewer::brewer.pal(7, "YlOrRd"))
      plot.igraph(edge_graph, ylim=c(-1,1), layout = layout, ylab = "Inferred time", vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                  vertex.pie.color=list(colours(length(levels(object@meta.data$Timepoints)))), pie.border=list(rep("white", 4)), vertex.frame.color="white", edge.arrow.size = 0.5, edge.width = 1.5, vertex.label.family="Arial",
                  vertex.label.color="black", edge.lty = E(edge_graph)$type, ...)
      axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
      legend("topleft", legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border="black")
    } else {
      colours <- brewer.pal(length(levels(object@meta.data$Timepoints)), "YlOrRd")
      plot.igraph(edge_graph, ylim=c(-1,1), ylab = "Inferred time", layout = layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), function(x) as.numeric(object@cluster.metadata[x,2:((length(levels(object@meta.data$Timepoints)))+1)])),
                  vertex.pie.color=list(colours), pie.border=list(rep("white", length(levels(object@meta.data$Timepoints)))), vertex.frame.color="white",
                  vertex.label.family="Arial", vertex.label.color="black", edge.lty = E(edge_graph)$type,...)
      legend("topleft", legend = levels(object@meta.data$Timepoints), fill=colours, bty = "n", border = "black")
      axis(side=2, at=c(-1,1), labels=c("Late","Early"), las=1)
    }
  }

  validObject(object)
  return(object)
}



