# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


###define S4 class
Tempora <- setClass(
  "Tempora",
  slots = c(
    data = "matrix",
    meta.data = "data.frame",
    cluster.metadata = "data.frame",
    cluster.pathways = "matrix",
    cluster.pathways.dr = "ANY",
    n.pcs = "numeric",
    trajectory = "data.frame",
    layouts = "matrix"
  ))
setValidity("Tempora", function(object)
  {
    if(nrow(object@meta.data) != ncol(object@data)){
      return("The numbers of cells in the expression matrix and metadata are different")
    }
    if (any(!rownames(object@meta.data) %in% colnames(object@data))){
      return("Cell names in the expression matrix and metadata are different")
    }
    return(TRUE)
  }
)


###accessors

setGeneric("data", function(x) standardGeneric("data"))
#' @export
setGeneric("data<-", function(x, matrix) standardGeneric("data<-"))
setMethod("data", "Tempora", function(x) x@data)
setMethod("data<-", "Tempora", function(x, matrix) {
  x@data <- matrix
  validObject(x)
  return(x)
})

setGeneric("meta.data", function(x) standardGeneric("meta.data"))
#' @export
setGeneric("meta.data<-", function(x, data.frame) standardGeneric("meta.data<-"))
setMethod("meta.data", "Tempora", function(x) x@meta.data)
setMethod("meta.data<-", "Tempora", function(x, data.frame) {
  x@meta.data <- data.frame
  validObject(x)
  return(x)
})

setGeneric("cluster.metadata", function(x) standardGeneric("cluster.metadata"))
#' @export
setGeneric("cluster.metadata<-", function(x, data.frame) standardGeneric("cluster.metadata<-"))
setMethod("cluster.metadata", "Tempora", function(x) x@cluster.metadata)
setMethod("cluster.metadata<-", "Tempora", function(x, data.frame) {
  x@cluster.metadata <- data.frame
  validObject(x)
  return(x)
})

setGeneric("cluster.pathways", function(x) standardGeneric("cluster.pathways"))
#' @export
setGeneric("cluster.pathways<-", function(x, matrix) standardGeneric("cluster.pathways<-"))
setMethod("cluster.pathways", "Tempora", function(x) x@cluster.pathway)
setMethod("cluster.pathways<-", "Tempora", function(x, matrix) {
  x@cluster.pathway <- matrix
  validObject(x)
  return(x)
})

setGeneric("cluster.pathways.dr", function(x) standardGeneric("cluster.pathways.dr"))
#' @export
setGeneric("cluster.pathways.dr<-", function(x, value) standardGeneric("cluster.pathways.dr<-"))
setMethod("cluster.pathways.dr", "Tempora", function(x) x@cluster.pathways.dr)
setMethod("cluster.pathways.dr<-", "Tempora", function(x, value) {
  x@cluster.pathways.dr <- value
  validObject(x)
  return(x)
})

setGeneric("n.pcs", function(x) standardGeneric("n.pcs"))
#' @export
setGeneric("n.pcs<-", function(x, value) standardGeneric("n.pcs<-"))
setMethod("n.pcs", "Tempora", function(x) x@cluster.pathways.dr)
setMethod("n.pcs<-", "Tempora", function(x, value) {
  x@n.pcs <- value
  validObject(x)
  return(x)
})

setGeneric("trajectory", function(x) standardGeneric("trajectory"))
#' @export
setGeneric("trajectory<-", function(x, data.frame) standardGeneric("trajectory<-"))
setMethod("trajectory", "Tempora", function(x) x@trajectory)
setMethod("trajectory<-", "Tempora", function(x, data.frame) {
  x@trajectory <- data.frame
  validObject(x)
  return(x)
})

setGeneric("layouts", function(x) standardGeneric("layouts"))
#' @export
setGeneric("layouts<-", function(x, data.frame) standardGeneric("layouts<-"))
setMethod("layouts", "Tempora", function(x) x@layouts)
setMethod("layouts<-", "Tempora", function(x, data.frame) {
  x@layouts <- data.frame
  validObject(x)
  return(x)
})


###Create new Tempora object from expression matrix
#' Create new Tempora object from a processed gene expression matrix
#' @param exprMatrix A normalized gene expression matrix containing cells from all timepoints of the time-series study. Batch effect correction is highly recommended before normalization.
#' @param meta.data A dataframe of meta data for all cells in the expression matrix. Each column is a feature and each row stores single cell information. At minimum, this dataframe should contain two columns: a "Clusters" column storing the clustering identity and a "Timepoints" column storing the timepoint when each cell comes from
#' @param timepoint_order An ordered vector of timepoint names from early to late
#' @param cluster_labels A vector of cluster annotations (cell types, cell states, cell cycles, etc.), ordered alphanumerically by cluster names. If NULL, cluster numbers will be used to label the trajectory plot
#' @export
#' @return A Tempora object containing the expression matrix and metadata
#' @examples tempora_dara <- CreateTemporaObject(exprMatrix, meta.data)
#'
CreateTemporaObject <- function(exprMatrix, meta.data, timepoint_order, cluster_labels=NULL){

  if (!'Timepoints' %in% colnames(meta.data)){
    stop("meta.data needs to contain a column named 'Timepoints' for downstream analyses")
  }
  if (!'Clusters' %in% colnames(meta.data)){
    stop("meta.data needs to contain a column named 'Clusters' for downstream analyses")
  }
  if (!is.numeric(exprMatrix)) {
    stop("Expression matrix is not numeric")
  }
  if (any(!meta.data$Timepoints %in% timepoint_order)){
    stop("List of timepoints does not match the timepoints in the data")
  }

  meta.data$Timepoints <- factor(meta.data$Timepoints, levels = timepoint_order)
  meta.data$Timescore <- as.integer(meta.data$Timepoints)

  clustmd <- meta.data[, c("Timepoints", "Clusters")]
  clustmd <- dcast(clustmd, Clusters~Timepoints, value.var = "Clusters", fun.aggregate = length)
  clustmd[, 2:ncol(clustmd)] <- t(apply(clustmd[, 2:ncol(clustmd)], 1, function(x) x/sum(x)))
  clustmd$Cluster_time_score <- apply(clustmd[, 2:ncol(clustmd)], 1,
                                      function(x) sum(mapply(function(t, y) t*y, as.numeric(x), sort(unique(meta.data$Timescore), decreasing = F))))
  colnames(clustmd)[1] <- "Id"
  if (!is.null(cluster_labels)){
    clustmd$label <- paste0("Cluster ", paste(rownames(testclustmd), cluster_labels, sep="-"))
  } else {
    clustmd$label <- paste("Cluster ", rownames(testclustmd))
  }

  tempora <- new("Tempora",
                 data = exprMatrix,
                 meta.data = meta.data,
                 cluster.metadata = clustmd)

  validObject(tempora)
  return(tempora)
}


#' Import data from a Seurat object
#'
#' Imports gene expression matrix and other metadata from a seurat object
#' @param seuratobj A Seurat object containing the normalized gene expression matrix and clustering result
#' @param clusters Name of the column in the meta.data dataframe containing the cluster identity of all cells in the dataset
#' @param timepoints Name of the column in the meta.data dataframe containing the collection time of all cells in the dataset
#' @param timepoint_order An ordered vector of timepoint names from early to late
#' @param cluster_labels A vector of cluster annotations (cell types, cell states, cell cycles, etc.), ordered alphanumerically by cluster names. If NULL, cluster numbers will be used to label the trajectory plot
#' @export
#' @return A Tempora object containing the expression matrix and metadata
#' @examples tempora_data <- ImportSeuratObject(seurat_object, clusters = "res.0.3", timepoints = "collection_time", timepoint_order = c("0H", "24H", "48H", "72H"), cluster_labels = c("Stem cells", "Differentiated cells"))

#' Import data from a Seurat object
ImportSeuratObject <- function(seuratobj, clusters, timepoints, timepoint_order, cluster_labels){
  if(class(seuratobj)[1]=='seurat'){
    requireNamespace("Seurat")
  } else {
    stop("Not a Seurat object. Tempora only supports importing Seurat objects at the moment. See ?Tempora::CreateTemporaObject to manually create a Tempora object from an expression matrix")
  }
  data <- seuratobj@data
  cat("Extracting data...")
  metadata <- seuratobj@meta.data
  cat("\nExtracting metadata...")
  colnames(metadata)[which(colnames(metadata)==clusters)] <- "Clusters"
  colnames(metadata)[which(colnames(metadata)==timepoints)] <- "Timepoints"
  cat("\nCreating Tempora object...")
  tempora_obj <- CreateTemporaObject(data, meta.data = metadata, timepoint_order =  timepoint_order, cluster_labels = cluster_labels)
  validObject(tempora_obj)
  return(tempora_obj)
}


### test ###
# load("~/Desktop/HSMM_seurat_aligned.RData")
# testobj <- ImportSeuratObject(eb1S, clusters = "res.0.6", timepoints = "orig.ident", timepoint_order=c("T0", "T24", "T48", "T72"), cluster_labels = c("A", "B", "C"))
#
#
# meta.data <- testobj@meta.data
# clustmd <- meta.data[, c("Timepoints", "Clusters")]
# clustmd <- dcast(clustmd, Clusters~Timepoints, value.var = "Clusters", fun.aggregate = length)
# clustmd[, 2:ncol(clustmd)] <- t(apply(clustmd[, 2:ncol(clustmd)], 1, function(x) x/sum(x)))
# clustmd$Cluster_time_score <- apply(clustmd[, 2:ncol(clustmd)], 1,
#                                     function(x) sum(mapply(function(t, y) t*y, as.numeric(x), sort(unique(meta.data$Timescore), decreasing = F))))
# colnames(clustmd)[1] <- "Id"
# clustmd <- merge(clustmd, cluster_time_score, by="Clusters")
#
# if (!is.null(cluster_labels)){
#   clustmd$label <- paste0("Cluster ", paste(rownames(testclustmd), cluster_labels, sep="-"))
# } else {
#   clustmd$label <- paste("Cluster ", rownames(testclustmd))
# }
