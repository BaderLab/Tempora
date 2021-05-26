
###define S4 class
#' Define a class of Tempora object
#'
#' A Tempora object contains the input gene expression matrix and metadata, as well as stores the meta data of each cluster,
#' the clusters' pathway enrichment profiles, the constructed trajectory as well as the Sugiyama layout for the trajectory plot
#'
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportClass Tempora
#'
#' @slot data A gene expression matrix (genes x cells), often aggregated from multiple time points
#' @slot meta.data A dataframe containing the metadata for the cells in the gene expression matrix, which at minimum includes the
#' collection timepoint and cluster identity of each cell
#' @slot cluster.metadata A dataframe containing the metadata for each cell cluster
#' @slot cluster.pathways A dataframe containing the pathway enrichment profile of each cluster as calculated by \code{\link{CalculatePWProfiles}}
#' @slot cluster.pathways.dr A prcomp object containing the PCA of the clusters' pathway enrichment profiles
#' @slot n.pcs The number of principal components to be used in trajectory construction
#' @slot trajectory A dataframe describing the inferred trajectory as inferred by \code{\link{BuildTrajectory}}
#' @slot layouts A matrix containing the Sugiyama layout of the trajectory to be used in \code{\link{PlotTrajectory}}
#' @slot varying.pws A list of temporally varying pathways identified by \code{\link{IdentifyVaryingPWs}}
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
    layouts = "matrix",
    varying.pws = "ANY",
    gams = "ANY"
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

#' Data method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod data
setGeneric("data", function(x) standardGeneric("data"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod data<-
setGeneric("data<-", function(x, value) standardGeneric("data<-"))
#' Extract data from Tempora object
#'
#' @rdname Tempora-class
#' @aliases data
#' @param x Tempora object
setMethod("data", "Tempora", function(x) x@data)

#' @rdname Tempora-class
#' @aliases data<-
#' @param value New value
setMethod("data<-", "Tempora", function(x, value) {
  x@data <- value
  validObject(x)
  return(x)
})

#' Metadata method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod meta.data
setGeneric("meta.data", function(x) standardGeneric("meta.data"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod meta.data<-
setGeneric("meta.data<-", function(x, value) standardGeneric("meta.data<-"))
#' Extract metadata from Tempora object
#'
#' @rdname Tempora-class
#' @aliases meta.data
setMethod("meta.data", "Tempora", function(x) x@meta.data)
#' Extract metadata from Tempora object
#'
#' @rdname Tempora-class
#' @aliases meta.data<-
setMethod("meta.data<-", "Tempora", function(x, value) {
  x@meta.data <- value
  validObject(x)
  return(x)
})

#' Cluster metadata method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod cluster.metadata
setGeneric("cluster.metadata", function(x) standardGeneric("cluster.metadata"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod cluster.metadata<-
setGeneric("cluster.metadata<-", function(x, value) standardGeneric("cluster.metadata<-"))
#' Extract cluster metadata from Tempora object
#'
#' @rdname Tempora-class
#' @aliases cluster.metadata
setMethod("cluster.metadata", "Tempora", function(x) x@cluster.metadata)
#' @rdname Tempora-class
#' @aliases cluster.metadata<-
setMethod("cluster.metadata<-", "Tempora", function(x, value) {
  x@cluster.metadata <- data.frame
  validObject(x)
  return(x)
})


#' Cluster pathway enrichment profile method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod cluster.pathways
setGeneric("cluster.pathways", function(x) standardGeneric("cluster.pathways"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod cluster.pathways<-
setGeneric("cluster.pathways<-", function(x, value) standardGeneric("cluster.pathways<-"))
#' Extract cluster pathway enrichment profiles from Tempora object
#'
#' @rdname Tempora-class
#' @aliases cluster.pathways
setMethod("cluster.pathways", "Tempora", function(x) x@cluster.pathways)
#' @rdname Tempora-class
#' @aliases cluster.pathways<-
setMethod("cluster.pathways<-", "Tempora", function(x, value) {
  x@cluster.pathway <- value
  validObject(x)
  return(x)
})

#' Dimension reduction of cluster pathway enrichment profile method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod cluster.pathways.dr
setGeneric("cluster.pathways.dr", function(x) standardGeneric("cluster.pathways.dr"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod cluster.pathways.dr<-
setGeneric("cluster.pathways.dr<-", function(x, value) standardGeneric("cluster.pathways.dr<-"))
#' Extract PCA of cluster pathway enrichment profiles from Tempora object
#'
#' @rdname Tempora-class
#' @aliases cluster.pathways.dr
setMethod("cluster.pathways.dr", "Tempora", function(x) x@cluster.pathways.dr)
#' @rdname Tempora-class
#' @aliases cluster.pathways.dr<-
setMethod("cluster.pathways.dr<-", "Tempora", function(x, value) {
  x@cluster.pathways.dr <- value
  validObject(x)
  return(x)
})

#' Varying pathways method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod varying.pws
setGeneric("varying.pws", function(x) standardGeneric("varying.pws"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod varying.pws<-
setGeneric("varying.pws<-", function(x, value) standardGeneric("varying.pws<-"))
#' Extract varying pathways from Tempora object
#'
#' @rdname Tempora-class
#' @aliases varying.pws
setMethod("varying.pws", "Tempora", function(x) x@varying.pws)
#' @rdname Tempora-class
#' @aliases varying.pws<-
setMethod("varying.pws<-", "Tempora", function(x, value) {
  x@varying.pws <- value
  validObject(x)
  return(x)
})


#' GAMs method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod gams
setGeneric("gams", function(x) standardGeneric("gams"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod gams<-
setGeneric("gams<-", function(x, value) standardGeneric("gams<-"))
#' Extract calculated GAMs from Tempora object
#'
#' @rdname Tempora-class
#' @aliases gams
setMethod("gams", "Tempora", function(x) x@gams)
#' @rdname Tempora-class
#' @aliases gams<-
setMethod("gams<-", "Tempora", function(x, value) {
  x@gams <- value
  validObject(x)
  return(x)
})


#' Number of PCs to use in trajectory construction method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod n.pcs
setGeneric("n.pcs", function(x) standardGeneric("n.pcs"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod n.pcs<-
setGeneric("n.pcs<-", function(x, value) standardGeneric("n.pcs<-"))
#' Extract number of PCA of cluster pathway enrichment profiles to use from Tempora object
#'
#' @rdname Tempora-class
#' @aliases n.pcs
setMethod("n.pcs", "Tempora", function(x) x@cluster.pathways.dr)
#' @rdname Tempora-class
#' @aliases n.pcs<-
setMethod("n.pcs<-", "Tempora", function(x, value) {
  x@n.pcs <- value
  validObject(x)
  return(x)
})

#' Trajectory method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod trajectory
setGeneric("trajectory", function(x) standardGeneric("trajectory"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod trajectory
setGeneric("trajectory<-", function(x, value) standardGeneric("trajectory<-"))
#' Extract the constructed trajectory from Tempora object
#'
#' @rdname Tempora-class
#' @aliases trajectory
setMethod("trajectory", "Tempora", function(x) x@trajectory)
#' @rdname Tempora-class
#' @aliases trajectory<-
setMethod("trajectory<-", "Tempora", function(x, value) {
  x@trajectory <- value
  validObject(x)
  return(x)
})

#' Trajectory layout method
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod layouts
setGeneric("layouts", function(x) standardGeneric("layouts"))
#' @name Tempora-class
#' @rdname Tempora-class
#' @exportMethod layouts<-
setGeneric("layouts<-", function(x, value) standardGeneric("layouts<-"))
#' Extract the layout of the trajectory from Tempora object
#'
#' @rdname Tempora-class
#' @aliases layouts
setMethod("layouts", "Tempora", function(x) x@layouts)
#' @rdname Tempora-class
#' @aliases layouts<-
setMethod("layouts<-", "Tempora", function(x, value) {
  x@layouts <- value
  validObject(x)
  return(x)
})


###Create new Tempora object from expression matrix
#' Create new Tempora object from a processed gene expression matrix
#' @param exprMatrix A normalized gene expression matrix containing cells from all timepoints of the time-series study. Batch effect correction is highly recommended before normalization.
#' @param meta.data A dataframe of meta data for all cells in the expression matrix. Each column is a feature and each row stores single cell information. At minimum, this dataframe should contain two columns: a "Clusters" column storing the clustering identity and a "Timepoints" column storing the timepoint when each cell comes from
#' @param timepoint_order An ordered vector of timepoint names from early to late
#' @param cluster_labels (Optional) A vector of cluster annotations (cell types, cell states, cell cycles, etc.), ordered alphanumerically by cluster names. If NULL and \code{cell_markers} is given, automatic cluster annotation using GSVA will be performed. If both are NULL, cluster numbers will be used to label the trajectory plot.
#' @param cell_markers (Optional) A list of possible cell types found in the dataset and their marker genes to be used for automatic cell type identification. If NULL and no \code{cluster_labels} is given, cluster numbers will be used to label the trajectory plot.
#' @param threshold (Optional) Threshold of GSVA score quantile. Cell types with GSVA scores in this quantile or higher compared to all other cell type scores for the same cluster would be included in cluster label. Numeric between 0-1, default to 0.9.
#' @export
#' @importFrom methods new validObject
#' @importFrom stats p.adjust prcomp screeplot
#' @importFrom reshape2 dcast
#' @return A Tempora object containing the expression matrix and metadata
#' @examples \dontrun{tempora_dara <- CreateTemporaObject(exprMatrix, meta.data)}
#'
CreateTemporaObject <- function(exprMatrix, meta.data, timepoint_order, cluster_labels=NULL, cell_markers=NULL, threshold=0.9){

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
  if (!is.null(cluster_labels) & !is.null(cell_markers)){
    clustmd$label <- paste0("Cluster ", paste(rownames(clustmd), cluster_labels, sep="-"))
  } else if (!is.null(cluster_labels)){
    # clustmd$label <- paste0("Cluster ", paste(rownames(clustmd), cluster_labels, sep="-"))
    clustmd$label <-  cluster_labels
  } else if (!is.null(cell_markers)) {
    cat("\nPerforming automated cluster labeling with GSVA...")
    cluster_number <- as.numeric(meta.data$Clusters)
    names(cluster_number) <- rownames(meta.data)
    cluster_labels <- IdentifyCellTypes(exprMatrix, cluster_labels=cluster_number, cell_markers=cell_markers, threshold=threshold)
    clustmd$label <- paste0("Cluster ", paste(rownames(clustmd), cluster_labels, sep="-"))
  } else {
    clustmd$label <- paste("Cluster ", rownames(clustmd))
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
#' @param seuratobj A Seurat or SingleCellExperiment object containing the normalized gene expression matrix and clustering result
#' @param assayType A length-one character vector representing the assay object
#'   in which the expression data is stored in the input object. For Seurat v1
#'   or v2 objects, set this to "". For Seurat v3 objects, this is often "RNA".
#'   For SingleCellExperiment objects, this is often "logcounts".
#' @param assaySlot An optional length-one character vector representing the
#'   slot of the Seurat v3 \code{\link[Seurat]{Assay}} object to use. In Seurat
#'   v3, normalized data is stored in the "data" slot, and counts in the
#'   "counts" slot.
#' @param clusters Name of the column in the meta.data dataframe containing the cluster identity of all cells in the dataset
#' @param timepoints Name of the column in the meta.data dataframe containing the collection time of all cells in the dataset
#' @param timepoint_order An ordered vector of timepoint names from early to late
#' @param cluster_labels (Optional) A vector of cluster annotations (cell types, cell states, cell cycles, etc.), ordered alphanumerically by cluster names. If NULL and \code{cell_markers} is given, automatic cluster annotation using GSVA will be performed. If both are NULL, cluster numbers will be used to label the trajectory plot.
#' @param cell_markers (Optional) A list of possible cell types found in the dataset and their marker genes to be used for automatic cell type identification. If NULL and no \code{cluster_labels} is given, cluster numbers will be used to label the trajectory plot.
#' @param threshold (Optional) Threshold of GSVA score quantile. Cell types with GSVA scores in this quantile or higher compared to all other cell type scores for the same cluster would be included in cluster label. Numeric between 0-1, default to 0.9.
#' @include dataAccess.R

#' @export
#' @importFrom methods new validObject
#' @importFrom stats p.adjust prcomp screeplot
#' @importFrom reshape2 dcast
#' @return A Tempora object containing the expression matrix and metadata
#' @examples \dontrun{tempora_data <- ImportSeuratObject(seurat_object, assayType="", clusters = "res.0.3", timepoints = "collection_time",
#' timepoint_order = c("0H", "24H", "48H", "72H"), cluster_labels = c("Stem cells", "Differentiated cells"))}


ImportSeuratObject <- function(seuratobj, assayType = "", assaySlot = NA,
                               clusters, timepoints, timepoint_order, cluster_labels=NULL, cell_markers=NULL, threshold=0.9){
  # if(class(seuratobj)[1]=='seurat'){
  #   requireNamespace("Seurat")
  # } else {
  #   stop("Not a Seurat object. Tempora only supports importing Seurat objects at the moment. See ?Tempora::CreateTemporaObject to manually create a Tempora object from an expression matrix")
  # }
  data <- getExpr(seuratobj,assayType,assaySlot)
  if (! is.numeric(data)) {
    data <- as.matrix(data)
    # often Seurat / SingleCellExperiment data matrices are stored as sparse
    # Matrix::dgCMatrix objects. Since the S4 class requires this object to be
    # numeric, we must convert them to traditional numeric R matrices, despite
    # the increased memory costs this entails. Or allow sparse matrices in the
    # S4 class, but that would involve some further debugging...

  }
  cat("Extracting data...")
  metadata <- getMD(seuratobj)
  cat("\nExtracting metadata...")
  colnames(metadata)[which(colnames(metadata)==clusters)] <- "Clusters"
  colnames(metadata)[which(colnames(metadata)==timepoints)] <- "Timepoints"
  cat("\nCreating Tempora object...")
  tempora_obj <- CreateTemporaObject(data, meta.data = metadata, timepoint_order =  timepoint_order, cluster_labels = cluster_labels, cell_markers = cell_markers, threshold=threshold)
  validObject(tempora_obj)
  return(tempora_obj)
}



