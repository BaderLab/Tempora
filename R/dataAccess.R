# Generics & methods for loading data from various single-cell data classes.

# Generics ----

# ^ getExpr ----

#' Get gene expression matrix from input data object
#'
#' Extract the gene expression matrix from a single-cell data object containing
#' the input data for scClustViz visualization.
#'
#' This is a wrapper function to the relevant class's normalized data slot
#' accessor method. Currently supported input object classes: \itemize{ \item
#' Class \code{\link[Seurat]{seurat}/\link[Seurat]{Seurat}} stored in
#' \code{x@data} or \code{x@assays[[assayType]]@assaySlot}, depending on Seurat
#' object version. \item Class
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} accessed by
#' \code{\link[SummarizedExperiment]{assay}(x,assayType)}. }
#' \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#' for other data objects here!}
#'
#' @param x The single-cell data object.
#' @param assayType A length-one character vector representing the assay object
#'   in which the expression data is stored in the input object. For Seurat v1
#'   or v2 objects, set this to "". For Seurat v3 objects, this is often "RNA".
#'   For SingleCellExperiment objects, this is often "logcounts". See Details
#'   for how this argument is used in the accessor functions for each class.
#' @param assaySlot An optional length-one character vector representing the
#'   slot of the Seurat v3 \code{\link[Seurat]{Assay}} object to use. In Seurat
#'   v3, normalized data is stored in the "data" slot, and counts in the
#'   "counts" slot. See Details for how this argument is used in the accessor
#'   functions for each class.
#' @name getExpr
#' @author Contributed by Brendan Innes from the BaderLab/scClustViz package
#' @export
#'
setGeneric("getExpr",function(x,assayType,assaySlot) standardGeneric("getExpr"))


# ^ getMD ----

#' Get metadata from input data object
#'
#' Extract the cell metadata data frame from a single-cell data object
#' containing the input data for scClustViz visualization.
#'
#' This is a wrapper function to the relevant class's cell metadata slot
#' accessor / assignment method. Currently supported input object classes:
#' \itemize{
#'   \item Class \code{\link[Seurat]{seurat}/\link[Seurat]{Seurat}} accessed by
#'     \code{x@data.info} or \code{x@meta.data},
#'     depending on Seurat object version.
#'   \item Class \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'     accessed by \code{\link[SingleCellExperiment]{colData}(x)}.
#' }
#' \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#' for other data objects here!}
#'
#' @param x The single-cell data object.
#' @name getMD
#' @author Contributed by Brendan Innes from the BaderLab/scClustViz package
#' @export
#'
setGeneric("getMD",function(x) standardGeneric("getMD"))


# Methods ----

# ^ seurat (v1/2) ----
suppressMessages(
  setMethod("getExpr","seurat",function(x) {
    slot(x,"data")
  })
)


suppressMessages(
  setMethod("getMD","seurat",function(x) {
    if (.hasSlot(x,"meta.data")) {
      slot(x,"meta.data")
    } else {
      slot(x,"data.info") #Seurat v1
    }
  })
)


# ^ Seurat (v3) ----
suppressMessages(
  setMethod("getExpr","Seurat",function(x,assayType,assaySlot) {
    if (missing(assayType)) {
      stop(paste(paste0("assayType must be specified."),
                 "The following assay data are available in this object:",
                 paste0(names(slot(x,"assays")),collapse=", "),sep="\n  "))
    }
    if (assayType %in% names(slot(x,"assays"))) {
      if (!missing(assaySlot)) {
        if (is.na(assaySlot) | assaySlot == "") {
          return(x@assays[[assayType]]@data)
        } else {
          return(slot(x@assays[[assayType]],assaySlot))
        }
      } else {
        return(x@assays[[assayType]]@data)
      }
    } else {
      stop(paste(paste0("assayType '",assayType,"' not found."),
                 "The following assay data are available in this object:",
                 paste0(names(slot(x,"assays")),collapse=", "),sep="\n  "))
    }
  })
)


suppressMessages(
  setMethod("getMD","Seurat",function(x) {
    return(slot(x,"meta.data"))
  })
)


# ^ SingleCellExperiment ----
suppressMessages(
  setMethod("getExpr","SingleCellExperiment",function(x,assayType) {
    if (missing(assayType)) {
      stop(paste(paste0("assayType must be specified."),
                 "The following assay data are available in this object:",
                 paste0(SummarizedExperiment::assayNames(x),collapse=", "),
                 sep="\n  "))
    }
    if (assayType %in% SummarizedExperiment::assayNames(x)) {
      return(SummarizedExperiment::assay(x,assayType))
    } else {
      stop(paste(paste0("assayType '",assayType,"' not found."),
                 "The following assay data are available in this object:",
                 paste0(SummarizedExperiment::assayNames(x),collapse=", "),
                 sep="\n  "))
    }
  })
)


suppressMessages(
  setMethod("getMD","SingleCellExperiment",
            function(x) data.frame(SingleCellExperiment::colData(x)))
)
