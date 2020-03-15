#' Calculate pathway enrichment profile
#'
#' Calculate cluster average gene expression profile and determine the pathway enrichment profile of each cluster
#' @param object A Tempora object containing a gene expression matrix and metadata (cluster identity and )
#' @param gmt_path Local path to database of pathways or genesets organized as a .gmt file. Genesets files in GMT format can be downloaded at http://baderlab.org/GeneSets. Please ensure
#' @param method Method used to estimate pathway enrichment profile per cluster. Can be "gsva", "ssgsea", "zscore" or "plage", default to "gsva". See ?gsva for more information.
#' @param min.sz Minimum size of the genesets used in enrichment estimation, set to 5 genes by default.
#' @param max.sz Maximum size of the genesets used in enrichment estimation, set to 200 genes by default.
#' @param parallel.sz Type of cluster architecture when using \code{snow}. If 1, no parallelization will be used. If 0, all available cores will be used.
#' @export
#' @importFrom methods new validObject
#' @importFrom stats p.adjust prcomp screeplot
#' @importFrom GSVA gsva
#' @importFrom GSEABase getGmt
#' @importFrom reshape2 dcast
#' @examples \dontrun{tempora_data <- CalculatePWProfiles(tempora_data, gmt_path="~/Human_AllPathways_September_01_2019_symbol.gmt", parallel.sz = detectCores()-2)}
#' @return An updated Tempora object containing the pathway enrichment profiles of each cluster, which can be accessed at \code{object@cluster.pathways}
#' CalculatePWProfiles
CalculatePWProfiles <- function(object, gmt_path, method="gsva", min.sz=5, max.sz=200, parallel.sz=1){
  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }

  cat("Calculating cluster average gene expression profile...")
  exprMatrix <- object@data
  exprMatrix_bycluster <- list()
  pathwaygmt <- GSEABase::getGmt(gmt_path)
  for (i in sort(unique(object@meta.data$Clusters))){
    exprMatrix_bycluster[[i]] <- rowMeans(exprMatrix[, which(colnames(exprMatrix) %in% rownames(object@meta.data)[which(object@meta.data$Clusters == i)])])
  }

  exprMatrix_bycluster <- do.call(cbind, exprMatrix_bycluster)
  colnames(exprMatrix_bycluster) <- sort(unique(object@meta.data$Clusters))
  rownames(exprMatrix_bycluster) <- rownames(exprMatrix)

  cat("\nCalculating cluster pathway enrichment profiles...\n")
  gsva_bycluster <- GSVA::gsva(as.matrix(exprMatrix_bycluster), pathwaygmt, method=method, min.sz=min.sz, max.sz=max.sz, parallel.sz=parallel.sz)

  colnames(gsva_bycluster) <- colnames(exprMatrix_bycluster)
  object@cluster.pathways <- gsva_bycluster

  gsva_bycluster_pca <- prcomp(t(gsva_bycluster), scale = T, center = T)
  screeplot(gsva_bycluster_pca, npcs=25, type="lines", main="PCA on pathway enrichment analysis result")
  object@cluster.pathways.dr <- gsva_bycluster_pca

  validObject(object)
  return(object)
}


















