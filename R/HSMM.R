#' Processed time-series scRNAseq data of \textit{"in vitro"} differentiation of human skeletal muscle myoblasts
#'
#' The dataset contains approximately 271 cells collected at 0, 24, 48 and 72 hours after the switch of human myoblast culture from growth to differentiation media.
#' Cells were sequenced using Fluidigm C1.
#' Cells from all timepoints were filtered to remove low-quality reads, normalized with \code{"scran"} and corrected for batch effect using \code{"Harmony"}.
#' The data was also filtered to remove non-cortical cells, as done in the original publication. All retained cells were then iteratively clustered.
#' Raw sequencing reads can be accessed in the Gene Expression Omnibus, accession number GSE52529.
#'
#' @docType data
#'
#' @usage data(HSMM)
#'
#' @format A \code{"Seurat"} object containing a processed expression matrix of all cells from four time-points (0, 24, 48 and 72 hours after the switch from growth to differentiation media), clustering results at multiple resolutions.
#'
#' @keywords HSMM, dataset
#'
#' @references Trapnell, Cole, et al. "The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells." Nature biotechnology 32.4 (2014): 381.
#'
#' @examples
#' data(HSMM)
#' hsmm_tempora <- ImportSeuratObject(HSMM, clusters = "res.1.5", timepoints = "Time_points", cluster_labels = c("Myoblasts", "Fibroblasts", "Myotubes","Undifferentiated","Intermediates"), timepoint_order = c("T0", "T24", "T48", "T72"))
