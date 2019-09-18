#' Processed time-series scRNAseq data of murine cerebral cortex development
#'
#' The dataset contains approximately 6,000 neural cells collected at embryonic days 11.5 (E11.5), E13.5, E15.5 and E17.58. Cells were sequenced using DropSeq.
#' These cells cover a wide spectrum of neuronal development, from the early precursors (apical and radial precursors) to intermediate progenitors
#' and differentiated cortical neurons. Cells from all timepoints were filtered to remove low-quality reads, normalized with \code{"scran"} and corrected for batch effect using \code{"Harmony"}.
#' The data was also filtered to remove non-cortical cells, as done in the original publication. Removed cells included cells expressing Aif1 (microglia), hemoglobin genes (blood cells), collagen genes (mesenchymal cells), as well as Dlx transcription factors and/or
#' interneuron genes (ganglionic eminence-derived cells). All retained cells were then iteratively clustered.
#' Raw sequencing reads can be accessed in the Gene Expression Omnibus, accession number GSE107122.
#'
#' @docType data
#'
#' @usage data(MouseCortex)
#'
#' @format A \code{"Seurat"} object containing a processed expression matrix of all cortical cells from four time-points (E11.5, E13.5, E15.5. E17.5), clustering results at multiple resolutions.
#'
#' @keywords mouse cortex, dataset
#'
#' @references Yuzwa, Scott A., et al. "Developmental emergence of adult neural stem cells as revealed by single-cell transcriptional profiling." Cell reports 21.13 (2017): 3970-3986. doi:https://doi.org/10.1016/j.celrep.2017.12.017.
#'
#' @examples
#' data(MouseCortex)
#' mousecortex_tempora <- ImportSeuratObject(MouseCortex, clusters = "res.0.6", timepoints = "Time_points",
#' cluster_labels = c("Neurons","Young neurons","APs/RPs","IPs","APs/RPs", "Young neurons", "IPs"), timepoint_order = c("e11", "e13", "e15", "e17"))
