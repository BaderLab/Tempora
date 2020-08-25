#' Calculate temporally changing pathways
#'
#' Identify the pathways that change over time by fitting a generalized additive model
#' @param object A Tempora object
#' @param pval_threshold P-value threshold to determine the significance of pathway enrichment over time. Default to 0.05.
#' @export
#' @importFrom mgcv gam anova.gam plot.gam
#' @importFrom methods new validObject
#' @importFrom stats p.adjust prcomp screeplot
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis legend par points text
#' @importFrom reshape2 dcast
#' @examples \dontrun{tempora_data <- IdentifyVaryingPWs(tempora_data, pval_threshold = 0.05)}

IdentifyVaryingPWs <- function(object, pval_threshold=0.05){

  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }
  if (is.null(object@n.pcs)){
    stop("BuildTrajectory has not been run. See ?Tempora::BuildTrajectory for details")
  }
  if (is.null(object@cluster.pathways)){
    stop("CalculatePWProfiles has not been run. See ?Tempora::CalculatePWProfiles for details")
  }
  gsva_bycluster <- object@cluster.pathways

  significant_pathways <- c()
  for (i in 1:object@n.pcs){
    genes_scaled <- scale(object@cluster.pathways.dr$rotation[,i])
    significant_pathways <- c(names(which(genes_scaled[,1] > 1.5 | genes_scaled[,1] < -1.5)), significant_pathways)
  }

  pca_pathways <- sub("%.*", "", significant_pathways)
  pca_pathways <- gsub("\\s*\\([^\\)]+\\)","",pca_pathways)
  pca_pathways_cleaned <- gsub("[[:punct:]]", "", pca_pathways)
  themes <- pca_pathways_cleaned

  cat("Fitting GAM models...")

  p_vals <- gams <- list()
  for (i in 1:length(themes)){
    print(i)
    if(length(grep(themes[i], rownames(gsva_bycluster))) == 0) {
      p_vals[[i]] <- 1
      gams[[i]] <- NA
      next
      }
    if (length(grep(themes[i], rownames(gsva_bycluster))) > 1){
      plot_df <- data.frame(cluster=colnames(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), value=colMeans(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ], na.rm=T))
    } else if (length(grep(themes[i], rownames(gsva_bycluster))) == 1){
      plot_df <- data.frame(cluster=names(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), value=gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]) }
    plot_df$time <- object@cluster.metadata$Cluster_time_score
    gams[[i]] <- mgcv::gam(value ~ s(time, k=3, bs='cr'), data=plot_df)
    temp_anova <- mgcv::anova.gam(gams[[i]])
    p_vals[[i]] <- temp_anova$s.pv
  }

  names(p_vals) <- names(gams) <- themes

  pval_threshold = pval_threshold
  p_vals_adj <- p.adjust(unlist(p_vals[which(unlist(p_vals) > 0)]), method = "BH")
  varying_pathways <- p_vals_adj[which(p_vals_adj < pval_threshold)]
  varying_pathways <- varying_pathways[!duplicated(names(varying_pathways))]

  if (length(varying_pathways)==0){
    cat("No temporally varying pathways detected. Please try running IdentifyVaryingPWs with a more relaxed p-value cutoff.")
  } else {
    object@varying.pws <- varying_pathways
    object@gams <- gams
    return(object)
  }
}


