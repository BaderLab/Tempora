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
    #eventhough the function was not successful return the object because in the vignette
    # this function call sets the original object to what is returned and if it is null
    # you loose all the processing you have done until now.
    return(object)
  } else {
    object@varying.pws <- varying_pathways
    object@gams <- gams
    return(object)
  }
}


#  ----
#' Calculate temporally changing pathways (parallel version)
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
#' @examples \dontrun{tempora_data <- IdentifyVaryingPWsParallel(tempora_data, pval_threshold = 0.05)}
IdentifyVaryingPWsParallel <- function(object, pval_threshold=0.05){
  require(BiocParallel)
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

  if (DEBUG) cat("Fitting GAM models...")

  # system.time({
  p_vals <- gams <- list()

  func = function(idx, object, themes, gsva_bycluster) {
    require(Tempora)
    require(Matrix)
    require(BiocGenerics)
    # cat(file = stderr(), paste(idx, "\n"))
    grp = BiocGenerics::grep(themes[idx], rownames(gsva_bycluster), ignore.case = T)
    crit = length(grp)
    if(crit == 0) {
      return(list(1, NA))
    }
    if (crit > 1){
      plot_df <- data.frame(
        cluster=colnames(gsva_bycluster[grp, ]),
        value=Matrix::colMeans(gsva_bycluster[grp, ],
                               na.rm=T))
    } else if (crit == 1){
      plot_df <- data.frame(cluster=names(gsva_bycluster[grp, ]),
                            value=gsva_bycluster[grp, ])
    } else {
      #should not happen
      return(list(1, NA))
    }
    plot_df$time <- object@cluster.metadata$Cluster_time_score
    gams <- mgcv::gam(value ~ s(time, k=3, bs='cr'), data=plot_df)
    temp_anova <- mgcv::anova.gam(gams)
    p_vals = temp_anova$s.pv
    list(p_vals, gams)
  }

  p_valsNew <- bplapply(1:length(themes), function(x) { func(x, object, themes, gsva_bycluster) })
  for (idx in 1:length(p_valsNew)){
    p_vals[[idx]] = p_valsNew[[idx]][[1]]
    gams[[idx]] = p_valsNew[[idx]][[2]]
  }
  names(p_vals) <- names(gams) <- themes

  pval_threshold = pval_threshold
  p_vals_adj <- p.adjust(unlist(p_vals[which(unlist(p_vals) > 0)]), method = "BH")
  varying_pathways <- p_vals_adj[which(p_vals_adj < pval_threshold)]
  varying_pathways <- varying_pathways[!duplicated(names(varying_pathways))]

  if (length(varying_pathways)==0){
    cat("No temporally varying pathways detected. Please try running IdentifyVaryingPWs with a more relaxed p-value cutoff.")
    #eventhough the function was not successful return the object because in the vignette
    # this function call sets the original object to what is returned and if it is null
    # you loose all the processing you have done until now.
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("No temporally varying pathways detected.", id = "temporaIdentifyVaryingPWs2", duration = 20)
    }
    return(object)
  } else {
    object@varying.pws <- varying_pathways
    object@gams <- gams
    return(object)
  }
  # to = object
}

