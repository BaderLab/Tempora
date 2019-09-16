#' Calculate pathway enrichment profile
#'
#' Calculate cluster average gene expression profile and determine the pathway enrichment profile of each cluster
#' @param object A Tempora object containing a gene expression matrix and metadata (cluster identity and )
#' @param pathwaygmt A database of pathways oirganized as a .gmt file
#' @param method Method used to calculate pathway enrichment profile. Can be "gsva", "ssgsea", "zscore" or "plage". See ?gsva for more information.
#' @param min.sz
#' @param max.sz
#' @param parallel.sz
#' @export
#' CalculatePWProfiles
#'
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


  devtools::use_package("mgcv")

  gsva_bycluster <- object@cluster.pathways

  significant_pathways <- c()
  for (i in 1:object@n.pcs){
    genes_scaled <- scale(object@cluster.pathways.dr$rotation[,i])
    significant_pathways <- c(names(which(genes_scaled[,1] > 1.5 | genes_scaled[,1] < -1.5)), significant_pathways)
  }

  pca_pathways <- sub("%.*", "", significant_pathways)
  pca_pathways_cleaned <- gsub("\\.|_|\\+", " ", pca_pathways)
  themes <- pca_pathways_cleaned

  cat("Fitting GAM models...")

  p_vals <- gams <- list()
  for (i in 1:length(themes)){
    if (length(grep(themes[i], rownames(gsva_bycluster))) > 1){
      plot_df <- data.frame(cluster=colnames(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), value=colMeans(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]))
    } else if  (length(grep(themes[i], rownames(gsva_bycluster))) == 1){
      plot_df <- data.frame(cluster=names(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), value=gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]) }
    plot_df$time <- object@cluster.metadata$Cluster_time_score
    gams[[i]] <- mgcv::gam(value ~ s(time, k=3, bs='cr'), data=plot_df)
    temp_anova <- mgcv::anova.gam(gams[[i]])
    p_vals[[i]] <- temp_anova$s.pv
  }

  names(p_vals) <- names(gams) <- themes

  p_vals_adj <- p.adjust(unlist(p_vals[which(unlist(p_vals) > 0)]), method = "BH")
  varying_pathways <- p_vals_adj[which(p_vals_adj < pval_threshold)]
  varying_pathways <- varying_pathways[!duplicated(varying_pathways)]

  cat("\nPlotting time-dependent pathways...")

  for (i in 1:length(varying_pathways)){
    if (length(grep(names(varying_pathways)[i], rownames(gsva_bycluster))) > 1){
      plot_df <- data.frame(cluster=colnames(gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ]), value=colMeans(gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ]))
      plot_df$time <- object@cluster.metadata$Cluster_time_score
    }
    else if (length(grep(names(varying_pathways)[i], rownames(gsva_bycluster))) == 1) {
      plot_df <- data.frame(cluster=names(gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ]), value=gsva_bycluster[grep(names(varying_pathways)[i], rownames(gsva_bycluster)), ])
      plot_df$time <- object@cluster.metadata$Cluster_time_score
    }
    id <- which(names(gams)==names(varying_pathways)[i])
    mgcv::plot.gam(gams[[id[1]]], main = paste0(names(varying_pathways)[i]), xlab = "Inferred time", ylab="Pathway expression level", bty="l",
             cex.main = 1, xaxt = "n", shade= F, se=3, scheme=1)
    xmin <- par("usr")[1]
    xmax <- par("usr")[2]
    points(x=plot_df$time, y=plot_df$value, pch=20, col="navy", cex=0.9)
    text(x=plot_df$time, y=plot_df$value, labels=plot_df[,1], pos = 4, cex = 1, col="navy")
    legend("topright", legend = "Cluster", pch = 20, col = "navy", bty="n", text.col="navy", cex=0.9)
    legend("topright", legend=paste0("\nAdjusted p-value = ", round(varying_pathways[[i]], 5)), bty="n", cex=0.9)
    axis(side=1, at=c(xmin, xmax), labels = c("Early", "Late"), tick=T)
    }
}

# IdentifyVaryingPWs(testobj, pval_threshold = 0.05)
