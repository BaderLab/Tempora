#' Plot temporally changing pathways
#'
#' Plot the expression of temporally changing pathways as identified by IdentifyVaryingPWs()
#' @param object A Tempora object
#' @export
#' @importFrom mgcv gam anova.gam plot.gam
#' @importFrom methods new validObject
#' @importFrom stats p.adjust prcomp screeplot
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis legend par points text
#' @importFrom reshape2 dcast
#' @examples \dontrun{tempora_data <- IdentifyVaryingPWs(tempora_data, pval_threshold = 0.05)}
PlotVaryingPWs <- function(object){

  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }
  if (is.null(object@varying.pws)){
    stop("IdentifyVaryingPWs has not been run or no temporally varying pathways were detected. Please run IdentifyVaryingPWs or re-run with a more relaxed p-value cutoff See ?Tempora::IdentifyVaryingPWs for details")
  }

varying_pathways <- object@varying.pws
gsva_bycluster <- object@cluster.pathways
gams <- object@gams

cat("\nPlotting time-dependent pathways...")

for (i in order(varying_pathways)){
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
