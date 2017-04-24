#' PCA of samples (if use of DESeq2)
#'
#' Principal Component Analysis of samples based on the 500 most variant features on VST- or rlog-counts (if use of DESeq2)
#'
#' @param counts.trans a matrix a transformed counts (VST- or rlog-counts)
#' @param group factor vector of the condition from which each sample belongs
#' @param n number of features to keep among the most variant
#' @param col colors to use (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named PCA.png in the figures directory with a pairwise plot of the three first principal components
#' @author Marie-Agnes Dillies and Hugo Varet

library(ggbiplot)
library(gridExtra)

PCAPlot <- function(counts.trans, group, n=min(500,nrow(counts.trans)), 
                    col=c("lightblue","orange","MediumVioletRed","SpringGreen"),
                    outfile=TRUE){
  # PCA on the 500 most variables features
  rv = apply(counts.trans, 1, var, na.rm=TRUE)
  pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE), ][1:n,]))
  prp <- pca$sdev^2 * 100 / sum(pca$sdev^2)
  prp <- round(prp[1:3],2)

  # create figure
  if (outfile) png(filename="figures/PCA.png",width=1800*2,height=1800,res=300)
    par(mfrow=c(1,2))
	# axes 1 et 2
  p1 <- ggbiplot::ggbiplot(pca, obs.scale = 1, var.scale = 1, 
                    groups = group, ellipse = TRUE, 
                    circle = FALSE, var.axes = FALSE) + 
                    fte_theme() +
                    theme(aspect.ratio=1) +
                    scale_colour_brewer(type = "qual", palette = 6)

	# axes 1 et 3
  p2 <- ggbiplot::ggbiplot(pca, choices = c(1,3), obs.scale = 1, var.scale = 1, 
                    groups = group, ellipse = TRUE, 
                    circle = FALSE, var.axes = FALSE) +
                    fte_theme() +
                    theme(aspect.ratio=1) +
                    scale_colour_brewer(type = "qual", palette = 6) +
                    theme(legend.position=c(.9,.9))
  
  grid.arrange(p1, p2, ncol=2)
  if (outfile) dev.off()

  return(invisible(pca$x))
}
