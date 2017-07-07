#' Clustering of the samples
#'
#' Clustering of the samples based on VST- or rlog-counts (if use of DESeq2) or cpm-counts (if use of edgeR)
#'
#' @param counts.trans a matrix a transformed counts (VST- or rlog-counts if use of DESeq2 or cpm-counts if use of edgeR)
#' @param group factor vector of the condition from which each sample belongs
#' @param outfile TRUE to export the figure in a png file
#' @return A file named cluster.png in the figures directory with the dendrogram of the clustering
#' @author Marie-Agnes Dillies and Hugo Varet

library(ggdendro)

clusterPlot <- function(counts.trans, group, outfile=TRUE){
  hc <- hclust(dist(t(counts.trans)), method="ward.D")
  if (is.numeric(group)) {
    palette <- 15
    type <- "seq"
  } else {
    palette <- 6
    type = "qual"
  }
  if (outfile) png(filename="figures/cluster.png",width=1800,height=1800,res=300) 
  ddata <- dendro_data(as.dendrogram(hc), type = "rectangle")
  ddata$labels <- full_join(ddata$labels, data.frame(label=colnames(counts.trans), group=group))
  print(ggplot(segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = ddata$labels, 
              aes(x = x, y = y, label = label, colour=factor(group)), angle = 90, hjust = 1, size=round(150/length(group))) +
    ylab("Height") +
    xlab ("Method: Euclidean distance - Ward criterion") +
    ggtitle("Cluster dendrogram") +
    fte_theme() +
    scale_colour_brewer(type = type, palette = palette) +
    expand_limits(y=-1500/length(group)) +
    #scale_y_continuous(limits=c(-100, 1500)) +
    theme(axis.text.x = element_blank()) +
    theme(legend.position=c(.9,.9))
  )
  if (outfile) dev.off()
}
