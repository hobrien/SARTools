#' Density plot of all samples
#'
#' Estimation the counts density for each sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the curves (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named densplot.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

densityPlot <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  if (outfile) png(filename="figures/densplot.png",width=1800,height=1800,res=300)
    labels <- data.frame(Sample=paste0('X', colnames(counts)), group=group)
    counts2 <- removeNull(counts)
    counts2 <- data.frame(counts2)
    counts2$Id<- rownames(counts2)
    counts2 <- gather(counts2, Sample, Count, -Id) %>% full_join(labels)
    if (is.numeric(group)) {
      palette <- 15
      type <- "seq"
    } else {
      palette <- 6
      type = "qual"
    }
    print(ggplot(counts2, aes(x=log2(Count+1), group=Sample, colour=factor(group))) + 
            geom_density() +
            ggtitle("Density of counts distribution") +
            fte_theme() +
            scale_colour_brewer(type=type, palette=palette) +
            theme(legend.position=c(.9,.9))
        )
  if (outfile) dev.off()
}
