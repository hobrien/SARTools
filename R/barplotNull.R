#' Percentage of null counts per sample
#'
#' Bar plot of the percentage of null counts per sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one color per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named barplotNull.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

barplotNull <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  if (outfile) png(filename="figures/barplotNull.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
    percentage <- apply(counts, 2, function(x){sum(x == 0)})*100/nrow(counts)
    percentage.allNull <- (nrow(counts) - nrow(removeNull(counts)))*100/nrow(counts)
    nullPecr <- data.frame(percentage=percentage, sample=colnames(counts))
    nullPecr$sample <- factor(nullPecr$sample,levels=unique(nullPecr$sample))
    if (is.numeric(group)) {
      palette <- 15
      type <- "seq"
    } else {
      palette <- 6
      type = "qual"
    }
    print(ggplot(nullPecr, aes(x=sample, y=percentage, fill=factor(group))) +
            geom_bar(stat='identity', position='dodge') +
            geom_abline(intercept = percentage.allNull, slope=0, linetype=2) +
            scale_y_continuous(limits=c(0,max(nullPecr$percent)+10)) +
            ylab("Proportion of null counts") +
            xlab ("") +
            ggtitle("Proportion of null counts per sample") +
            fte_theme() +
            scale_fill_brewer(type=type, palette=palette) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, size=round(400/nrow(nullPecr)))) +
            theme(legend.position=c(.9,.9))
          
    )    
  if (outfile) dev.off()
}
