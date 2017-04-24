#' Total number of reads per sample
#'
#' Bar plot of the total number of reads per sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one color per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named barplotTotal.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

barplotTotal <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  if (outfile) png(filename="figures/barplotTotal.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
  libsize <- data.frame(sample=colnames(counts), group=group, mil_reads=colSums(counts)/1e6)
  libsize$sample <- factor(libsize$sample,levels=unique(libsize$sample))
  print(ggplot(libsize, aes(x=sample, y=mil_reads, fill=group)) +
          geom_bar(stat='identity', position='dodge') +
          scale_y_continuous(limits=c(0,max(libsize$mil_reads)*1.1)) +
          ylab("Total read count (million)") +
          xlab ("") +
          ggtitle("Total read count per sample (million)") +
          fte_theme() +
          scale_fill_brewer(type = "qual", palette = 6) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size=round(400/nrow(libsize)))) +
          theme(legend.position=c(.9,.9))
  )    
  if (outfile) dev.off()
}
