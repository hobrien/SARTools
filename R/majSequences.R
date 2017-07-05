#' Most expressed sequences per sample
#'
#' Proportion of reads associated with the three most expressed sequences per sample
#'
#' @param counts \code{matrix} of counts
#' @param n number of most expressed sequences to return
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A \code{matrix} with the percentage of reads of the three most expressed sequences and a file named majSeq.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

majSequences <- function(counts, n=3, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){

  seqnames <- apply(counts, 2, function(x){x <- sort(x, decreasing=TRUE); names(x)[1:n]})
  seqnames <- unique(unlist(as.character(seqnames)))

  sum <- apply(counts,2,sum)
  counts2 <- counts[seqnames,]
  sum <- matrix(sum,nrow(counts2),ncol(counts2),byrow=TRUE)
  p <- round(100*counts2/sum,digits=3)
  maj <- data.frame(apply(p, 2, max))
  seqname <- rownames(p)[apply(p, 2, which.max)]
  maj_seq <- data.frame(Sample=rownames(maj), Id=seqname, prop=maj[,1], group=group)
  if (is.numeric(group)) {
    palette <- 15
    type <- "seq"
  } else {
    palette <- 6
    type = "qual"
  }
  if (outfile) png(filename="figures/majSeq.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
  print(ggplot(maj_seq, aes(x=Sample, y=prop, fill=factor(group))) + 
    geom_bar(stat='identity') +
    geom_text(aes(label=Id), angle=90, size=2) +
    scale_fill_brewer(type=type, palette=palette) +
    ggtitle("Proportion of reads from most expressed sequence") +
    ylab("Proportion of reads") +
    fte_theme() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=round(400/length(group)))) +
      theme(legend.position=c(.9,.9))
  )
  if (outfile) dev.off()
  
  return(invisible(p))
}
