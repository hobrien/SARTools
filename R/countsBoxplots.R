library(gridExtra)
#' Box-plots of (normalized) counts distribution per sample
#'
#' Box-plots of raw and normalized counts distributions per sample to assess the effect of the normalization
#'
#' @param object a \code{DESeqDataSet} object from DESeq2 or a \code{DGEList} object from edgeR
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the boxplots (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named countsBoxplots.png in the figures directory containing boxplots of the raw and normalized counts
#' @author Marie-Agnes Dillies and Hugo Varet

countsBoxplots <- function(object, group, col = c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  if (class(object)=="DESeqDataSet"){
    counts <- counts(object)
    counts <- removeNull(counts)
    norm.counts <- counts(object, normalized=TRUE)
    norm.counts <- removeNull(norm.counts)  
  } else{
    counts <- object$counts
    counts <- removeNull(counts)
	tmm <- object$samples$norm.factors
    N <- colSums(object$counts)
    f <- tmm * N/mean(tmm * N)
    norm.counts <- scale(object$counts, center=FALSE, scale=f)
    norm.counts <- removeNull(norm.counts)    
  }

  if (outfile) png(filename="figures/countsBoxplots.png",width=2*min(2200,1800+800*ncol(norm.counts)/10),height=1800,res=300)
	# raw counts
    varInt <- data.frame(varInt=group, sample=colnames(counts))
    counts <- counts %>% as.data.frame() %>% 
      mutate(ID=rownames(counts)) %>% 
      gather("sample", "count", -ID) %>% 
      full_join(varInt)
    norm.counts <- norm.counts %>% as.data.frame() %>% 
      mutate(ID=rownames(norm.counts)) %>% 
      gather("sample", "count", -ID) %>% 
      full_join(varInt)
    if (is.numeric(group)) {
      palette <- 15
      type <- "seq"
    } else {
      palette <- 6
      type = "qual"
    }
    countPlot <- ggplot(counts, aes(x=sample, y=log2(count+1), fill=factor(varInt), colour=factor(varInt))) + 
        geom_boxplot() +
      scale_colour_brewer(type = type, palette = palette) +
      scale_fill_brewer(type = type, palette = palette) +
      xlab ("") +
        fte_theme() +
        theme(axis.text.x = element_text(angle=90, size=round(400/nrow(varInt))),
              legend.position = "none")
    normPlot <- ggplot(norm.counts, aes(x=sample, y=log2(count+1), fill=factor(varInt), colour=factor(varInt))) + 
      geom_boxplot() +
      scale_colour_brewer(type = type, palette = palette) +
      scale_fill_brewer(type = type, palette = palette) +
      xlab ("") +
      fte_theme() +
      theme(legend.position="right") +
      theme(axis.text.x = element_text(angle=90, size=round(400/nrow(varInt))))
    grid.arrange(countPlot, normPlot, ncol=2)
    if (outfile) dev.off()
}
