#' Wrapper to run DESeq2
#'
#' Wrapper to run DESeq2: create the \code{DESeqDataSet}, normalize data, estimate dispersions, statistical testing...
#'
#' @param counts \code{matrix} of raw counts
#' @param target target \code{data.frame} of the project
#' @param varInt name of the factor of interest (biological condition)
#' @param batch batch effect to take into account (\code{NULL} by default)
#' @param locfunc \code{"median"} (default) or \code{"shorth"} to estimate the size factors
#' @param fitType mean-variance relationship: "parametric" (default) or "local"
#' @param pAdjustMethod p-value adjustment method: \code{"BH"} (default) or \code{"BY"} for instance
#' @param cooksCutoff outliers detection threshold (TRUE to let DESeq2 choosing it or FALSE to disable the outliers detection)
#' @param independentFiltering \code{TRUE} or \code{FALSE} to perform the independent filtering or not
#' @param alpha significance threshold to apply to the adjusted p-values
#' @param ... optional arguments to be passed to \code{nbinomWaldTest()}
#' @return A list containing the \code{dds} object (\code{DESeqDataSet} class), the \code{results} objects (\code{DESeqResults} class) and the vector of size factors
#' @author Hugo Varet


run.DESeq2 <- function(counts, target, varInt, batch=NULL, interact=NULL, num_sva=0, kallisto=FALSE,
                       locfunc="median", fitType="parametric", pAdjustMethod="BH",
		       cooksCutoff=TRUE, independentFiltering=TRUE, alpha=0.05, ...){
  # building dds object
  if (kallisto) {
    dds <- DESeqDataSetFromTximport(counts, target, 
                                  formula(paste("~", varInt, 
                                                       ifelse(!is.null(interact), paste(c("", interact), collapse = " * "), ""),
                                                       ifelse(!is.null(batch), paste(c("", batch), collapse = " + "), "")
                                  )))
    dds<-DESeq(dds)
    
  } else {
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=target, 
                                design=formula(paste("~", varInt, 
                                                     ifelse(!is.null(interact), paste(c("", interact), collapse = " * "), ""),
                                                     ifelse(!is.null(batch), paste(c("", batch), collapse = " + "), "")
                                )))

  reduced=formula(paste("~ 1",  ifelse(!is.null(interact), paste(c("", interact), collapse = " + "), ""),
                        ifelse(!is.null(batch), paste(c("", batch), collapse = " + "), "")
  ))
  # normalization
  dds <- estimateSizeFactors(dds,locfunc=eval(as.name(locfunc)))
  cat("\nNormalization factors:\n")
  print(sizeFactors(dds))

  
  # surrogate variable analysis
  if (num_sva > 0){
    library(sva)
    dat<-counts(dds, normalized=TRUE)
    idx<-rowMeans(dat) > 1
    dat<-dat[idx,]
    svseq<-svaseq(dat, model.matrix(design(dds), colData(dds)), model.matrix(reduced, colData(dds)), n.sv=num_sva)
    surVar <- as.data.frame(svseq$sv)
    colData(dds)<-cbind(colData(dds), surVar)
    design(dds) <- formula(paste("~ ", varInt, 
                                   ifelse(!is.null(interact), paste(c("", interact), collapse = " * "), ""),
                                   ifelse(!is.null(batch), paste(c("", batch), collapse = " + "), ""),
                                   paste(c("", colnames(surVar)), collapse = " + ")
    ))
  }  

  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)),collapse=" "),"\n")					  
  
    # estimating dispersions
  dds <- estimateDispersions(dds, fitType=fitType)
  
  # statistical testing: perform all the comparisons between the levels of varInt
  dds <- nbinomWaldTest(dds, ...)
}
  results <- list()
  for (comp in combn(nlevels(colData(dds)[,varInt]), 2, simplify=FALSE)){
    levelRef <- levels(colData(dds)[,varInt])[comp[1]]
    levelTest <- levels(colData(dds)[,varInt])[comp[2]]
    results[[paste0(levelTest,"_vs_",levelRef)]] <- results(dds, contrast=c(varInt, levelTest, levelRef),
                                                            pAdjustMethod=pAdjustMethod, cooksCutoff=cooksCutoff,
                                                            independentFiltering=independentFiltering, alpha=alpha)
    cat(paste("Comparison", levelTest, "vs", levelRef, "done\n"))
  }
  
  return(list(dds=dds,results=results,sf=sizeFactors(dds)))
}
