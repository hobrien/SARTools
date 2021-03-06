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
#' @param testMethod test method: \code{"Wald"} (default) or \code{"LRT"} for instance
#' @param cooksCutoff outliers detection threshold (TRUE to let DESeq2 choosing it or FALSE to disable the outliers detection)
#' @param independentFiltering \code{TRUE} or \code{FALSE} to perform the independent filtering or not
#' @param alpha significance threshold to apply to the adjusted p-values
#' @param ... optional arguments to be passed to \code{nbinomWaldTest()}
#' @return A list containing the \code{dds} object (\code{DESeqDataSet} class), the \code{results} objects (\code{DESeqResults} class) and the vector of size factors
#' @author Hugo Varet

run.DESeq2.LRT <- function(counts, target, varInt, batch=NULL, interact=NULL, reduced=NULL,
                       locfunc="median", fitType="parametric", pAdjustMethod="BH", kallisto=TRUE,
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
  
  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)),collapse=" "),"\n")					  

  cat("Design of the reduced model:\n")
  cat(paste(as.character(reduced),collapse=" "),"\n")					  
  
  # normalization
  dds <- estimateSizeFactors(dds,locfunc=eval(as.name(locfunc)))
  cat("\nNormalization factors:\n")
  print(sizeFactors(dds))
  
  # estimating dispersions
  dds <- estimateDispersions(dds, fitType=fitType)
  
  dds <- nbinomLRT(dds, reduced=reduced)
  }
  results <- list()
  results[[paste0("drop_", varInt)]]<- results(dds, pAdjustMethod=pAdjustMethod, cooksCutoff=cooksCutoff,
                                                            independentFiltering=independentFiltering, alpha=alpha, name=varInt)
  
  return(list(dds=dds,results=results,sf=sizeFactors(dds)))
}
