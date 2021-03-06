#' Wrapper to run edgeR
#'
#' Wrapper to run edgeR: create the \code{dge} object, normalize data, estimate dispersions, statistical testing...
#'
#' @param counts \code{matrix} of counts
#' @param target target \code{data.frame} of the project
#' @param varInt name of the factor of interest (biological condition)
#' @param condRef reference biological condition
#' @param batch batch effect to take into account (\code{NULL} by default)
#' @param cpmCutoff counts-per-million cut-off to filter low counts
#' @param normalizationMethod normalization method: \code{"TMM"} (default), \code{"RLE"} (DESeq) or \code{"upperquartile"}
#' @param pAdjustMethod p-value adjustment method: \code{"BH"} (default) or \code{"BY"}
#' @param ... optional arguments to be passed to \code{glmFit()}
#' @return A list containing the \code{dge} object and the \code{results} object
#' @author Hugo Varet

run.voom <- function(counts, target, varInt, condRef, batch=NULL, cpmCutoff=1, 
                      normalizationMethod="TMM", pAdjustMethod="BH", ...){
  
  # filtering: select features which contain at least 
  # minReplicates (smallest number of replicates) with
  # at least cpmCutoff counts per million
  minReplicates <- min(table(target[,varInt]))
  fcounts <- counts[rowSums(cpm(counts) >= cpmCutoff) >= minReplicates,]
  cat("Number of features discarded by the filtering:\n")
  cat(nrow(counts)-nrow(fcounts),"\n")
  
  # building dge object
  design <- formula(paste("~", varInt, 
                       ifelse(!is.null(batch), paste(c("", batch), collapse = " + "), "")
                        )
                    )
                       
  dge <- DGEList(counts=fcounts, remove.zeros=TRUE)
  dge$design <- model.matrix(design, data=target)
  cat("\nDesign of the statistical model:\n")
  cat(paste(as.character(design),collapse=" "),"\n")					  
  
  # normalization
  dge <- calcNormFactors(dge, method=normalizationMethod)
  cat("\nNormalization factors:\n")
  print(dge$samples$norm.factors)

  v <- voom(dge, design=dge$design, plot=TRUE)
  fit <- lmFit(v)
  fit <- eBayes(fit)
  dir.create("tables")
  AllResults <-topTable(fit, coef=2, number=5000, sort.by="P", p.value=.1)
  AllResults$Id<-rownames(AllResults)
  MaleUp <- filter(AllResults, logFC>0)
  write_delim(MaleUp, "tables/MaleUp.txt", delim='\t')
  FemaleUp <-filter(AllResults, logFC<0)
  write_delim(FemaleUp, "tables/FemaleUp.txt", delim='\t')
  cat(paste("Coefficients of the model:",paste(colnames(fit$design),collapse="  ")),"\n")
  colsToTest <- grep(varInt,colnames(fit$design))
  namesToTest <- paste0(gsub(varInt,"",colnames(fit$design)[colsToTest]),"_vs_",condRef)
  results <- list()
  # testing coefficients individually (tests againts the reference level)
  for (i in 1:length(colsToTest)){
    cat(paste0("Comparison ",gsub("_"," ",namesToTest[i]),": testing coefficient ",colnames(fit$design)[colsToTest[i]]),"\n")
    lrt <- glmLRT(fit, coef=colsToTest[i])
    results[[namesToTest[i]]] <- topTags(lrt,n=nrow(dge$counts),adjust.method=pAdjustMethod,sort.by="none")$table
  }
  # defining contrasts for the other comparisons (if applicable)
  if (length(colsToTest)>=2){
    colnames <- gsub(varInt,"",colnames(fit$design))
    for (comp in combn(length(colsToTest),2,simplify=FALSE)){ 
  	contrast <- numeric(ncol(dge$design))
  	contrast[colsToTest[comp[1:2]]] <- c(-1,1)
  	namecomp <- paste0(colnames[colsToTest[comp[2]]],"_vs_",colnames[colsToTest[comp[1]]])
  	cat(paste0("Comparison ",gsub("_"," ",namecomp),": testing contrast (",paste(contrast,collapse=", "),")"),"\n")
      lrt <- glmLRT(fit, contrast=contrast)
      results[[namecomp]] <- topTags(lrt,n=nrow(dge$counts),adjust.method=pAdjustMethod,sort.by="none")$table
    }
  }
  
  return(list(dge=dge,results=results))
}
