## calculate relatedness functions ---------------------------------------------

#' Estimate relatedness values
#' 
#' @param data A GENIND class object
#' @param input_allele_f A vector of allele frequencies (optional)
#' 
#' @return A data.frame with \eqn{(n * (n - 1)) \ 2} (where \code{n} is the 
#'          number of genotypes) rows of pairwise relatedness values for a 
#'          number of relatedness indices.
#'  
#'  \itemize{
#'    \item{\code{ind1, ind2}: Pair of individuals for which relatedness values
#'        refer to.}
#'    \item{\code{n_missing_loci}: number of missing loci in the two individuals.}
#'    \item{\code{QG89_xy}: Queller and Goodnight (1989) index with \code{ind1}
#'                          as the numerator and \code{ind2} as the denominator.#'                          
#'                          \deqn{
#'                            \begin{equation}
#'                              r_{QG89_yx} = \sum_{l}\frac{(\sum_{alleles}(p_{xm} - \bar{p_{m}}))}
#'                                                        {(\sum_{alleles}(p_{ym} - \bar{p_{m}}))}

#'                            \end{equation}
#'                          }{
#'                          rQG89_xy = 
#'                              sum_{l}((sum(p_xm - p_m)) / sum_{l}(sum(p_ym - p_m)))
#'                          }
#'                          }
#'    \item{\code{QG89_yx}: Queller and Goodnight (1989) index with \code{ind2}
#'                          as the numerator and \code{ind1} as the denominator.
#'                          \deqn{
#'                            \begin{equation}
#'                              r_{QG89_yx} = \sum_{l}\frac{(\sum_{alleles}(p_{ym} - \bar{p_{m}}))}
#'                                                        {(\sum_{alleles}(p_{xm} - \bar{p_{m}}))}
#'                            \end{equation}
#'                          }{
#'                          rQG89_yx =                          
#'                              sum_{l}((sum(p_ym - p_m)) / sum_{l}(sum(p_xm - p_m)))
#'                          }
#'                          }
#'    \item{\code{QG89_avg}: the mean of \code{QG89_xy} and \code{QG89_yx}.}
#'    \item{\code{QG89_rsxy}: sum of the numerators of \code{QG89_xy} and 
#'                          \code{QG89_yx} divided by the sum of the denominators 
#'                          of \code{QG89_xy} and \code{QG89_yx}.}
#'    \item{\code{LR99_avg}: the mean of the Lynch and Ritland (1999) index 
#'                          calculated with \code{ind1} as the reference individual
#'                          , and  \code{ind2} as the reference.}
#'    \item{\code{W02_unc}: the uncorrected Wang (2002) index.}
#'    \item{\code{W02_cor}: the corrected Wang (2002) index.}
#'    \item{\code{prop_allele_shared}: proportion of alleles that are shared
#'                          between \code{ind1} and \code{ind2}.}
#'    \item{\code{prop_loci_shared}: proportion of loci that are that have at 
#'                          least one allele in common
#'                          between \code{ind1} and \code{ind2}.}
#'    \item{\code{missing_loci}: id of missing loci.}
#'  }
#' 
#' @references Queller and Goodnight. 1989.
#' @references Lynch and Ritland. 1999.
#'
#' @useDynLib irelr estimateRel
#'
#' @export
#' 
#' @examples
#' library(adegenet)
#' data(nancycats)
#' estimate_rel(data = nancycats)

estimate_rel<-function(data,input_allele_f = NULL){
  genoTable = data@tab  
  nInd = as.integer(length(data@ind.names))
  inds = as.character(data@ind.names)
  if(is.null(input_allele_f)){
    alleleFreq = as.numeric(apply(genoTable,2,sum,na.rm=T)/
                              apply(genoTable,2,function(allele)
                                sum(!is.na(allele))))
  } else {
    alleleFreq = as.numeric(input_allele_f)
  }
  nAllPerLocus = as.integer(data@loc.nall)
  nLoc = as.integer(length(nAllPerLocus))
  locNames = as.character(data@loc.names)
  nAlleles = as.integer(sum(nAllPerLocus))
  totalDyads = as.integer((nInd*(nInd-1)/2))
  heteroz = sapply(lapply(seploc(data), function(loc) apply(loc@tab, 1, function(g) {if(any(is.na(g))) {NA} else if (any(g == 0.5)) {1} else {0}})), function(h) sum(h, na.rm = T) / sum(!is.na(h)))
  relMatrix = matrix(0,ncol=(13+nLoc),nrow=totalDyads)
  out<-.Call("estimateRel",genoTable,nInd,alleleFreq,nAllPerLocus,nLoc,nAlleles,relMatrix,totalDyads,heteroz)
  locusMissing = relMatrix[,13:(13+nLoc)]
  locusMissingF = apply(locusMissing,1,function(row) paste(locNames[!as.logical(row)],collapse='.'))
  relDF = data.frame(relMatrix[,1:13],stringsAsFactors = F)
  names(relDF)<-c('ind1', 'ind2','n_missing_loci',"QG89_xy","QG89_yx","QG89_avg","QG89_rsxy","LR99_avg",'W02_unc','W02_cor','HK08','prop_alleles_shared','prop_loci_shared')
  relDF[,"ind1"] <- inds[relDF$ind1]
  relDF[,"ind2"] <- inds[relDF$ind2]
  relDF$missing_loci <- locusMissingF
  return(relDF)
}
