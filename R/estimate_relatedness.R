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
#'    \item{\code{QG89_xy}: Queller and Goodnight (1989) index with \code{ind1} as the numerator and \code{ind2} as the denominator.}
#'    \item{\code{QG89_yx}: Queller and Goodnight (1989) index with \code{ind2} as the numerator and \code{ind1} as the denominator.}
#'    \item{\code{QG89_avg}: the mean of \code{QG89_xy} and \code{QG89_yx}.}
#'    \item{\code{QG89_rsxy}: sum of the numerators of \code{QG89_xy} and \code{QG89_yx} divided by the sum of the denominators of \code{QG89_xy} and \code{QG89_yx}.}
#'    \item{\code{LR99_avg}: the mean of the Lynch and Ritland (1999) index calculated with \code{ind1} as the reference individual, and  \code{ind2} as the reference.}
#'    \item{\code{W02_unc}: the uncorrected Wang (2002) index.}
#'    \item{\code{W02_cor}: the corrected Wang (2002) index.}
#'    \item{\code{prop_allele_shared}: proportion of alleles that are shared between \code{ind1} and \code{ind2}.}
#'    \item{\code{prop_loci_shared}: proportion of loci that are that have at least one allele in common between \code{ind1} and \code{ind2}.}
#'    \item{\code{missing_loci}: id of missing loci.}
#'  }
#'
#' @details 
#' 
#' Queller and Goodnight (1989)
#' 
#' \deqn{r_{QG89_{xy}} = \sum_{l}\frac{(\sum_{alleles}(p_{xm} - \bar{p_{m}}))}{(\sum_{alleles}(p_{ym} - \bar{p_{m}}))}}{r_{QG89_{xy}} = sum_{l}((sum(p_xm - p_m)) / sum_{l}(sum(p_ym - p_m)))}
#'  
#' \deqn{r_{QG89_{yx}} = \sum_{l}\frac{(\sum_{alleles}(p_{ym} - \bar{p_{m}}))}{(\sum_{alleles}(p_{xm} - \bar{p_{m}}))}}{r_{QG89_{yx}} = sum_{l}((sum(p_ym - p_m)) / sum_{l}(sum(p_xm - p_m)))}
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
  geno_table = adegenet::tab(data, freq = T)
  n_ind = adegenet::nInd(data)
  inds = adegenet::indNames(data)
  if(is.null(input_allele_f)){
    allele_freq = as.numeric(apply(geno_table,2,sum,na.rm=T)/
                              apply(geno_table,2,function(allele)
                                sum(!is.na(allele))))
  } else {
    allele_freq = as.numeric(input_allele_f)
  }
  n_all_per_locus = as.integer(adegenet::nAll(data))
  n_loc = as.integer(adegenet::nLoc(data))
  loc_names = as.character(adegenet::locNames(data))
  n_alleles = as.integer(sum(n_all_per_locus))
  total_dyads = as.integer((n_ind*(n_ind-1)/2))
  heteroz = adegenet::summary(data, verbose = F)$Hobs
  rel_matrix = matrix(0,ncol=(13+n_loc),nrow=total_dyads)
  out<-.Call("estimateRel", 
             geno_table, 
             n_ind,allele_freq, 
             n_all_per_locus, 
             n_loc, 
             n_alleles, 
             rel_matrix, 
             total_dyads, 
             heteroz)
  locus_missing = rel_matrix[,13:(13+n_loc)]
  locus_missing_f = apply(locus_missing,1,function(row) paste(loc_names[!as.logical(row)],collapse='.'))
  rel_df = data.frame(rel_matrix[,1:13], 
                      stringsAsFactors = F)
  names(rel_df)<-c('ind1', 
                   'ind2',
                   'n_missing_loci', 
                   "QG89_xy", 
                   "QG89_yx", 
                   "QG89_avg", 
                   "QG89_rsxy", 
                   "LR99_avg", 
                   "W02_unc", 
                   "W02_cor", 
                   "HK08", 
                   "prop_alleles_shared", 
                   "prop_loci_shared")
  rel_df[,"ind1"] <- inds[rel_df$ind1]
  rel_df[,"ind2"] <- inds[rel_df$ind2]
  rel_df$missing_loci <- locus_missing_f
  return(rel_df)
}
