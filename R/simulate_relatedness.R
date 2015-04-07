## simulation functions --------------------------------------------------------

#' Simulate relatedness values
#'
#'@description This function takes as input a GENIND object, and simulates
#'      \code{reps} dyads related according to the specifications in the 
#'      the \code{k_vector}. It outputs a data.frame with relatedness values
#'      for each simulated dyad for each of the indices described in the 
#'      \code{estimate_rel} function. A vector of allele frequencies can be
#'      optionally supplied. If none is supplied, allele frequencies are 
#'      estimated from the data.
#'      
#'@param data An object of class GENIND containing genotypes for a single population.
#'@param reps The number of dyads to simulate.
#'@param k_vector A vector of three numbers that must add to 1.0 that gives the
#'                expected proportion of loci that are under one of the three
#'                possible IBD modes for non-inbred populations (see Details).
#'@param allele_frq A vector of size \code{n} alleles (total number of alleles
#'                  across all loci), containing the frequency of each allele in
#'                  the population. Should be ordered as in the \code{tab} slot
#'                  of the GENIND object. This parameter is optional (see Details).
#'                  
#'@details The goal of this function is to take population allele frequency data
#'          and calculate relatedness values calculated from simulated pairs of 
#'          individuals (dyads) that have any given 
#'          relatedness level (assuming an outbred population). The distribution
#'          of relatedness values can then be used to help classify observed
#'          pairs of sampled individuals into relatedness categories (e.g.,
#'          Blouin et al. 1996).
#'
#'@references Blouin M, Parsons M, Lacaille V, Lotz S (1996) Use of microsatellite 
#'                  loci to classify individuals by relatedness. 
#'                  Molecular Ecology, 5, 393–401.
#'@references Wang JL (2011) Unbiased relatedness estimation in structured 
#'                  populations. Genetics, 187, 887–901.
#'
#'@export
sim_rel <- function(data, reps, k_vector = c(1.0,0.0,0.0), allele_frq = NULL){
  #test for genind
  if(class(data)!="genind"){
    stop("Data must be of class genind. 
         Please look at the package adegenet 
         to transform the data to the 
         appropriate format.")
  }
  #determine number of individuals
  n_ind = length(data$ind.names)
  #determine number of loci
  n_loc = length(data@loc.names)
  #number of alleles per locus
  n_als_per_locus = data@loc.nall
  #total number of alleles
  n_alleles = sum(n_als_per_locus)
  #calculate the allele frequency based on the data 
  # - corrected for missing data
  # or if supplied by user, check that the data
  # match what is expected from the data
  if(is.null(allele_frq)){
    allele_frq = apply(data@tab,2,sum,na.rm=T)/
      apply(data@tab,2,function(allele)
        sum(!is.na(allele)))
  } else {
    if(length(allele_frq != n_alleles)){
      stop("Number of alleles provided in the allele_frq is 
           different from the number of alleles 
           in your data")
    }
    }
  if(sum(k_vector)!=1.0){
    stop("k-vector values must sum to 1.0.")
  }
  if(length(k_vector)!= 3){
    stop("k-vector must have three values.")
  }
  heteroz = sapply(lapply(seploc(data), function(loc) apply(loc@tab, 1, function(g) {if(any(is.na(g))) {NA} else if (any(g == 0.5)) {1} else {0}})), function(h) sum(h, na.rm = T) / sum(!is.na(h)))
  res = sim_rvalues(n_ind,n_loc,n_alleles,allele_frq,n_als_per_locus,k_vector,reps, heteroz)
  colnames(res)<-c("QG89_xy","QG89_yx","QG89_avg","QG89_rsxy","LR99_avg",'W02_unc','W02_cor', 'hk08','propAllelesShared','propLociShared')
  return(res)
    }

#' @useDynLib irelr sim_relC
sim_rvalues <-function(n_ind, n_loc,n_alleles,allele_frq,n_alleles_per_locus,k_values,reps, heteroz){
  n_ind = as.integer(n_ind)
  n_loc = as.integer(n_loc)
  n_alleles = as.integer(n_alleles)
  allele_frq = as.numeric(allele_frq)
  n_alleles_per_locus = as.integer(n_alleles_per_locus)
  k_values = as.numeric(k_values)
  reps = as.integer(reps)
  heteroz = as.numeric(heteroz)
  res = matrix(numeric(0),ncol=10,nrow=reps)
  out <- .Call("sim_relC",n_ind,allele_frq,n_alleles_per_locus,k_values,n_loc,n_alleles,res,reps,heteroz)
  return(res)
}

sim_rel_dataframe = function(data,relCategories = c('po', 'fs', 'hs', 'un'),numbIter = 10000, lociToKeep = 'all', all_freqs = NULL){
  relFactor = factor(rep(relCategories,each=9*numbIter),levels=relCategories)
  if(lociToKeep != 'all'){
    missLoc = strsplit(lociToKeep,"\\.")[[1]]
    data = data[loc=levels(data@loc.fac)[!locNames(data) %in% missLoc]]
  }
  simData = sapply(relCategories,function(cat) 
    melt(sim_rel(data=data,reps=numbIter,k_vector=defaultKVectors[[cat]], allele_frq = all_freqs))
  )
  simDataDF=data.frame(Categories = relFactor, Index = factor(unlist(simData[2,])),Value = unlist(simData[3,]))
  return(simDataDF)
}