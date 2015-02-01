
## simulation functions --------------------------------------------------------

#'
#'
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
  res = sim_rvalues(n_ind,n_loc,n_alleles,allele_frq,n_als_per_locus,k_vector,reps)
  colnames(res)<-c("QG89_xy","QG89_yx","QG89_avg","QG89_rsxy","LR99_avg",'W02_unc','W02_cor','propAllelesShared','propLociShared')
  return(res)
    }

#' @useDynLib irelr sims
sim_rvalues <-function(n_ind, n_loc,n_alleles,allele_frq,n_alleles_per_locus,k_values,reps){
  n_ind = as.integer(n_ind)
  n_loc = as.integer(n_loc)
  n_alleles = as.integer(n_alleles)
  allele_frq = as.numeric(allele_frq)
  n_alleles_per_locus = as.integer(n_alleles_per_locus)
  k_values = as.numeric(k_values)
  reps = as.integer(reps)
  res = matrix(numeric(0),ncol=9,nrow=reps)
  out <- .Call("sims",n_ind,allele_frq,n_alleles_per_locus,k_values,n_loc,n_alleles,res,reps)
  return(res)
}


#'
#'
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