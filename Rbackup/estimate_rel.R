#### estimate relatedness indeces ####

#for testing
library(adegenet)
setwd("~/Dropbox/projects/irel/www_irel/")
data(nancycats)

estimate_rel<-function(data,input_allele_f = NULL){
  
  dyn.load("src/estimate_relR.so")
  
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
  
  relMatrix = matrix(0,ncol=(12+nLoc),nrow=totalDyads)
  
  out<-.Call("estimateRel",genoTable,nInd,alleleFreq,nAllPerLocus,nLoc,nAlleles,relMatrix,totalDyads)
  
  dyn.unload("src/estimate_relR.so")
  
 locusMissing = relMatrix[,13:(12+nLoc)]
 
 locusMissingF = apply(locusMissing,1,function(row) paste(locNames[!as.logical(row)],collapse='.'))
  
  relDF = data.frame(relMatrix[,1:12],stringsAsFactors = F)
  names(relDF)<-c('ind1', 'ind2','numbMissingLoci',"QG89_xy","QG89_yx","QG89_avg","QG89_rsxy","LR99_avg",'W02_unc','W02_cor','propAllelesShared','propLociShared')
  relDF[,"ind1"] <- inds[relDF$ind1]
  relDF[,"ind2"] <- inds[relDF$ind2]
  relDF$missingLoci <- locusMissingF

  return(relDF)
}

