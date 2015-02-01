### data summary ###

#alleles per locus
alleleCountPlot <- function(data){

  lociDF = data.frame(locusNames = as.character(data$loc.names), 
                      numbAlleles = as.numeric(data$loc.nall)
  )
  
 p = ggplot(lociDF,aes(x=locusNames, y=numbAlleles))+geom_bar(stat='identity',fill='#008cba') +
    xlab("Loci") +
    ylab("Count of alleles") +
    theme(axis.title = element_text(size=18),
          axis.text = element_text(size=14,
                                   colour='black'),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1))
 return(p)
}

#missing data

missingDataPlot <- function(data){
  
  bar_col = function(countMissingData,numbInd,threshold = 0.05){
    maxMissingData = numbInd * threshold
    isOverMax = sapply(countMissingData,function(count) if(count < maxMissingData) {'#008cba'} else {'#f04124'})
    return(isOverMax)
  }
  
  locusData=seploc(data)
  
  countMissingData = sapply(locusData, function(loc) sum(apply(loc$tab,1,function(r) any(is.na(r)))))
  
  missingDF = data.frame(locusNames = as.character(data$loc.names),
                         missingData = countMissingData)
  
  barColors = bar_col(countMissingData,length(data$ind.names),0.05)
  
  p = ggplot(missingDF,aes(x=locusNames,y=missingData))+geom_bar(stat="identity",fill=barColors)+
    ylim(c(0,length(data$ind.names))) +
    xlab("Loci") +
    ylab("Number of individuals") +
    theme(axis.title = element_text(size=18),
          axis.text = element_text(size=14,
                                   colour='black'),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1))
  return(p)
}  

