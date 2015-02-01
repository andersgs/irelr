#blouin classifier
# dat = nancycats
# 
# est = estimate_rel(dat)
# 
# relCat = c("un","hs","fs","po")
# iters = 10000
# missL = sort(unique(est$missingLoci))
# if(missL[1] != ''){
# 	missL = append(missL, 'all', after = 0)
# } else {
# 	missL[1] = 'all'
# }
# 
# sims = lapply(missL,function(sim) {
# 	sim = sim_dataframe(dat,relCat,iters,sim)
# 	sim = data.table(sim)
# 	setkeyv(sim,c("Categories","Index"))
# 	return(sim)
# 	}
# )
# 
# names(sims) <- missL
# 
# indices = c("QG89_xy","QG89_yx","QG89_avg","QG89_rsxy","LR99_avg",'W02_unc','W02_cor')
# 
# combCats = combn(relCat,2)
# 
# classifyBlouin = function(d_relEstimates, d_sims = NULL, d_data = NULL, d_alFreq = NULL ,l_kLowOrderVectors = NULL ,l_kFirstOrderVectors = NULL, nIters = 10000){
# 	
# 	#estimate relatedness categories if needed
# 	if(is.null(d_relEstimates)){
# 		d_relEstimates = estimate_rel(d_data,d_alFreq)
# 	}
# 	
# 	#figure out missing loci
# 	missLoci = sort(unique(d_relEstimates$missingLoci))
# 	
# 	#figure out what thresholds to calculate
# 	relCats = c(names(l_kLowOrderVectors),names(l_kFirstOrderVectors))
# 	nRelCats = length(relCats)
# 	if(nRelCats == 4){
# 		nThresholds = 4
# 		classLabels = c("un","hs","fo")
# 	} else if (nRelCats == 3){
# 		if(sum(relCats %in% c('fs','po'))==2){
# 			nThresholds = 3
# 			classLabels = c(xx,"fo")
# 		} else {
# 			nThresholds = 4
# 			classLabels = relCats
# 		}
# 	} else if (nRelCats == 2){
# 		if(sum(relCats %in% c('fs','po'))==2){
# 		
# 		} else if(sum(relCats %in% c('fs','po'))==2){
# 		
# 		} else if(sum(relCats %in% c('fs','po'))==2){
# 		
# 		}
# 		
# 	}
# 
# 	boundaries = matrix(0,cols = (nRelCats), rows = length(missLoci) )  # plus 1 cols because there are nRelCats-1 cols plus two for -1 and 1, the lower most and upper most boundaries	
# 	for(l in 1:length(l_missLocCats)){
# 		for(i in 1:nRelCats){
# 				boundaries[l,] = -1
# 				boundaries[l,nRelCats+2] = 1
# 				boundaries[l,i+1] = findBoundary(d_sims,orderedRelCats[i],orderedRelCats[i+1], index, l_missLocCats[l])
# 		}
# 	}
# } 
# 
# orderRelCats <- function(l_kVectors){
# 	#order cats by expected number of IBD alleles
# 	expectedNumbIBD = lapply(l_kVectors,function(x) sum(x*c(0,1,2)))
# 	firstOrderCats = names(expectedNumbIBD)[expectedNumbIBD >= 1]
# 	lowerOrderCats = names(expectedNumbIBD)[expectedNumbIBD < 1]
# 	orderFirstOrderCats = order(firstOrderCats)
# 	orderLowerOrderCats = order(as.numeric(expectedNumbIBD[lowerOrderCats]))
# 	return(list(lower=lowerOrderCats[orderLowerOrderCats],upper=firstOrderCats[orderFirstOrderCats]))
# }
# 
# findBoundary = function(d_sims,cat1,cat2, index, locCat){
# 	#find boundary that minimises overlap between distributions of simulated relatedeness indices
# 	relSims1 = d_sims[[locCat]][.(cat1,index)]
# 	relSims2 = d_sims[[locCat]][.(cat2,index)]
# 	totalSims = dim(relSims1)[1]
# 	
# 	boundary=optimise(findMaximalOverlap,
# 		relSims1=relSims1,relSims2=relSims2,
# 		totalSims=totalSims,interval=c(0,0.5),maximum = T)
# 	
# 	return(boundary$maximum)
# 
# }
# 
# for (i in indices){
# 	for (j in 1:(dim(combCats)[2]-1)){
# 		for (l in missL){
# 		
# 			cat1 = combCats[1,j]
# 			cat2 = combCats[2,j]
# 			ix = i
# 			locCat = l
# 			
# 			#estimate bloin's limits
# 			relSims1 = sims[[locCat]][.(cat1,ix)]
# 			relSims2 = sims[[locCat]][.(cat2,ix)]
# 			totalSims = dim(relSims1)[1]
# 			
# 			boundary=optimise(findMaximalOverlap,
# 				relSims1=relSims1,relSims2=relSims2,
# 				totalSims=totalSims,interval=c(0,0.5),maximum = T)
# 			missClass1 = sum(relSims1[,Value] >= boundary$maximum)/totalSims
# 			missClass2 = sum(relSims2[,Value] < boundary$maximum)/totalSims
# 			
# 			print(c(ix,locCat,cat1,cat2,boundary$maximum,missClass1))
# 			print(c(ix,locCat,cat2,cat1,boundary$maximum,missClass2))
# 
# 		}
# 	}
# }
# ggplot(sims[['all']][.(c("un","po"),"QG89_rsxy")],aes(Value,fill=Categories))+geom_histogram(alpha=0.5,aes(y = ..density..), position = 'identity')+geom_vline(xintercept=boundary$maximum+0.01666667)
# +geom_point(data = est,aes(x=QG89_rsxy,y=countMissingLoci/5,shape=as.factor(countMissingLoci)),fill=I("red"),guide=F)
# 
# findMaximalOverlap = function(relValue,relSims1,relSims2,totalSims){
#   pRel1 = sum(relSims1[,Value]>=relValue)/totalSims
#   pRel2 = sum(relSims2[,Value]>=relValue)/totalSims
#   return(pRel2-pRel1)
# }
