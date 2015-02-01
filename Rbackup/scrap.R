library(inline)
library(adegenet)
library(reshape2)
library(ggplot2)
library(data.table)
source("R/estimate_rel.R")
source("R/sim_rel.R")

#old code
# test <- cfunction(c(a = "data.frame",ind = 'numeric'), '
#   int i, pos, indIx,countNA;
#   countNA = 0;
#   indIx = asInteger(ind) - 1;
#   int nRows = Rf_nrows(a);
#   int nCols = Rf_ncols(a);
#   SEXP out = PROTECT(allocVector(REALSXP, nCols));
#   for(i=0;i<nCols;i++){
#       pos = indIx+i*nRows;
#       if(ISNA(REAL(a)[pos])){
#         countNA = countNA + 1;
#       }
#       REAL(out)[i] = REAL(a)[pos];
#   }
#   Rprintf("%d\\n",countNA);
#   UNPROTECT(1);
#   return out;
# ')
# 
# q89 <-cfunction(c(ind1 = 'numeric', ind2= 'numeric',al_frq= 'numeric', al_per_loc= 'numeric',n_loc= 'numeric',n_alls= 'numeric'), '
#     //Calculate a one way Queller and Goodnight 1989 index.
#     int marker, i, j, k;
#     marker = 0;
#     
#     double *ind1P, *ind2P, *alFrqP;
#     int *alPerLoc;
#     
#     ind1P = REAL(ind1);
#     ind2P = REAL(ind2);
#     alFrqP = REAL(al_frq);
#     alPerLoc = INTEGER(al_per_loc);
# 
#     int nLoc = asInteger(n_loc);
#     int nAlleles = asInteger(n_alls);
#     
#     SEXP num, dem, countMissingData, out;
#     num = PROTECT(allocVector(REALSXP,1));
#     dem = PROTECT(allocVector(REALSXP,1));
#     countMissingData = PROTECT(allocVector(REALSXP,1));
#     out = PROTECT(allocVector(VECSXP,6));
# 
#     memset(REAL(num),0,sizeof(double)*1);
#     memset(REAL(dem),0,sizeof(double)*1);
#     memset(REAL(countMissingData),0,sizeof(double)*1);
# 
#     //accounting for missing data
#     double *newInd1, *newInd2, *newAllFrq;
#     int newNLoc,newNAll,newMarker;
#     int *newNAllPerLoc;
#     
#     newNLoc = 0;
#     newNAll = 0;
#     newMarker = 0;
# 
#     newInd1 = (double *)malloc(sizeof(double)*nAlleles);
#     newInd2 = (double *)malloc(sizeof(double)*nAlleles);
#     newAllFrq = (double *)malloc(sizeof(double)*nAlleles);
#     newNAllPerLoc = (int *)malloc(sizeof(int)*nLoc);
# 
#     k=0;
#     for (i=0; i<nLoc; i++){
#       if(ISNA(ind1P[marker]) || ISNA(ind2P[marker])){
#           marker = marker + alPerLoc[i];
#           continue;
#       } else {
#           newNLoc = newNLoc + 1;
#           newNAllPerLoc[(newNLoc-1)] = alPerLoc[i];
#           newNAll = newNAll + alPerLoc[i];
#       for(j=marker; j<(marker+alPerLoc[i]); j++){
#           newAllFrq[k] = alFrqP[j];
#           newInd1[k] = ind1P[j];
#           newInd2[k] = ind2P[j];
#           k++;
#         }
#        marker = marker + alPerLoc[i];
#       }
#   }
#   Rprintf("New number of loci %d\\n",newNLoc);
#   Rprintf("New number of alleles %d\\n",newNAll);
#   for(i=0;i<newNLoc;i++){
#     for(j=0;j<newNAll;j++)
#       Rprintf("%d,%f,%f,%f\\n",i,newInd1[j],newInd2[j],newAllFrq[j]);
#   }
#   UNPROTECT(4);
# return R_NilValue;
# ')
# 
# 
# ind1 = nancycats@tab[1,]
# ind2 = nancycats@tab[2,]
# al_per_loc = nancycats@loc.nall
# allele_frq = apply(nancycats@tab,2,sum,na.rm=T)/
# apply(nancycats@tab,2,function(allele)
# sum(!is.na(allele)))
# n_loc = length(al_per_loc)
# n_alls = sum(al_per_loc)
# 
# q89(ind1,ind2,allele_frq,al_per_loc,n_loc,n_alls)
# 


# data
data(nancycats)
data = nancycats

inDex = c("QG89_xy","QG89_yx","QG89_avg","QG89_rsxy","LR99_avg",'W02_unc','W02_cor')
relCats = c("un","hs","fs")
nRelCats = length(relCats)

if(nRelCats == 2) {
  nGroups = 1
} else if (nRelCats > 2) {
  nGroups = nRelCats - 1
} else {
  stop("Not enough relatedeness categories!")
}

est=data.frame(estimate_rel(data))
est = data.table(est)
est = melt(est,id.vars=c(1,2,3,13))
setkeyv(est,c("missingLoci","variable"))


#groups by missing data
missD = unique(est[,missingLoci])
simData = list()
for(i in 1:length(missD)){
  #sim data by missing loci
  mD = missD[i]
  for(boot in 1:1){
    mD = missD[i]
    if(mD == ""){
      sims = sim_dataframe(data = data, relCategories = relCats,10000)
      sims = data.table(sims)
      setkeyv(sims,c("Categories","Index"))
      simData[["all"]] = sims
      missD[i] = "all"
      mD = missD[i]
    } else {
      sims = sim_dataframe(data = data, relCategories = relCats,10000,lociToKeep = mD)
      sims = data.table(sims)
      setkeyv(sims,c("Categories","Index"))
      simData[[mD]] = sims
    }
  
    #thresholds list
    thresholds = list()
    for(j in 1:length(inDex)){
      ix = inDex[j]
      
      thresholds[[ix]] <- matrix(0,ncol=nGroups, nrow = length(missD))
      tmp = thresholds[[ix]]
      #start threshold calculation loop here
      for(k in 1:nGroups){
        c1 = relCats[k]
        c2 = relCats[(k+1)]
        #estimate bloin's limits
        relSims1 = simData[[mD]][.(c1,ix)][,Value]
        relSims2 = simData[[mD]][.(c2,ix)][,Value]
        totalSims = length(relSims1)
        
        findMaximalOverlap = function(relValue,relSims1,relSims2,totalSims){
          pRel1 = sum(relSims1>=relValue)/totalSims
          pRel2 = sum(relSims2>=relValue)/totalSims
          return(pRel2-pRel1)
        }
        
        tmp[i,k] = optimise(findMaximalOverlap,relSims1=relSims1,relSims2=relSims2,totalSims=totalSims,interval=c(0,0.5),maximum = T)$maximum
        # finish threshold calculation loop here
      }
      tmp = c(-1.01,tmp[i,],1.01)
      s = simData[[mD]][.(relCats,ix)]
      s_class = cut(s[,Value],breaks = tmp, labels = relCats)
      missMatch = (s_class == s[,Categories])
      meanError = sum(missMatch)/nrow(s)
      print(c(mD,ix,meanError))
    }
  }
  
  e = est[.(mD,ix),]
  e_class = cut(e[,value],breaks = tmp, labels = relCats)
  print(table(e_class))
}

ggplot(simData[[mD]][.(relCats,ix)],aes(Value,fill=Categories))+geom_histogram(alpha=0.5,aes(y = ..density..), position = 'identity')+geom_vline(xintercept=tmp+0.02)
+geom_point(data = est,aes(x=ix,y=countMissingLoci/5,shape=as.factor(countMissingLoci)),fill=I("red"),guide=F)

missClass1 = sum(relSims1[,Value]>=boundary$maximum)/totalSims
missClass2 = sum(relSims2[,Value]<boundary$maximum)/totalSims

####
ind1 = c(1,1,0,0)/2
ind2 = c(0,0,2,0)/2
sum((ind1 - ind2)^2)


mean_h = sapply(lapply(seploc(nancycats), function(loc) apply(loc@tab, 1, function(g) {if(any(is.na(g))) {NA} else if (any(g == 0.5)) {1} else {0}})), function(h) sum(h, na.rm = T) / sum(!is.na(h)))


n_inds = length(nancycats@ind.names)

pairs = combn(x = 1:n_inds, m = 2)

calc_dij = function(pair){
  pairs_by_loc = seploc(pair)
  dij2 = sapply(pairs_by_loc, function(g) {
    genos = g@tab
    if (any(is.na(genos))) {
      NA
    } else {
      sum((genos[1,] - genos[2,])^2)
    }
  })
  not_miss = which(!is.na(dij2))
  mean(dij2[not_miss]) / mean(mean_h[not_miss])
}

mean_h = sapply(lapply(seploc(nancycats), function(loc) apply(loc@tab, 1, function(g) {if(any(is.na(g))) {NA} else if (any(g == 0.5)) {1} else {0}})), function(h) sum(h, na.rm = T) / sum(!is.na(h)))

nancy_rhk = apply(pairs, 2, function(inds) 1 - mean(calc_dij(nancycats[inds,]), na.rm = T))

##
kingroup_test = read.csv("data//kingroup_test.csv", header = T, stringsAsFactors = F)
kingroup_test$loc1 = gsub(pattern = "^\\s", replacement = "", kingroup_test$loc1)
kingroup_test$loc2 = gsub(pattern = "\\s?fast", replacement = "100", kingroup_test$loc2)
kingroup_test$loc2 = gsub(pattern = "\\s?slow1", replacement = "102", kingroup_test$loc2)
kingroup_test$loc2 = gsub(pattern = "\\s?slow2", replacement = "104", kingroup_test$loc2)
kingroup_test$loc2 = gsub(pattern = "\\s?slow3", replacement = "106", kingroup_test$loc2)
kingroup_test$loc2 = gsub(pattern = "\\s?slow4", replacement = "108", kingroup_test$loc2)
kingroup_test$loc3 = gsub(pattern = "\\s?big", replacement = "100", kingroup_test$loc3)
kingroup_test$loc3 = gsub(pattern = "\\s?small", replacement = "102", kingroup_test$loc3)
kingroup_test$loc4 = gsub(pattern = "\\s?a", replacement = "100", kingroup_test$loc4)
kingroup_test$loc4 = gsub(pattern = "\\s?b", replacement = "102", kingroup_test$loc4)
kingroup_test$loc4 = gsub(pattern = "\\s?c", replacement = "104", kingroup_test$loc4)
kingroup_test$loc4 = gsub(pattern = "\\s?d", replacement = "106", kingroup_test$loc4)
kingroup_test$loc4 = gsub(pattern = "\\s?e", replacement = "108", kingroup_test$loc4)
kingroup_test$loc4 = gsub(pattern = "\\s?f", replacement = "110", kingroup_test$loc4)
kingroup_test$loc4 = gsub(pattern = "\\s?g", replacement = "112", kingroup_test$loc4)
kingroup_test$loc5 = gsub(pattern = "\\s?all1", replacement = "100", kingroup_test$loc5)
kingroup_test$loc5 = gsub(pattern = "\\s?all2", replacement = "102", kingroup_test$loc5)
kingroup_test$loc6 = gsub(pattern = "\\s?com", replacement = "100", kingroup_test$loc6)
kingroup_test$loc6 = gsub(pattern = "\\s?rare", replacement = "102", kingroup_test$loc6)
kingroup_test$loc7 = gsub(pattern = "^\\s", replacement = "", kingroup_test$loc7)
kingroup_genind = df2genind(kingroup_test[,3:9], sep = "/", loc.names = names(kingroup_test)[3:9])

mean_h = sapply(lapply(seploc(kingroup_genind), function(loc) apply(loc@tab, 1, function(g) {if(any(is.na(g))) {NA} else if (any(g == 0.5)) {1} else {0}})), function(h) sum(h, na.rm = T) / sum(!is.na(h)))

n_inds = length(kingroup_genind@ind.names)

pairs = combn(x = 1:n_inds, m = 2)

apply(pairs[,1:10], 2, function(inds) 1 - mean(calc_dij(kingroup_genind[inds,]), na.rm = T))

test_rel = estimate_rel(data = kingroup_genind)

test_sim = irelr::sim_rel(data = kingroup_genind, reps = 100, k_vector = c(1,0,0))

#mariana
dat_mariana = read.genepop(file = "~/Downloads/inputGenepop.gen")
dat_mariana@ind.names <- paste("ind", 1:length(dat_mariana@ind.names), sep = "")
mean_h = sapply(lapply(seploc(dat_mariana), function(loc) apply(loc@tab, 1, function(g) {if(any(is.na(g))) {NA} else if (any(g == 0.5)) {1} else {0}})), function(h) sum(h, na.rm = T) / sum(!is.na(h)))

n_inds = length(dat_mariana@ind.names)

pairs = combn(x = 1:n_inds, m = 2)

mariana_rhk = apply(pairs, 2, function(inds) 1 - mean(calc_dij(dat_mariana[inds,]), na.rm = T))
