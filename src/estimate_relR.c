//
//  estimate_relR.c
//  iRel2
//
//  Created by Anders Goncalves da Silva on 17/04/13.
//
//

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "calc_r.h"
#include "rngs.h"


// function declarations
SEXP estimateRel(SEXP geno_table, SEXP n_ind, SEXP al_frq, SEXP n_al_per_l, SEXP n_loc, SEXP n_alls, SEXP out, SEXP totalDyads);

//functions

SEXP estimateRel(SEXP geno_table, SEXP n_ind, SEXP al_frq, SEXP n_al_per_l, SEXP n_loc, SEXP n_alls, SEXP out, SEXP totalDyads)
{

//declare counters
  int i, j, k, l, m, marker;

//translate from SEXP to C	
   double *genotable = REAL(geno_table);
   double *af = REAL(al_frq);
   int *napl = INTEGER(n_al_per_l);
   int nl = INTEGER(n_loc)[0];
   int na = INTEGER(n_alls)[0];
   int ni = INTEGER(n_ind)[0];
  
  //accounting for missing data
  double *newInd1, *newInd2, *newAllFrq;
  int newNLoc,newNAll,newMarker,pos1,pos2;
  int *newNAllPerLoc;
  int *missingLoci; // added on 24/7/14 to keep track of which loci are missing for each dyad
  
  //allocate memory
  newInd1 = (double *)malloc(sizeof(double)*na);
  newInd2 = (double *)malloc(sizeof(double)*na);
  newAllFrq = (double *)malloc(sizeof(double)*na);
  newNAllPerLoc = (int *)malloc(sizeof(int)*nl);
  missingLoci = (int *)malloc(sizeof(int)*nl);
  
  //counting total shared alleles and total shared alleles per locus
  // added on 24/7/14 to calculate the above quantities
  double propSharedAllelesPerLocus = 0.0;
  int countTotalSharedAlleles = 0;
  int countLocusSharedAlleles = 0;
  int flag = 0; // a flag to make sure a locus with shared alleles only gets counted once
  
  //variables for wang 2002
    //prep for uncorrected wang 2002
    double *u_unc, *a2_unc, *a3_unc, *a4_unc, *a22_unc, *weights_unc;
    int uncorrect = 0;
    
    u_unc = (double *)malloc(nl*sizeof(double));
    a2_unc = (double *)malloc(nl*sizeof(double));
    a3_unc = (double *)malloc(nl*sizeof(double));
    a4_unc = (double *)malloc(nl*sizeof(double));
    a22_unc = (double *)malloc(nl*sizeof(double));
    weights_unc = (double *)malloc(nl*sizeof(double));
    
    //prep for corrected wang 2002
    double *u, *a2, *a3, *a4, *a22, *weights;
    int correct = 1;
    
    u = (double *)malloc(nl*sizeof(double));
    a2 = (double *)malloc(nl*sizeof(double));
    a3 = (double *)malloc(nl*sizeof(double));
    a4 = (double *)malloc(nl*sizeof(double));
    a22 = (double *)malloc(nl*sizeof(double));
    weights = (double *)malloc(nl*sizeof(double));
	
  //keep track of the number of dyads
  int dyad, totaldyads;
  dyad = 0;
  totaldyads = asInteger(totalDyads);
  
  // dyads are formed by paring individual i with individual j
  for(i=0;i<(ni-1);i++){
    for(j=(i+1);j<ni;j++){
//       Rprintf("%d,%d\n",i,j);
      //set local variables for the current dyad
      newNLoc = 0;
      newNAll = 0;
      newMarker = 0;
      marker = 0;
      m = 0;
      pos1=pos2=0;
      countTotalSharedAlleles = 0;
      countLocusSharedAlleles = 0;
      
      //scrub memory
      memset(newInd1,0,sizeof(double)*na);
      memset(newInd2,0,sizeof(double)*na);
      memset(newAllFrq,0,sizeof(double)*na);
      memset(newNAllPerLoc,0,sizeof(int)*nl);
      memset(missingLoci,0,sizeof(int)*nl);
      //wang 2002 uncorrect parameters
      memset(u_unc,0,sizeof(double)*nl);
      memset(a2_unc,0,sizeof(double)*nl);
      memset(a3_unc,0,sizeof(double)*nl);
      memset(a4_unc,0,sizeof(double)*nl);
      memset(a22_unc,0,sizeof(double)*nl);
      memset(weights_unc,0,sizeof(double)*nl);
      //wang 2002 corrected parameters
      memset(u,0,sizeof(double)*nl);
      memset(a2,0,sizeof(double)*nl);
      memset(a3,0,sizeof(double)*nl);
      memset(a4,0,sizeof(double)*nl);
      memset(a22,0,sizeof(double)*nl);
      memset(weights,0,sizeof(double)*nl);

      
      //checking if pairs have missing loci, and making changes accordingly.
      // k iterates over loci
      for (k=0; k<nl; k++){
        pos1 = i+marker*ni;
        pos2 = j+marker*ni;
        flag = 0;
        if(ISNA(genotable[pos1]) || ISNA(genotable[pos2])){
            missingLoci[k] = 0;
            marker = marker + napl[k];
            continue;
        } else {
        	missingLoci[k] = 1;
            newNLoc = newNLoc + 1;
            newNAllPerLoc[(newNLoc-1)] = napl[k];
            newNAll = newNAll + napl[k];
        for(l=marker; l<(marker+napl[k]); l++){
            pos1 = i+l*ni;
            pos2 = j+l*ni;
            newAllFrq[m] = af[l];
            newInd1[m] = genotable[pos1];
            newInd2[m] = genotable[pos2];
            if(newInd1[m] > 0 && newInd2[m] > 0){
            	countTotalSharedAlleles += 1;
                if(flag == 0){
                    countLocusSharedAlleles += 1;
                    flag = 1;
                }
            }
            m++;
          }
         marker = marker + napl[k];
        }
    }
      int z; // to be used below to append data on presence/absence of a locus
//      lr99(newInd1, newInd2, newAllFrq,newNAllPerLoc,newNLoc,newNAll);
//      Rprintf("*********************\n");
//      Rprintf("%d,%d\n",newNLoc,newNAll);
//      for(z=0;z<newNAll;z++){
//        Rprintf("%f,%f,%f\n",newInd1[z],newInd2[z],newAllFrq[z]);
//      }
//      Rprintf("\n");
//      for(z=0;z<newNLoc;z++){
//        Rprintf("%d ",newNAllPerLoc[z]);
//      }
//      Rprintf("\n");
//      Rprintf("QG89 %d,%d,%d,%f\n",i,j,ni,lr99(newInd1, newInd2, newAllFrq,newNAllPerLoc,newNLoc,newNAll));
//      Rprintf("%d,%d\n",i,j);
//      Rprintf("*********************\n");
     
      //continue;
      //prep for wang 2002
      prep_wang02(newAllFrq, newNAllPerLoc, newNLoc, newNAll,ni,u_unc,a2_unc,a3_unc,a4_unc,a22_unc,weights_unc,uncorrect);
      prep_wang02(newAllFrq, newNAllPerLoc, newNLoc, newNAll,ni,u,a2,a3,a4,a22,weights,correct);
      
//      Rprintf("%f,%f,%f,%f\n",u[0],a2[0],a3[0],a4[0],a22[0]);
//      Rprintf("QG89 %f\n",qr89(newInd1, newInd2, newAllFrq,newNAllPerLoc, newNLoc, newNAll));
      //calculate indices for dyad
      REAL(out)[dyad+0*totaldyads] = i+1;
      REAL(out)[dyad+1*totaldyads] = j+1;
      REAL(out)[dyad+2*totaldyads] = nl-newNLoc;
      REAL(out)[dyad+3*totaldyads] = qr89(newInd1, newInd2, newAllFrq,newNAllPerLoc, newNLoc, newNAll);
      REAL(out)[dyad+4*totaldyads] = qr89(newInd2, newInd1, newAllFrq,newNAllPerLoc, newNLoc, newNAll);
      REAL(out)[dyad+5*totaldyads] = (REAL(out)[dyad+3*totaldyads]+REAL(out)[dyad+4*totaldyads])/2.0;
      REAL(out)[dyad+6*totaldyads] = qr89_sum(newInd2, newInd1, newAllFrq,newNAllPerLoc, newNLoc, newNAll);
      REAL(out)[dyad+7*totaldyads] = lr99(newInd1, newInd2, newAllFrq,newNAllPerLoc,newNLoc,newNAll);
      REAL(out)[dyad+8*totaldyads] = wang02(newInd1, newInd2, newAllFrq, newNAllPerLoc, newNLoc, a2_unc, a3_unc, a4_unc, a22_unc, weights_unc);
      REAL(out)[dyad+9*totaldyads] = wang02(newInd1, newInd2, newAllFrq, newNAllPerLoc, newNLoc, a2, a3, a4, a22, weights);
      REAL(out)[dyad+10*totaldyads] = (double)countTotalSharedAlleles/(2.0*(double)newNLoc);
      REAL(out)[dyad+11*totaldyads] = (double)countLocusSharedAlleles/(double)newNLoc;
     for(z=0; z<nl;z++){
      	REAL(out)[dyad+(12+z)*totaldyads] = missingLoci[z];
      }
      dyad = dyad + 1;
    }
  }
  return R_NilValue;
}