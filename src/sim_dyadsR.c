//
//  sim_dyads.c
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

int i, j, z, max, rep;
int start,end;
int a, b, c, d, tmp; // tmp allele positions
double indicator, prob1, prob2;

int share[3] = {0};
// int counts[13] = {0};
// int total = 0;


// function declarations
int al_pos(double *al_frq, int start, int end);
SEXP sims(SEXP n_ind, SEXP al_frq, SEXP n_al_per_l, SEXP k_array, SEXP n_loc, SEXP n_alls, SEXP rep_out, SEXP reps, SEXP heteroz);
void sim_dyad(double *al_frq, int *n_al_per_l, double *k_array, int n_loc, double *ind1, double *ind2);
void sim_noIBD(double *al_frq, int start, int end, double *ind1, double *ind2);
void sim_oneIBD(double *al_frq, int start, int end, double *ind1, double *ind2);
void sim_twoIBD(double *al_frq, int start, int end, double *ind1, double *ind2);


//functions

SEXP sims(SEXP n_ind, SEXP al_frq, SEXP n_al_per_l, SEXP k_array, SEXP n_loc, SEXP n_alls, SEXP rep_out, SEXP reps, SEXP heteroz)
{
	
//translate to SEXP to C	
  	double *af = REAL(al_frq);
  	int *napl = INTEGER(n_al_per_l);
  	double *k = REAL(k_array);
	  double *h = REAL(heteroz);
    int nl = INTEGER(n_loc)[0];
	  int na = INTEGER(n_alls)[0];
    int ni = INTEGER(n_ind)[0];
	
	
    double *sim_ind1,*sim_ind2;
    sim_ind1 = (double *)malloc(sizeof(double)*na);
    sim_ind2 = (double *)malloc(sizeof(double)*na);
    SelectStream(0);
    PutSeed(time(NULL));

	// a variable to count shared alleles
	int countTotalSharedAlleles;
    int countLociShared;
    int flag;
    int i_loc;
    int i_al;

    //prep for uncorrected wang 2002
    double *u_unc, *a2_unc, *a3_unc, *a4_unc, *a22_unc, *weights_unc;
    int uncorrect = 0;
    
    u_unc = (double *)malloc(nl*sizeof(double));
    memset(u_unc,0,sizeof(double)*nl);
    a2_unc = (double *)malloc(nl*sizeof(double));
    memset(a2_unc,0,sizeof(double)*nl);
    a3_unc = (double *)malloc(nl*sizeof(double));
    memset(a3_unc,0,sizeof(double)*nl);
    a4_unc = (double *)malloc(nl*sizeof(double));
    memset(a4_unc,0,sizeof(double)*nl);
    a22_unc = (double *)malloc(nl*sizeof(double));
    memset(a22_unc,0,sizeof(double)*nl);
    weights_unc = (double *)malloc(nl*sizeof(double));
    memset(weights_unc,0,sizeof(double)*nl);

    prep_wang02(af, napl, nl, na,ni,u_unc,a2_unc,a3_unc,a4_unc,a22_unc,weights_unc,uncorrect);

    //prep for uncorrected wang 2002
    double *u, *a2, *a3, *a4, *a22, *weights;
    int correct = 1;
    
    u = (double *)malloc(nl*sizeof(double));
    memset(u,0,sizeof(double)*nl);
    a2 = (double *)malloc(nl*sizeof(double));
    memset(a2,0,sizeof(double)*nl);
    a3 = (double *)malloc(nl*sizeof(double));
    memset(a3,0,sizeof(double)*nl);
    a4 = (double *)malloc(nl*sizeof(double));
    memset(a4,0,sizeof(double)*nl);
    a22 = (double *)malloc(nl*sizeof(double));
    memset(a22,0,sizeof(double)*nl);
    weights = (double *)malloc(nl*sizeof(double));
    memset(weights,0,sizeof(double)*nl);

    prep_wang02(af, napl, nl, na,ni,u,a2,a3,a4,a22,weights,correct);


    for(rep=0; rep<INTEGER(reps)[0]; rep++){
        memset(sim_ind1,0,sizeof(double)*na);
        memset(sim_ind2,0,sizeof(double)*na);
        countTotalSharedAlleles = 0;
        countLociShared = 0;
        z=0;
        sim_dyad(af, napl, k, nl, &sim_ind1[0], &sim_ind2[0]);
        
        for (i_loc=0; i_loc<nl; i_loc++){
            flag = 0;
                for(i_al=z; i_al<(z+napl[i_loc]); i_al++){
                    if(sim_ind1[i_al] > 0 && sim_ind2[i_al] > 0){
                        countTotalSharedAlleles += 1;
                        if(flag == 0){
                            countLociShared += 1;
                            flag = 1;
                        }
                    }
                }
                z = z + napl[i_loc];
            }
//    Rprintf("%d, %d\n",countTotalSharedAlleles, countLociShared);
    
//        for(z=0; z<na ; z++){
//        	if(sim_ind1[z] > 0 && sim_ind2[z] > 0){
//        		countTotalSharedAlleles += 1;
//        	}
//        }
        REAL(rep_out)[rep+0*INTEGER(reps)[0]] = qr89(sim_ind1, sim_ind2, af,napl, nl, na);
        REAL(rep_out)[rep+1*INTEGER(reps)[0]] = qr89(sim_ind2, sim_ind1, af,napl, nl, na);
        REAL(rep_out)[rep+2*INTEGER(reps)[0]] = (REAL(rep_out)[rep+0*INTEGER(reps)[0]]+REAL(rep_out)[rep+1*INTEGER(reps)[0]])/2.0;
        REAL(rep_out)[rep+3*INTEGER(reps)[0]] = qr89_sum(sim_ind2, sim_ind1, af,napl, nl, na);
        REAL(rep_out)[rep+4*INTEGER(reps)[0]] = lr99(sim_ind1, sim_ind2, af,napl,nl,na);
        REAL(rep_out)[rep+5*INTEGER(reps)[0]] = wang02(sim_ind1, sim_ind2, af, napl, nl, a2_unc, a3_unc, a4_unc, a22_unc, weights_unc);
        REAL(rep_out)[rep+6*INTEGER(reps)[0]] = wang02(sim_ind1, sim_ind2, af, napl, nl, a2, a3, a4, a22, weights);
        REAL(rep_out)[rep+7*INTEGER(reps)[0]] = hk08(sim_ind1, sim_ind2, napl, h);
        REAL(rep_out)[rep+8*INTEGER(reps)[0]] = (double)countTotalSharedAlleles/(2.0*(double)nl);
        REAL(rep_out)[rep+9*INTEGER(reps)[0]] = (double)countLociShared/((double)nl);
    }
//        Rprintf("%d, %d, %d\n",share[0], share[1], share[2]);
        
        return R_NilValue;

}


void sim_dyad(double *al_frq, int *n_al_per_l, double *k_array, int n_loc, double *ind1, double *ind2)
/* main simulation function
Takes as input a memory address for two individual vectors of length n_allels, and simulates their genotypes in accordance with k_matrix
*/

{
    //Rprintf("Got to here\n");
    int count = 0;
    for(i=0; i<n_loc; i++){
        indicator = Random();
        start = count;
        count = count + n_al_per_l[i];
        end = count - 1;
        if(indicator <= k_array[0]){
            share[0] = share[0]+1;
            sim_noIBD(al_frq, start, end, ind1, ind2);
            continue;
        }
        if(indicator > k_array[0] && indicator <= (k_array[0]+k_array[1])){
            share[1] = share[1]+1;
            sim_oneIBD(al_frq, start, end, ind1, ind2);
            continue;
        }
        if(indicator > (k_array[0]+k_array[1])){
            share[2] = share[2]+1;
            sim_twoIBD(al_frq, start, end, ind1, ind2);
            continue;
        }
    }
}

void sim_noIBD(double *al_frq, int start, int end, double *ind1, double *ind2){
    a = al_pos(al_frq, start, end);
    b = al_pos(al_frq, start, end);
    c = al_pos(al_frq, start, end);
    d = al_pos(al_frq, start, end);
    if (a == b){
        *(ind1+a) = 1.0;
    } else {
        *(ind1+a) = 0.5;
        *(ind1+b) = 0.5;
    }
    if (c == d){
        *(ind2+c) = 1.0;
    } else {
        *(ind2+c)= 0.5;
        *(ind2+d) = 0.5;
    }
}

void sim_oneIBD(double *al_frq, int start, int end, double* ind1, double* ind2){
    a = al_pos(al_frq, start, end);
    b = al_pos(al_frq, start, end);
    c = a;
    d = al_pos(al_frq, start, end);
    if (a == b){
        *(ind1+a) = 1.0;
    } else {
        *(ind1+a) = 0.5;
        *(ind1+b) = 0.5;
    }
    if (c == d){
        *(ind2+c) = 1.0;
    } else {
        *(ind2+c)= 0.5;
        *(ind2+d) = 0.5;
    }
}

void sim_twoIBD(double *al_frq, int start, int end, double *ind1, double *ind2){
    a = al_pos(al_frq, start, end);
    b = al_pos(al_frq, start, end);
    c = a;
    d = b;
    if (a == b){
        ind1[a] = 1.0;
    } else {
        ind1[a] = 0.5;
        ind1[b] = 0.5;
    }
    if (c == d){
        ind2[c] = 1.0;
    } else {
        ind2[c]= 0.5;
        ind2[d] = 0.5;
    }
}

int al_pos(double *al_frq, int start, int end)
{
    tmp = start;
    prob1 = Random();
    prob2 = al_frq[start];
    while(prob1>prob2){
        tmp = tmp + 1;
        prob2 = prob2 + al_frq[tmp];
    }
//     counts[tmp] = counts[tmp] + 1;
//     total = total + 1;
    
//    printf("%d, %d, %d, %f, %f\n",start, tmp, end, prob1, prob2);
    return tmp;
}
