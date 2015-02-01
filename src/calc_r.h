//
//  calc_r.h
//  iRel2
//
//  Created by Anders Goncalves da Silva on 17/04/13.
//
//

#ifndef iRel2_calc_r_h
double qr89(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, int n_alls);
double qr89_sum(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, int n_alls);
double lr99(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, int n_alls);
void prep_wang02(double *al_frq, int *al_per_loc, int n_loc, int n_alls, int n_ind, double *u, double *a2, double *a3, double *a4, double *a22, double *weights, int correct);
double wang02(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, double *a2, double *a3, double *a4, double *a22, double *weights);

#define iRel2_calc_r_h



#endif
