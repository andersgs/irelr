//
//  calc_r.c
//  iRel2
//
//  Created by Anders Goncalves da Silva on 17/04/13.
//
//

//#define N 778
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int i, j, reps, marker, count_x, count_y;

double qr89(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, int n_alls);
double qr89_sum(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, int n_alls);
double lr99(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, int n_alls);
void prep_wang02(double *al_frq, int *al_per_loc, int n_loc, int n_alls, int n_ind, double *u, double *a2, double *a3, double *a4, double *a22, double *weights, int correct);
double wang02(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, double *a2, double *a3, double *a4, double *a22, double *weights);
void wang_phi_delta(double p1, double p2, double p3, double b, double c, double d, double e, double f, double g, double *phi, double *delta);
double hk08(double *ind1, double *ind2, int *al_per_loc, int n_loc, double *h);

double qr89(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, int n_alls){
    //Calculate a one way Queller and Goodnight 1989 index.
    marker = 0;
    double num=0.0, dem=0.0;
    
    for (i=0; i<n_loc; i++){
        for(j=marker; j<(marker+al_per_loc[i]); j++){
            if(ind1[j]){
              //Rprintf("In QG %f,%f,%f\n", ind1[j], ind2[j], al_frq[j]);
                num = num + (ind2[j] - al_frq[j]);
                dem = dem + (ind1[j] - al_frq[j]);
            }
        }
        marker = marker + al_per_loc[i];
    }
    return num/dem;
}

double qr89_sum(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, int n_alls){
    //Calculate a one way Queller and Goodnight 1989 index.
    marker = 0;
    double num=0.0, dem=0.0;
    
    for (i=0; i<n_loc; i++){
        for(j=marker; j<(marker+al_per_loc[i]); j++){
            if(ind1[j]){
                num = num + (ind2[j] - al_frq[j]);
                dem = dem + (ind1[j] - al_frq[j]);
            }
        }
        marker = marker + al_per_loc[i];
    }
    
    marker = 0;
    for (i=0; i<n_loc; i++){
        for(j=marker; j<(marker+al_per_loc[i]); j++){
            if(ind2[j]){
                num = num + (ind1[j] - al_frq[j]);
                dem = dem + (ind2[j] - al_frq[j]);
            }
        }
        marker = marker + al_per_loc[i];
    }

    return num/dem;
}


double lr99(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, int n_alls){
    
    int tmpa_x, tmpb_x, tmpc_y, tmpd_y;
    
    double pa, pb;
    short int Sab, Sbc, Sbd, Sac, Sad;
    double rXY=0.0, rYX=0.0, sum_weights_XY = 0.0, sum_weights_YX = 0.0, estimate = 0.0;
    double *loc_r_XY, *loc_weight_XY, *loc_r_YX, *loc_weight_YX;
    loc_r_XY = (double *)malloc(sizeof(double)*n_loc);
    loc_weight_XY = (double *)malloc(sizeof(double)*n_loc);
    
    loc_r_YX = (double *)malloc(sizeof(double)*n_loc);
    loc_weight_YX = (double *)malloc(sizeof(double)*n_loc);
    
    
    marker = 0;
    
    
    for(i=0; i<n_loc; i++){
        count_x = count_y = 0;
        for(j=marker; j<(marker+al_per_loc[i]);j++){
            if(ind1[j]){
                if(count_x == 0){
                    if(ind1[j] == 1){
                        tmpa_x = j;
                        tmpb_x = j;
                    } else {
                        tmpa_x = j;
                        count_x = count_x + 1;
                    }
                } else {
                    tmpb_x = j;
                }
            }
            if(ind2[j]){
                if(count_y == 0){
                    if(ind2[j] == 1){
                        tmpc_y = j;
                        tmpd_y = j;
                    } else {
                        tmpc_y = j;
                        count_y = count_y + 1;
                    }
                } else{
                    tmpd_y = j;
                }
            }
        }
        marker = marker + al_per_loc[i];
        //X as reference
        Sab = tmpa_x == tmpb_x ? 1 : 0;
        Sac = tmpa_x == tmpc_y ? 1 : 0;
        Sad = tmpa_x == tmpd_y ? 1 : 0;
        Sbc = tmpb_x == tmpc_y ? 1 : 0;
        Sbd = tmpb_x == tmpd_y ? 1 : 0;
        pa = al_frq[tmpa_x];
        pb = al_frq[tmpb_x];
        
        loc_r_XY[i] = (pa*(Sbc+Sbd) + pb*(Sac+Sad) - 4*pa*pb)/((1+Sab)*(pa+pb) - 4*pa*pb);
        loc_weight_XY[i] = ((1+Sab)*(pa+pb) - 4*pa*pb)/(2*pa*pb);
        
        //Y as reference
        Sab = tmpc_y == tmpd_y ? 1 : 0;
        Sac = tmpc_y == tmpa_x ? 1 : 0;
        Sad = tmpc_y == tmpb_x ? 1 : 0;
        Sbc = tmpd_y == tmpa_x ? 1 : 0;
        Sbd = tmpd_y == tmpb_x ? 1 : 0;
        pa = al_frq[tmpc_y];
        pb = al_frq[tmpd_y];
        
        loc_r_YX[i] = (pa*(Sbc+Sbd) + pb*(Sac+Sad) - 4*pa*pb)/((1+Sab)*(pa+pb) - 4*pa*pb);
        loc_weight_YX[i] = ((1+Sab)*(pa+pb) - 4*pa*pb)/(2*pa*pb);
    }
    
    for(i=0; i<n_loc; i++){
        rXY = rXY + loc_r_XY[i]*loc_weight_XY[i];
        sum_weights_XY = sum_weights_XY + loc_weight_XY[i];
        rYX = rYX + loc_r_YX[i]*loc_weight_YX[i];
        sum_weights_YX = sum_weights_YX + loc_weight_YX[i];
    }
    
    free(loc_r_XY);
    free(loc_weight_XY);
    free(loc_r_YX);
    free(loc_weight_YX);
    
    estimate = ((rXY/sum_weights_XY)+(rYX/sum_weights_YX))/2.0;
    
    //with few loci, the average can be slightly larger than 1 or smaller than -1
    
    if(estimate <= 1.0 & estimate >= -1){
        return estimate;
    }
    else if(estimate > 1.0){
        return 1.0;
    }
    else {
        return -1.0;
    }
}

void prep_wang02(double *al_frq, int *al_per_loc, int n_loc, int n_alls, int n_ind, double *u, double *a2, double *a3, double *a4, double *a22, double *weights, int correct)
{
    int N = 2*n_ind;
    marker = 0;
    double u_total = 0.0;
    for(i=0; i<n_loc;i++){
        for(j=marker; j<(marker+al_per_loc[i]);j++){
            a2[i] = a2[i] + pow(al_frq[j],2.0);
            a3[i] = a3[i] + pow(al_frq[j],3.0);
            a4[i] = a4[i] + pow(al_frq[j],4.0);
        }
        marker = marker + al_per_loc[i];
    }
    if(correct){
        for(i=0; i<n_loc; i++){
            a2[i] = (N*a2[i] - 1)/(N - 1);
            a3[i] = ((pow(N,2.0)*a3[i]) - 3*(N-1)*a2[i] - 1)/((N-1)*(N-2));
            a4[i] = ((pow(N,3.0)*a4[i]) - (6*(N-1)*(N-2)*a3[i]) - (7*(N-1)*a2[i]) - 1)/(pow(N,3.0) - 6*pow(N,2.0) + (11*N) - 6);
        }
    }

    for(i=0; i<n_loc; i++){
        u[i] = 2*a2[i] - a3[i];
        u_total = u_total + 1/u[i];
    }
    
    for(i=0; i<n_loc; i++){
        weights[i] = (1/(u_total*u[i]));
    }
    
    for(i=0; i<n_loc; i++){
        a22[i] = pow(a2[i],2.0);
    }
}

double wang02(double *ind1, double *ind2, double *al_frq, int *al_per_loc, int n_loc, double *a2, double *a3, double *a4, double *a22, double *weights)
{
    double p1, p2, p3;
    p1 = p2 = p3 = 0;
    double b, c, d, e, f, g;
    b = c = d = e = f = g = 0;
    double a2_bar, a3_bar, a4_bar, a22_bar;
    a2_bar = a3_bar = a4_bar = a22_bar = 0;
    double phi, delta;
    int sim;
    marker = 0;
    double total_w = 0;
    
    for(i=0; i<n_loc; i++){
        sim = 0;
        for(j=marker; j<(marker+al_per_loc[i]); j++){
            if(ind1[j] == 0 && ind2[j] == 0){
                continue;
            }
            if(ind1[j] == 1.0 && ind2[j] == 1.0){
                sim = 4;
                break;
            }
            if(ind1[j] == 1.0 && ind2[j] == 0.5){
                sim = 3;
                break;
            }
            if(ind1[j] == 0.5 && ind2[j] == 1.0){
                sim = 3;
                break;
            }
            if(ind1[j] == 0.5 && ind2[j] == 0.5){
                sim = sim + 2;
                continue;
            }
        }
        marker = marker + al_per_loc[i];

        if(sim == 4){
            p1 = p1 + weights[i];
        }
        if(sim == 3){
            p2 = p2 + weights[i];
        }
        if(sim == 2){
            p3 = p3 + weights[i];
        }
        a2_bar = a2_bar + a2[i]*weights[i];
        a3_bar = a3_bar + a3[i]*weights[i];
        a4_bar = a4_bar + a4[i]*weights[i];
        a22_bar = a22_bar + a22[i]*weights[i];
        total_w = total_w + weights[i];
    }
    
    p1 = p1/total_w;
    p2 = p2/total_w;
    p3 = p3/total_w;
    a2_bar = a2_bar/total_w;
    a3_bar = a3_bar/total_w;
    a4_bar = a4_bar/total_w;
    a22_bar = a22_bar/total_w;

    b = 2*a22_bar - a4_bar;
    c = a2_bar - 2*a22_bar + a4_bar;
    d = 4*(a3_bar - a4_bar);
    e = 2*(a2_bar - 3*a3_bar + 2*a4_bar);
    f = 4*(a2_bar - a22_bar - 2*a3_bar + 2*a4_bar);
    g = (1 - 7*a2_bar + 4*a22_bar + 10*a3_bar - 8*a4_bar);
    
    wang_phi_delta(p1,p2,p3,b,c,d,e,f,g,&phi,&delta);
    
    return ((phi/2)+delta);
}

void wang_phi_delta(double p1, double p2, double p3, double b, double c, double d, double e, double f, double g, double *phi, double *delta)
{
    double V;
    
    V = ((pow((1-b),2.0)*((pow(e,2.0)*f) + d*pow(g,2.0)))) - ((1-b)*(pow((e*f - d*g),2.0))) + (2*c*d*f*(1-b)*(g+e)) + (pow(c,2.0)*d*f*(d+f));
    
    *phi = ((d*f*((e+g)*(1-b) + c*(d+f)))*(p1-1) + (d*(1-b)*((g*(1-b-d) + f*(c+e))))*(p3) + (f*(1-b)*((e*(1-b-f) + d*(c+g))))*(p2))/V;
    
    *delta = ((c*d*f*(e+g)*(p1+1-(2*b))) + (((1-b)*(f*pow(e,2.0)+d*pow(g,2.0))) - pow((e*f - d*g),2.0))*(p1-b) + c*(d*g - e*f)*(d*p3 - f*p2) - pow(c,2.0)*d*f*(p3 + p2 - d - f) - c*(1-b)*(d*g*p3 + e*f*p2))/V;

}

double hk08(double *ind1, double *ind2, int *al_per_loc, int n_loc, double *h) {
    marker = 0;
    double dij2 = 0.0;
    double total_h = 0.0;
    double mean_dij, mean_h;
    
    for (i=0; i<n_loc; i++){
        for(j=marker; j<(marker+al_per_loc[i]); j++){
          dij2 = dij2 + pow((ind1[j] - ind2[j]), 2.0);
        }
        marker = marker + al_per_loc[i];
        total_h = total_h + h[i];
    } 
    mean_dij = dij2 / n_loc;
    mean_h = total_h / n_loc;
  return  1 - (mean_dij / mean_h) ;
}
