#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cont_quad_knapsack.h"

double initial_multiplier_fixed_point(cqk_problem *restrict p, 
                                double *restrict abp,
                                double *restrict bbp,
                                double *restrict lr,
                                double *restrict ur,
                                double *restrict bl,
                                double *restrict bu) {

    double d1 = 0;
    double d2 = 0;
    double lambda = 0;
    for (unsigned i = 0; i < p->n; ++i) {
        lr[i] = (p->a[i] - p->d[i]*p->low[i])/p->b[i];
        ur[i] = (p->a[i] - p->d[i]*p->up[i])/p->b[i];
        abp[i] = (p->a[i]*p->b[i])/p->d[i];
        bbp[i] = (p->b[i]*p->b[i])/p->d[i];
        bl[i] = p->b[i]*p->low[i];
        bu[i] = p->b[i]*p->up[i];
        d1 += abp[i];
        d2 += bbp[i];
    }
    lambda = (d1 - p->r)/d2;
    return lambda;
}

void xis_quad(cqk_problem *restrict p, double lambda,
                double *restrict lr,double *restrict ur,
                double *restrict bl,double *restrict bu,
                double *restrict abp,double *restrict bbp,
                 double *Sbl, double *Sbu, double *Sbap, double *Sbbp){

    *Sbl = 0.0;
    *Sbu = 0.0;
    *Sbap = 0.0; 
    *Sbbp = 0.0;

    for (unsigned i = 0; i < p->n; i++){
        if (lambda >= lr[i]){
            *Sbl += bl[i];
        }
        else{
            if (lambda <= ur[i]){
                *Sbu += bu[i];
            }
            else{
                *Sbap = *Sbap + abp[i];
                *Sbbp = *Sbbp + bbp[i];
            }
        }
    }
}

void mount_x(cqk_problem *restrict p, double *restrict x, double lambda,
            double *restrict lr,double *restrict ur){

    for (unsigned i = 0; i < p->n; i++){
        if (lambda >= lr[i]){
            x[i] = p->low[i];
        }
        else{
            if (lambda <= ur[i]){
                x[i] = p->up[i];
            }
            else{
                x[i] = ( p->a[i] - lambda* p->b[i] ) / p->d[i];
            }
        }
    }
}

int fixed_point(cqk_problem *restrict p, double *x) {

    unsigned n_iters = 1;
    double lambda;
    double antLambda = 0;
    double Sbl;
    double Sbu;
    double Sbap;
    double Sbbp;
    
    unsigned n = p->n;

    double *restrict abp = (double *) malloc(n*sizeof(double));
    double *restrict bbp = (double *) malloc(n*sizeof(double));
    double *restrict lr = (double *) malloc(n*sizeof(double));
    double *restrict ur = (double *) malloc(n*sizeof(double));
    double *restrict bl = (double *) malloc(n*sizeof(double));
    double *restrict bu = (double *) malloc(n*sizeof(double));

    lambda = initial_multiplier_fixed_point(p,abp,bbp,lr,ur,bl,bu);
    while (n_iters <= MAXITERS) {
        xis_quad(p, lambda, lr,ur,bl,bu,abp,bbp, &Sbl, &Sbu, &Sbap, &Sbbp);
        
        antLambda = lambda;
        lambda = (Sbl + Sbu + Sbap - p->r)/(Sbbp);
        
        if(lambda - antLambda == 0){
            mount_x(p,x, lambda, lr,ur);
            break;
        }

        n_iters++;
    }

    free(abp);
    free(bbp);
    free(lr);
    free(ur);
    free(bl);
    free(bu);

    return n_iters;
}
