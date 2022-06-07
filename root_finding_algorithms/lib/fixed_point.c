#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fixed_point.h"

/*
 * Function Name: allocate_cqk_problem
 *
 * Description: Allocate the memory necessary hold a CQN problem in
 *     dimension n.
 *
 * Input:
 *     unsigned n: The dimension of the problem.
 * 
 * Output: cqk_problem *p: pointer to the cqk_problem whose memory is
 *     being allocated.
 */
void allocate_cqk_problem(unsigned n, cqk_problem *restrict p) {
    p->n = n;
    p->d = (double *) malloc(p->n*sizeof(double));
    p->a = (double *) malloc(p->n*sizeof(double));
    p->b = (double *) malloc(p->n*sizeof(double));
    p->low = (double *) malloc(p->n*sizeof(double));
    p->up = (double *) malloc(p->n*sizeof(double));
    if (!p->d || !p->a || !p->b || !p->low || !p->up) {
        fprintf(stderr, "Memory allocation error, line %d, file %s\n",
                __LINE__, __FILE__);
        exit(1);
    }
}

/*
 * Function Name: free_cqk_problem
 *
 * Description: Free the memory associated to a cqk_problem.
 *
 * Input/Output: 
 *     cqk_problem *p: at entry a cqk_problem with allocated
 *         memory, at output all the memory will be freed and the pointers
 *         in the cqk_problem structure will point to NULL.
 */
void free_cqk_problem(cqk_problem *restrict p) {
    p->n = 0;
    free(p->d);
    p->d = NULL;
    free(p->a);
    p->a = NULL;
    free(p->b);
    p->b = NULL;
    free(p->low);
    p->low = NULL;
    free(p->up);
    p->up = NULL;
}

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
