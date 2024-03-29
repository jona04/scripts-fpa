/*
 * Module Name: fixed_point
 *
 * Description: Implements the Fixed point method described in
 *    the manuscript 
 *
 *    A. Alves, J.O.L. Silva, L.C. Matioli, P.S.M. Santos, S.S. Souza. 
 *    "A fixed-point algorithm for solving quadratic convex separable 
 *     knapsack problems".
 *
 *
 * Copyright: Jonatas O. L. da Silva <jonatas.iw@gmail.com> 2022.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cont_quad_knapsack.h"

/********** Fixed point method Separable Quadratic Knapsack problem
            described in the manuscript

    A. Alves, J.O.L. Silva, L.C. Matioli, P.S.M. Santos, S.S. Souza. 
    "A fixed-point algorithm for solving quadratic convex separable 
    knapsack problems".

**********/ 

/*
 * Function Name: initial_multiplier_fixed_point
 *
 * Description: Estimate an initial multiplier from an approximate
 *     solution.
 *
 * Input:
 *     cqk_problem *restrict p: description of the cont. quad. knapsack.
 *
 * Output:
 *      double *restrict abp: terms resultant of the operation (a*b)/d.
 *      double *restrict bbp: terms resultant of the operation (b*b)/d.
 *      double *restrict lr: terms resultant of the operation (a-d*low)/b.
 *      double *restrict ur: terms resultant of the operation (a-d*up)/b.
 *      double *restrict bl: terms resultant of the operation b*low.
 *      double *restrict bu: terms resultant of the operation b*up.
 * 
 * Return value: the desired estimate.
 * 
 * Obs: The multiplier is computed as suggested by the variable fixing
 * methods.
 */
double initial_multiplier_fixed_point_new(cqk_problem *restrict p, 
                                double *r,
                                double *restrict abp,
                                double *restrict bbp,
                                double *restrict ur,
                                double *restrict lr,
                                double *restrict bu) {

    double d2 = 0.0;
    double d1 = 0.0;
    double lambda = 0.0;
    double pre = 0.0;
    double aa, uu;
    for (unsigned i = 0; i < p->n; ++i) {

        aa = p->a[i] - (p->d[i] * p->low[i]);
        uu = p->up[i] - p->low[i];
        ur[i] = (aa - p->d[i] * uu)/p->b[i];
        lr[i] = aa/p->b[i];
        pre = p->b[i] / p->d[i];
        abp[i] = aa * pre;
        bbp[i] = p->b[i] * pre;
        bu[i] = p->b[i] * uu;
        d1 += abp[i];
        d2 += bbp[i];
        *r = *r - (p->b[i] * p->low[i]);
    }
    lambda = (d1 - *r)/d2;
    return lambda;
}

/*
 * Function Name: xis_quad
 *
 * Description: Compute sum of the term bl, bu, abp, bbp presented
 *              in its respectives sets necessary to calculate the new multiplier.
 *
 * Input:
 *     cqk_problem *restrict p: description of the cont. quad. knapsack.
 *     double lambda: current multiplier.     
 *     double *restrict lr: terms contained the upper bounds for the multiplier.
 *     double *restrict ur: terms contained the lower bounds for the multiplier. 
 *
 * Output:
 *      double *Sbl: sum of terms bl presented in set I_l.
 *      double *Sbu: sum of terms bu presented in set I_u.
 *      double *Sbap: sum of terms abp presented in set I_eq.
 *      double *Sbbp: sum of terms bbp presented in set I_eq.
 */
void xis_quad_new(unsigned n,
                double *restrict lambda,
                double *restrict ur,  double *restrict lr,
                double *restrict bu,
                double *restrict abp,double *restrict bbp,
                double *Sbu, double *Sbap, double *Sbbp){

    *Sbu = 0.0;
    *Sbap = 0.0; 
    *Sbbp = 0.0;

    for (unsigned i = 0; i < n; i++){
        
        if (*lambda <= ur[i]){
            *Sbu += bu[i];
            continue;
        }
        
        if (*lambda <= lr[i]){
            *Sbap += abp[i];
            *Sbbp += bbp[i];
            continue;
        }

    }
}

/*
 * Function Name: mount_x
 *
 * Description: Mount the final result x based on the final multiplier and
 *              box constraints.
 *
 * Input:
 *      cqk_problem *restrict p: description of the cont. quad. knapsack.
 *      double lambda: final multiplier.     
 *      double *restrict lr: terms contained the upper bounds for the multiplier.
 *      double *restrict ur: terms contained the lower bounds for the multiplier.
 *
 * Output:
 *     double *x: solution obtained.
 */
void mount_x_new(cqk_problem *restrict p, double *restrict x, double lambda,
            double *restrict ur, double *restrict lr){
    
    for (unsigned i = 0; i < p->n; i++){

        if (lambda >= lr[i])
            x[i] = p->low[i];
        else if (lambda <= ur[i])
            x[i] = p->up[i];
            else 
                x[i] = (( p->a[i] - lambda * p->b[i] ) / p->d[i]);
    }
}

/*
 * Function Name: fixed_point
 *
 * Description: Fixedpoint method for the sep. quad. knapsack problem.
 *
 * Input:
 *     cqk_problem *restrict p: the problem description.
 * 
 * Output:
 *     double *x: solution obtained.
 * 
 * Return value: the number of iterations required to compute the
 *     solution, -1 in case of failure, -2 in case of infeasible
 *     problem.
 */
int fixed_point_new(cqk_problem *restrict p, double *x) {

    unsigned n_iters = 1;
    double lambda, antLambda;
    double Sbu, Sbap, Sbbp;

    unsigned n = p->n;
    double r = p->r;

    double *restrict abp = (double *) malloc(n*sizeof(double));
    double *restrict bbp = (double *) malloc(n*sizeof(double));
    double *restrict ur = (double *) malloc(n*sizeof(double));
    double *restrict lr = (double *) malloc(n*sizeof(double));
    double *restrict bu = (double *) malloc(n*sizeof(double));

    lambda = initial_multiplier_fixed_point_new(p,&r,abp,bbp,ur,lr,bu);
    while (n_iters <= MAXITERS) {
        xis_quad_new(n, &lambda,ur,lr,bu,abp,bbp, &Sbu, &Sbap, &Sbbp);
        
        antLambda = lambda;
        lambda = (Sbu + Sbap - r)/(Sbbp);

        if(fabs(lambda - antLambda) <= PREC){
            mount_x_new(p,x, lambda,ur,lr);
            break;
        }

        n_iters++;
    }

    free(abp);
    free(bbp);
    free(ur);
    free(lr);
    free(bu);

    return n_iters;
}
