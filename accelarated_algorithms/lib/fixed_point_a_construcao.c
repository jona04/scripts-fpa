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
#include "../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h"

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
double initial_multiplier_fixed_point_a(cqk_problem *restrict p, 
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
void xis_quad_a(unsigned n,
                double lambda,
                double *restrict ur,  double *restrict lr,
                double *restrict bu,
                double *restrict abp,double *restrict bbp,
                double *Sbu, double *Sbap, double *Sbbp){

    *Sbu = 0.0;
    *Sbap = 0.0; 
    *Sbbp = 0.0;

    for (unsigned i = 0; i < n; i++){
        
        if (lambda < ur[i]){
            *Sbu += bu[i];
            continue;
        }
        
        if (lambda < lr[i]){
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
void mount_x_a(cqk_problem *restrict p, double *restrict x, double lambda,
            double *restrict ur, double *restrict lr){
    
    for (unsigned i = 0; i < p->n; i++){

        if (lambda > lr[i])
            x[i] = p->low[i];
        else if (lambda < ur[i])
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
int fixed_point_a(cqk_problem *restrict p, double *x) {

    unsigned k = 0;
    double gammaK, x0, xk, gxk, tempvec, A;
    double residual = 1e9;
    double Sbu;
    double Sbap;
    double Sbbp;
    
    unsigned mk;
    int m = 1;

    unsigned n = p->n;
    double r = p->r;

    double *restrict Fk = (double *) malloc(100*sizeof(double));
    double *restrict alpha = (double *) malloc(100*sizeof(double));
    double *restrict g_mk_bulk = (double *) malloc(100*sizeof(double));

    double *restrict abp = (double *) malloc(n*sizeof(double));
    double *restrict bbp = (double *) malloc(n*sizeof(double));
    double *restrict ur = (double *) malloc(n*sizeof(double));
    double *restrict lr = (double *) malloc(n*sizeof(double));
    double *restrict bu = (double *) malloc(n*sizeof(double));

    x0 = initial_multiplier_fixed_point_a(p,&r,abp,bbp,ur,lr,bu);
    
    xis_quad_a(n, x0,ur,lr,bu,abp,bbp, &Sbu, &Sbap, &Sbbp);
    xk = (Sbu + Sbap - r)/(Sbbp);

    while (k <= MAXITERS) {
        mk = MIN(k, m);

        xis_quad_a(n, xk,ur,lr,bu,abp,bbp, &Sbu, &Sbap, &Sbbp);
        gxk = (Sbu + Sbap - r)/(Sbbp);

        residual = xk - gxk;

        if (fabs(residual) <= PREC){
            mount_x_a(p,x, gxk,ur,lr);
            break;
        }

        if (k > m){
            g_mk_bulk[k] = gxk;
            Fk[k] = residual;
            for (unsigned i =0; i < m; i++){
                Fk[i] = Fk[i+1];
                g_mk_bulk[i] = g_mk_bulk[i+1];
            }
        }else {
            Fk[k] = residual;
            g_mk_bulk[k] = gxk;
        }

        tempvec = Fk[mk] * -1.0;
        if (mk == 0){
            A = Fk[0];
            alpha[0] = 1;
        }else{
            
        }
        gammaK = gg[k] / Gk[0];
        result2 = (Xk[0] + Gk[0]) * gammaK;
        xx[k + 1] = xx[k] + gg[k] - result2;
        
        
        
        // printf("teste = %f - %f - %f - %f - %f -%f - %f - %f - %f \n", lambda0, xx[k + 1], xx[k], gg[k], result2, Xk[0], Gk[0], gammaK, Sbbp);
        if (fabs(lambda - lambdaAnt) <= PREC){
            mount_x_a(p,x, lambda,ur,lr);
            break;
        }
        
        gg[k + 1] = lambda - xx[k + 1];
        
        Xk[k] = xx[k + 1] - xx[k];
        Gk[k] = gg[k + 1] - gg[k];

        if (k >= m){
            for (unsigned i =0; i <= m; i++){
                Xk[i] = Xk[i+1];
                Gk[i] = Gk[i+1];
            }
        }

        // Xk = xx[k + 1] - xx[k];
        // Gk = gg[k + 1] - gg[k];

        lambdaAnt = lambda;
        k++;
    }

    free(abp);
    free(bbp);
    free(lr);
    free(ur);
    free(bu);

    return k+1;
}
