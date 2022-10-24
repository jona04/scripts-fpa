#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h"

void findIntervalS(cqk_problem *restrict p, double *intervalA, double *intervalB,
                    double *ga, double *gb){

    for (double k = 0; k < 100; k++){

        double xb = 0;
        double xa = 0;
        *ga = 0;
        *gb = 0;

        for (int i = 0; i < p->n; i++){
            
            xa = (p->a[i] - (-k)*p->b[i] ) / p->d[i] ;
            if(xa < p->low[i]){
                xa = p->low[i];
            }else{
                if(xa > p->up[i]){
                    xa = p->up[i];
                }
            }
            *ga += p->b[i]*xa;

            xb = (p->a[i] - (k)*p->b[i] ) / p->d[i] ;
            if(xb <= p->low[i]){
                xb = p->low[i];
            }else{
                if(xb > p->up[i]){
                    xb = p->up[i];
                }
            }
            *gb += p->b[i]*xb;

        }
        *ga = *ga - p->r;
        *gb = *gb - p->r;
        if ((*ga)*(*gb) < 0.00000){
            *intervalA = -k;
            *intervalB = k;
            return;
        }
    }

}

int secant_method(cqk_problem *restrict p, double *x, double eps){

    double intervalA;
    double intervalB;
    double ga;
    double gb;

    findIntervalS(p, &intervalA, &intervalB, &ga, &gb);
    double l0 = intervalA;
    double l1 = intervalB;
    double l2 = 0;

    int k_aux = 0;
    for (int k = 0; k < 100; k++){

        double gl0 = 0;
        double gl1 = 0;
        double gl2 = 0;
        double xl0 = 0;
        double xl1 = 0;
        double xl2 = 0;
        
        for (int i = 0; i < p->n; i++){
            xl0 = (p->a[i] - (l0)*p->b[i] ) / p->d[i] ;
            if(xl0 <= p->low[i]){
                xl0 = p->low[i];
            }else{
                if(xl0 > p->up[i]){
                    xl0 = p->up[i];
                }
            }
            gl0 += p->b[i]*xl0;

            xl1 = (p->a[i] - (l1)*p->b[i] ) / p->d[i] ;
            if(xl1 < p->low[i]){
                xl1 = p->low[i];
            }else{
                if(xl1 > p->up[i]){
                    xl1 = p->up[i];
                }
            }
            gl1 += p->b[i]*xl1;

        }

        gl0 = gl0 - p->r;
        gl1 = gl1 - p->r;
        l2 = l1 - ((l1 - l0)/(gl1 - gl0))*gl1;
        
        for (int i = 0; i < p->n; i++){
            xl2 = (p->a[i] - (l2)*p->b[i] ) / p->d[i] ;
            if(xl2 <= p->low[i]){
                xl2 = p->low[i];
            }else{
                if(xl2 > p->up[i]){
                    xl2 = p->up[i];
                }
            }
            gl2 += p->b[i]*xl2;
        }
        gl2 = gl2 - p->r;
        if(fabs(gl2) < eps){
            return k + k_aux;
        }
        l0 = l1;
        l1 = l2;
    }
    return -1;
}
