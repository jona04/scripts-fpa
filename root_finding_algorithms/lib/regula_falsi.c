#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h"

void findIntervalR(cqk_problem *restrict p, double *intervalA, double *intervalB,
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

int regula_falsi(cqk_problem *restrict p, double *x, double eps){

    double intervalA;
    double intervalB;
    double gaa;
    double gbb;

    findIntervalR(p, &intervalA, &intervalB, &gaa, &gbb);
    for (int k = 0; k < 100; k++){

        double ga = 0;
        double gb = 0;
        double gl = 0;
        double xa = 0;
        double xb = 0;
        double xl = 0;
        double lambda = 0;
        
        for (int i = 0; i < p->n; i++){
            xa = (p->a[i] - intervalA*p->b[i] ) / p->d[i] ;
            if(xa <= p->low[i]){
                xa = p->low[i];
            }else{
                if(xa > p->up[i]){
                    xa = p->up[i];
                }
            }
            ga += p->b[i]*xa;

            xb = (p->a[i] - intervalB*p->b[i] ) / p->d[i] ;
            if(xb <= p->low[i]){
                xb = p->low[i];
            }else{
                if(xb > p->up[i]){
                    xb = p->up[i];
                }
            }
            gb += p->b[i]*xb;
        }
        ga = ga - p->r;
        gb = gb - p->r;

        lambda = ( (intervalA * gb) - (intervalB * ga) )/(gb - ga);

        for (int i = 0; i < p->n; i++){
            xl = (p->a[i] - lambda*p->b[i] ) / p->d[i] ;
            if(xl <= p->low[i]){
                xl = p->low[i];
            }else{
                if(xl > p->up[i]){
                    xl = p->up[i];
                }
            }
            gl += p->b[i]*xl;
        }
        gl = gl - p->r;
        if (fabs(gl) < eps){
            return k;
        }
        if ( (gl>0.0000000 && ga>0.0000000) || (gl<0.0000000 && ga<0.0000000) ){
            intervalA = lambda;
        }
        else{
            intervalB = lambda;
        }
    }

    return -1;
}
