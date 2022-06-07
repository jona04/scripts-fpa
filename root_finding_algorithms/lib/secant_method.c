#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "secant_method.h"


int secant_method(cqk_problem *restrict p, double *x, double eps){

    double intervalA;
    double intervalB;
    double ga;
    double gb;

    findInterval(p, &intervalA, &intervalB, &ga, &gb);
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
