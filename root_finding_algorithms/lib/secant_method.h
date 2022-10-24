#ifndef SECANT_METHOD_H
#define SECANT_METHOD_H

#include "../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h"

#include <float.h>

void findIntervalS(cqk_problem *restrict p, double *intervalA, double *intervalB,
                    double *ga, double *gb);

int secant_method(cqk_problem *restrict p, double *x, double eps);


#endif
