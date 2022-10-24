#ifndef REGULA_FALSI_H
#define REGULA_FALSI_H

#include "../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h"

#include <float.h>

void findIntervalR(cqk_problem *restrict p, double *intervalA, double *intervalB,
                    double *ga, double *gb);

int regula_falsi(cqk_problem *restrict p, double *x , double eps);


#endif
