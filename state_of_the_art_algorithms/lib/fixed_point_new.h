#ifndef FIXED_POINT_NEW_H
#define FIXED_POINT_NEW_H

#include "cont_quad_knapsack.h"

double initial_multiplier_fixed_point_new(cqk_problem *restrict p, 
                                double *restrict r, 
                                double *restrict a,
                                double *restrict u,
                                double *restrict abp,
                                double *restrict bbp,
                                double *restrict ur,
                                double *restrict bl,
                                double *restrict bu);

int fixed_point_new(cqk_problem *restrict prob, double *restrict x);

#endif
