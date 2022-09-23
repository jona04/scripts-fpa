#ifndef FIXED_POINT_H
#define FIXED_POINT_H

#include "cont_quad_knapsack.h"

double initial_multiplier_fixed_point(cqk_problem *restrict p,
                                double *restrict abp,
                                double *restrict bbp,
                                double *restrict lr,
                                double *restrict ur,
                                double *restrict bl,
                                double *restrict bu,
                                double *restrict alpp);

int fixed_point(cqk_problem *restrict prob, double *restrict x);

#endif
