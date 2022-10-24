#ifndef FIXED_POINT_AI_H
#define FIXED_POINT_AI_H

#include "../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h"

double initial_multiplier_fixed_point_ai(cqk_problem *restrict p, 
                                double *restrict r, 
                                double *restrict a,
                                double *restrict u,
                                double *restrict abp,
                                double *restrict bbp,
                                double *restrict ur,
                                double *restrict bl,
                                double *restrict bu);

int fixed_point_ai(cqk_problem *restrict prob, double *restrict x);

#endif
