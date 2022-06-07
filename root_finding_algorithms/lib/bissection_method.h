#ifndef BISSECTION_METHOD_H
#define BISSECTION_METHOD_H

#include "fixed_point.h"

#include <float.h>

void findInterval(cqk_problem *restrict p, double *intervalA, double *intervalB,
                    double *ga, double *gb);

int bissection_method(cqk_problem *restrict p, double *x ,double eps);

#endif
