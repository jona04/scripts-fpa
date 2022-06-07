#ifndef REGULA_FALSI_H
#define REGULA_FALSI_H

#include "bissection_method.h"

#include <float.h>

int regula_falsi(cqk_problem *restrict p, double *x , double eps);


#endif
