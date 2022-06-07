#ifndef FIXED_POINT_H
#define FIXED_POINT_H

#include "fixed_point.h"

#include <float.h>

/* 
 * Structure to describe a continuous quadratic knapsack problem in the form
 * 
 * min  1/2 x'Dx - a'x
 * s.t. b'x = r
 *      low <= x <= up,
 *
 * where D is a positive diagonal matrix.
 *
 * IMPORTANT: As in the manuscript the code assumes that b > 0. It is
 * the user responsibility to change the problem to fit this format.
 */
typedef struct cqk_problem {
    unsigned n;           /* Dimension of the problem. */
    double *restrict d;   /* D positive diagonal (D = diag(d)). */
    double *restrict a;   /* The a in the objective funtion definition. */
                          /* If the problem is interpreted as projection */
                          /* in the norm induced by D, the point being */  
                          /* projected is D^-1a. */
    double *restrict b;   /* Slopes (positive) that define the plane. */
    double r;             /* Righthand side of the plane equation. */
    double *restrict low; /* lower and */
    double *restrict up;  /* upper bounds for the constraints. */
} cqk_problem;

/********** Interface of a cqn_structure *********/

/*
 * Function Name: allocate_cqk_problem
 *
 * Description: Allocate the memory necessary hold a CQN problem in
 *     dimension n.
 *
 * Input:
 *     unsigned n: The dimension of the problem.
 * 
 * Output: cqk_problem *p: pointer to the cqk_problem whose memory is
 *     being allocated.
 */
void allocate_cqk_problem(unsigned n, cqk_problem *restrict p);

/*
 * Function Name: free_cqk_problem
 *
 * Description: Free the memory associated to a cqk_problem.
 *
 * Input/Output: 
 *     cqk_problem *p: at entry a cqk_problem with allocated
 *         memory, at output all the memory will be freed and the pointers
 *         in the cqk_problem structure will point to NULL.
 */
void free_cqk_problem(cqk_problem *restrict p);

/* Maximum number of iterations for each method. */
#define MAXITERS 100000

double initial_multiplier_fixed_point(cqk_problem *restrict p, 
                                double *restrict abp,
                                double *restrict bbp,
                                double *restrict lr,
                                double *restrict ur,
                                double *restrict bl,
                                double *restrict bu);    

int fixed_point(cqk_problem *restrict prob, double *restrict x);

#endif
