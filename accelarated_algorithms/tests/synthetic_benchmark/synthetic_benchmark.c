/*
 * Module Name: synthetic_benchmark
 *
 * Description: Implements a set of synthetically generated random
 *     tests for methods to solve continuous quadratic knapsack
 *     problems. The first group of tests are described in Kiwiel's
 *     papers, like for example
 *
 *     Kiwiel KC. Variable Fixing Algorithms for the Continuous
 *     Quadratic Knapsack Problem. Journal of Optimization Theory and
 *     Applications. 2008;136(3):445-458. Available at:
 *     http://www.springerlink.com/index/10.1007/s10957-007-9317-7.
 *
 *     We also implement a test described by Dai and Fletcher that
 *     should capture the structure of some real world multicommodity
 *     network flow and logistics in the paper
 *
 *     Dai Y-H, Fletcher R. New algorithms for singly linearly
 *     constrained quadratic programs subject to lower and upper
 *     bounds. Mathematical Programming. 2005;106(3):403-421. Available
 *     at: http://www.springerlink.com/index/10.1007/s10107-005-0595-2.
 *
 * Copyright: Paulo J. S. Silva <pjssilva@gmail.com> 2012.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../../../state_of_the_art_algorithms/lib/fixed_point.h"
#include "../../../state_of_the_art_algorithms/lib/fixed_point_new.h"
#include "../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h"
#include "../../../state_of_the_art_algorithms/lib/third_party_methods.h"
#include "../../lib/fixed_point_a.h"
#include "../../lib/fixed_point_ai.h"
#include "../../../third_party/dSFMT-src-2.2.2/dSFMT.h"

/* Possible benchmark types */
typedef enum {
    correlated,
    weakly_correlated,
    uncorrelated,
    flow,
    invalid
} problem_type;

/* Global random number generator instance */
static dsfmt_t dsfmt;

/*
 * Function Name: random_uniform
 *
 * Description: Naive random number generator with uniform distribution.
 *
 * Return value: a pseudo random number uniformly in [0, 1)
 */
double random_uniform () {
    return dsfmt_genrand_close_open(&dsfmt);;
}

/*
 * Function Name: generate_cqk_problem
 *
 * Description: Create a test instace for the projection problems
 *     following the description in Kiwiel papers.
 *
 * Input:
 *     problem_type type: problem category (uncorrelated, weakly
 *         correlated, correlated) Input/Output
 *
 * Input/Output:
 *     cqk_problem *restrict p: on input p->n has the dimension of the
 *        desired problem. At output it points to a new random
 *        problem.
 */
void generate_cqk_problem(problem_type type, cqk_problem *restrict p) {
    double
        temp,
        lb = 0.0, /* Receives p->b^t p->low. */
        ub = 0.0; /* Receives p->b^t p->up. */

    /* Generate the problem */
    for (unsigned i = 0; i < p->n; ++i) {
        p->b[i] = 10.0 + 15.0*random_uniform();
        p->low[i] = 1.0 + 14.0*random_uniform();
        p->up[i] = 1.0 + 14.0*random_uniform();
        if (p->low[i] > p->up[i]) {
            temp = p->low[i];
            p->low[i] = p->up[i];
            p->up[i] = temp;
        }
        lb += p->b[i]*p->low[i];
        ub += p->b[i]*p->up[i];
        switch (type) {

        case uncorrelated:
            p->d[i] = 10.0 + 15.0*random_uniform();
            p->a[i] = 10.0 + 15.0*random_uniform();
            break;

        case weakly_correlated:
            p->d[i] = p->b[i] - 5.0 + 10.0*random_uniform();
            p->a[i] = p->b[i] - 5.0 + 10.0*random_uniform();
            break;

        case correlated:
        default:
            p->d[i] = p->b[i] + 5.0;
            p->a[i] = p->b[i] + 5.0;
        }
    }
    p->r = lb + (ub - lb)*random_uniform();
}

/*
 * Function Name: generate_flow_problem
 *
 * Description: Generate a random instance related to multicommodity
 *     network flow and logistics problems, as suggested in Dai and
 *     Fletcher paper.
 *
 * Input/Output:
 *     cqk_problem *restrict p: on input p->n has the dimension of the
 *         desired problem. At output it points to a new random
 *         problem.
 */
void generate_flow_problem(cqk_problem *restrict p) {
    double
        lb = 0.0, /* Receives p->b^t p->low. */
        ub = 0.0; /* Receives p->b^t p->up. */

    for (unsigned i = 0; i < p->n; ++i) {
        p->d[i] = 1.0 + 1.0e+4*random_uniform();
        p->a[i] = -1000.0 + 2000.0*random_uniform();
        p->low[i] = 0.0;
        p->up[i] = 1000*random_uniform();
        p->b[i] = 1.0;
        lb += p->b[i]*p->low[i];
        ub += p->b[i]*p->up[i];
    }
    p->d[0] = 1.0;
    p->d[p->n - 1] = 1.0e+4;
    p->r = lb + (ub - lb)*random_uniform();
}

/*
 * Function Name: residual_dist
 *
 * Description: Given a point that conform to the problem bounds,
 *     compute the relative residual of the linear contraint to see if
 *     it is really feasible. Also compute the D-distance to D^{-1}a,
 *     the point that should be projected on the feasible set when
 *     viewing the problema as a projection problem in D-norm. This
 *     values may be used to assess the quality of an approximate
 *     solution, for example comparing different methods to see if
 *     they compute the same point.
 *
 * Input:
 *     cqk_problem *restrict p: cont. quad. knapsack problem.
 *     double *restrict x: approximate solution vector.
 *
 * Output:
 *     double *retrict residual: relative residual of the linear constraint.
 *     double *restrict dist: D-distance to D^{-1}a
 */
void residual_dist(cqk_problem *restrict p, double *restrict x,
                   double *restrict residual, double *restrict dist) {
    *dist = 0.0;
    *residual = 0.0;
    for (unsigned i = 0; i < p->n; ++i) {
        *dist += (x[i] - p->a[i]/p->d[i])*(x[i] - p->a[i]/p->d[i]);
        *residual += p->b[i]*x[i];
    }
    *residual = (*residual - p->r)/fmax(1.0, p->r);
}

/*
 * Function Name: get_problem_type
 *
 * Description: Maps a string describing the problem type its code.
 *
 * Input:
 *     char *restrict option: string description of a problem type
 *
 * Return value: code of the problem.
 */
problem_type get_problem_type(char *restrict option) {
    if (strcmp(option, "correlated") == 0)
        return correlated;
    else if (strcmp(option, "weakly_correlated") == 0)
        return weakly_correlated;
    else if (strcmp(option, "uncorrelated") == 0)
        return uncorrelated;
    else if (strcmp(option, "flow") == 0)
        return flow;
    else
        return invalid;
}

/*
 * Function Name: run_a_test
 *
 * Description: Run a solver for a CQN problem and save performance
 *     statistics.
 *
 * Input:
 *     cqk_problem *restrict p: problem description.
 *     char *restrict method: method name (newton, secant, fix, or search).
 *     unsigned nreps: number of repetitions to get reasonable timing
 *          statistics.
 *     FILE *out: file to save statistics.
 *
 * Output:
 *     double *restrict x: computed solution vector.
 *
 * Side effects: The file out will receive the statistics.
 */
void run_a_test(cqk_problem *restrict p, char *restrict method, unsigned nreps,
                FILE *restrict out, double *restrict x) {
    /* Run the test and measure its running time. */
    int status = -1;
    double start, duration;
    printf("Start %s...\n", method);
    if (!strcmp(method, "fixed_point")) {
        start = clock();
        for (unsigned i = 0; i < nreps; ++i){
            status = fixed_point(p, x);
        }
        duration = (double) clock() - start;
    }else if (!strcmp(method, "fixed_point_new")) {
        start = clock();
        for (unsigned i = 0; i < nreps; ++i){
            status = fixed_point_new(p, x);
        }
        duration = (double) clock() - start;
    }else if (!strcmp(method, "aitken")) {
        start = clock();
        for (unsigned i = 0; i < nreps; ++i){
            status = fixed_point_ai(p, x);
        }
        duration = (double) clock() - start;
    }else if (!strcmp(method, "anderson")) {
        start = clock();
        for (unsigned i = 0; i < nreps; ++i){
            status = fixed_point_a(p, x);
        }
        duration = (double) clock() - start;
    }
    duration /= CLOCKS_PER_SEC;

    /* Print and save statistics */
    printf("Projection time=%.4e %d\n\n", duration/nreps, status);
    double residual, dist;
    if (status >= 0) {
        residual_dist(p, x, &residual, &dist);
    } else {
        dist = residual = -1.0;
        duration = -duration;
    }
    fprintf(out, "%d\t%e\t%e\t%e\t", status, duration/nreps, dist, residual);
}

/*
 * Function Name: main
 *
 * Description: Main function that runs a bunch of random test of
 *     given type.
 *
 * Usage: synthetic <test type> <problem dimension>
 *
 *        where test type is in {correlated, weakly_correlated,
 *        uncorrelated, flow}
 *
 * Input:
 *     int argc: number of command line paramters.
 *     char *argv[]: command line paramters.
 */
int main(int argc, char *argv[]) {
    /* List of possible methods */
    char *methods[4] = {"aitken","anderson","fixed_point_new","fixed_point"};
    const unsigned n_methods = 4;

    problem_type ptype; /* Code for the problem type */
    unsigned n;         /* Problem dimension */

    /* Get the problem dimensions and define number of test
       repetitions. */
    if (argc != 3 || (n = atoi(argv[2])) == 0 ||
        ((ptype = get_problem_type(argv[1])) == invalid)) {
        fprintf(stderr, "Usage: synthetic <test type> <problem dimension>\n");
        exit(1);
    }
    printf("Projection type:%d\n", ptype);

    /* Total number of random tests to do. */
    unsigned ntests = 5;
    /* Number of repetions for each test */
    unsigned nreps = 3;
    nreps = nreps > 0 ? nreps : 1;

    /* Problem data */
    cqk_problem p;
    allocate_cqk_problem(n, &p);

    /* Solution vector */
    double *x = (double *) malloc(p.n*sizeof(double));

    /* Initialize the random number generator with known seed so that
       the same 'random' problems are generated always. */
    dsfmt_init_gen_rand(&dsfmt, 12345);

    /* File to hold the statistics */
    FILE *out = fopen("times.dat", "w");

    for (unsigned repetitions = 0; repetitions < ntests; ++repetitions) {

        if (ptype != flow){
            generate_cqk_problem(ptype, &p);
        }
        else{
            generate_flow_problem(&p);
        }
        for (unsigned met = 0; met < n_methods; ++met)
            run_a_test(&p, methods[met], nreps, out, x);
        fprintf(out, "\n");
    }
    free_cqk_problem(&p);
    return 0;
}
