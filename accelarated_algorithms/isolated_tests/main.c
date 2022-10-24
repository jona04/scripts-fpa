#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// #include "../lib/polinomial.c"
#include "../../state_of_the_art_algorithms/lib/fixed_point.c"
#include "../lib/fixed_point_ai.c"
#include "../lib/fixed_point_a.c"
#include "../../state_of_the_art_algorithms/lib/fixed_point_new.c"
#include "../../state_of_the_art_algorithms/lib/cont_quad_knapsack.c"

// #include "../lib/cont_quad_knapsack.h"
// #include "../lib/third_party_methods.h"
// #include "../../third_party/dSFMT-src-2.2.2/dSFMT.h"

/* Possible benchmark types */
typedef enum {
    correlated,
    weakly_correlated,
    uncorrelated,
    flow,
    invalid
} problem_type;

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
// int RAND_MAX = 32767;
double r2()
{
    return (double)rand() / (double)RAND_MAX ;
}

void printRandoms(int lower, int upper,  
                             int count) 
{ 
    int i; 
    for (i = 0; i < count; i++) { 
        double num = (rand() % 
           (upper - lower + 1)) + lower; 
        // printf("%f ", num); 
    } 
} 

/* Global random number generator instance */
// static dsfmt_t dsfmt;

// double random_uniform () {
//     return dsfmt_genrand_close_open(&dsfmt);;
// }

double generatorNumber(int lower, int upper) 
{ 
    return (double)(rand() % (upper - lower + 1)) + lower; 
} 

// void residual_dist(cqk_problem *restrict p, double *restrict x,
//                    double *restrict residual, double *restrict dist) {
//     *dist = 0.0;
//     *residual = 0.0;
//     for (unsigned i = 0; i < p->n; ++i) {
//         *dist += (x[i] - p->a[i]/p->d[i])*(x[i] - p->a[i]/p->d[i]);
//         *residual += p->b[i]*x[i];
//     }
//     *residual = (*residual - p->r)/fmax(1.0, p->r);
// }

void limit(cqk_problem *restrict p, double *restrict x,
                   double *restrict limit, double *obj) {
    *limit = 0.0;
    *obj = 0.0;
    for (unsigned i = 0; i < p->n; ++i) {
        *limit += p->b[i]*x[i];
        *obj += p->d[i]*(x[i]*x[i]) - p->a[i]*x[i];
    }
}

void generate_cqk_problem(problem_type type, cqk_problem *restrict p) {
    double
        temp,
        lb = 0.0, /* Receives p->b^t p->low. */
        ub = 0.0; /* Receives p->b^t p->up. */
    int lower = 10, upper = 25; 
    /* Generate the problem */

    // Use current time as  
    // seed for random generator 
    srand(time(0)); 

    for (unsigned i = 0; i < p->n; ++i) {


        p->b[i] = 10 + 14*r2();
        p->low[i] = 1 + 14*r2();
        p->up[i] = 1 + 14*r2();
        if (p->low[i] > p->up[i]) {
            temp = p->low[i];
            p->low[i] = p->up[i];
            p->up[i] = temp;
        }
        lb += p->b[i]*p->low[i];
        ub += p->b[i]*p->up[i];
        switch (type) {

        case uncorrelated:
            p->d[i] = generatorNumber(10,25);
            p->a[i] = generatorNumber(10,25);
            break;

        case weakly_correlated:
            p->d[i] = p->b[i] - generatorNumber(10,25);
            p->a[i] = p->b[i] - generatorNumber(10,25);
            break;

        case correlated:
        default:
            p->d[i] = p->b[i] + 5.0;
            p->a[i] = (p->b[i] + 5.0);
        }
    }
    p->r = lb + (ub - lb)*0.7;
}

// void generate_cqk_problem(problem_type type, cqk_problem *restrict p) {
//     double
//         temp,
//         lb = 0.0, /* Receives p->b^t p->low. */
//         ub = 0.0; /* Receives p->b^t p->up. */

//     /* Generate the problem */
//     for (unsigned i = 0; i < p->n; ++i) {
//         p->b[i] = 10.0 + 15.0*random_uniform();
//         p->low[i] = 1.0 + 14.0*random_uniform();
//         p->up[i] = 1.0 + 14.0*random_uniform();
//         if (p->low[i] > p->up[i]) {
//             temp = p->low[i];
//             p->low[i] = p->up[i];
//             p->up[i] = temp;
//         }
//         lb += p->b[i]*p->low[i];
//         ub += p->b[i]*p->up[i];
//         switch (type) {

//         case uncorrelated:
//             p->d[i] = 10.0 + 15.0*random_uniform();
//             p->a[i] = 10.0 + 15.0*random_uniform();
//             break;

//         case weakly_correlated:
//             p->d[i] = p->b[i] - 5.0 + 10.0*random_uniform();
//             p->a[i] = p->b[i] - 5.0 + 10.0*random_uniform();
//             break;

//         case correlated:
//         default:
//             p->d[i] = p->b[i] + 5.0;
//             p->a[i] = p->b[i] + 5.0;
//         }
//     }
//     p->r = lb + (ub - lb)*random_uniform();
// }


void result_vars(cqk_problem *restrict p, double *restrict x) {
    for(unsigned i = 0; i < p->n; i++) {
        printf( " %f ", p->low[i] );
    }
    printf("\n");
    for(unsigned i = 0; i < p->n; i++) {
        printf( " %f ", p->up[i] );
    }
    printf("\n");
    for(unsigned i = 0; i < p->n; i++) {
        printf( " %f ", x[i] );
    }
    printf("\n");
}


void gerar(cqk_problem *restrict p){
    printf("d = %f - %f - %f - %f - %f \n", p->d[0], p->d[1], p->d[2], p->d[3], p->d[4]);
    printf("a = %f - %f - %f - %f - %f \n", p->a[0], p->a[1], p->a[2], p->a[3], p->a[4]);
    printf("b = %f - %f - %f - %f - %f \n", p->b[0], p->b[1], p->b[2], p->b[3], p->b[4]);
    printf("up = %f - %f - %f - %f - %f \n", p->up[0], p->up[1], p->up[2], p->up[3], p->up[4]);
    printf("low = %f - %f - %f - %f - %f \n", p->low[0], p->low[1], p->low[2], p->low[3], p->low[4]);
    printf("r = %f \n", p->r);
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
void run_a_test1(cqk_problem *restrict p, double *restrict x) {
    /* Run the test and measure its running time. */
    double status = -1;

    status = fixed_point(p,x);
    printf("\nresult 1 = %f\n",status);
    double limit_var, obj;
    limit(p, x, &limit_var, &obj);
    // result_vars(p,x);
    // printf("x = %f - %f - %f \n", x[0], x[1], p->d[0]);
    printf("limit1 = %f | obj = %f \n", limit_var, obj);
} 

void run_a_test11(cqk_problem *restrict p, double *restrict x) {
    /* Run the test and measure its running time. */
    double status = -1;

    status = fixed_point_new(p,x);
    printf("\nresult 2 = %f\n",status);
    double limit_var, obj;
    limit(p, x, &limit_var, &obj);
    // result_vars(p,x);
    // printf("x = %f - %f \n", x[0], x[1]);
    printf("limit2 = %f | obj = %f \n", limit_var, obj);
} 

void run_a_test2(cqk_problem *restrict p, double *restrict x) {
    /* Run the test and measure its running time. */
    double status = -1;

    status = newton(p, NULL, x);
    printf("\nresult 1 = %f\n",status);
    double limit_var, obj;
    limit(p, x, &limit_var, &obj);
    // result_vars(p,x);
    printf("\n limit = %f \n", limit_var);
} 

void run_a_test3(cqk_problem *restrict p, double *restrict x) {
    /* Run the test and measure its running time. */
    double status = -1;

    status = fixed_point_ai(p, x);
    printf("\nresult 1 = %f\n",status);
    double limit_var, obj;
    limit(p, x, &limit_var, &obj);
    // result_vars(p,x);
    printf("\n limit = %f \n", limit_var);
    printf("limit3 = %f | obj = %f \n", limit_var, obj);
} 

void run_a_test4(cqk_problem *restrict p, double *restrict x) {
    /* Run the test and measure its running time. */
    double status = -1;

    status = fixed_point_a(p, x);
    printf("\nresult 1 = %f\n",status);
    double limit_var, obj;
    limit(p, x, &limit_var, &obj);
    // result_vars(p,x);
    printf("\n limit = %f \n", limit_var);
    printf("limit4 = %f | obj = %f \n", limit_var, obj);
} 

int main(int argc, char *argv[]) {

    clock_t start1, end1;
    double cpu_time_used1;

    problem_type ptype; /* Code for the problem type */
    unsigned n;         /* Problem dimension */ 
    n = 2e6;
    ptype = uncorrelated;

    /* Problem data */
    cqk_problem p;
    allocate_cqk_problem(n, &p);

    /* Solution vector */
    double *x = (double *) malloc(p.n*sizeof(double));

    generate_cqk_problem(ptype, &p);

    // p.d[0] = 17.276016; p.d[1] = 20.366308; p.d[2] = 11.126102; p.d[3] = 13.357732; p.d[4] = 21.254665; 
    // p.a[0] = 16.100559; p.a[1] = 24.946266; p.a[2] = 18.833744; p.a[3] = 16.180804; p.a[4] = 11.500267; 
    // p.b[0] = 14.552062; p.b[1] = 18.132834; p.b[2] = 10.820949; p.b[3] = 17.423883; p.b[4] = 24.874358; 
    // p.up[0] = 8.590115; p.up[1] = 8.170043; p.up[2] = 12.415994; p.up[3] = 10.676436; p.up[4] = 6.141831; 
    // p.low[0] = 3.696003; p.low[1] = 4.329827; p.low[2] = 9.009884; p.low[3] = 8.129925; p.low[4] = 1.120140;
    // p.r = 647.768350;
    // p.n = 5;


    // gerar(&p);

    start1 = clock();
    run_a_test1(&p,x);
    end1 = clock();
    cpu_time_used1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
    printf("fun() took %f seconds to execute \n", cpu_time_used1);

    // start1 = clock();
    // run_a_test11(&p,x);
    // end1 = clock();
    // cpu_time_used1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
    // printf("fun() took %f seconds to execute \n", cpu_time_used1);


    // start1 = clock();
    // run_a_test2(&p,x);
    // end1 = clock();
    // cpu_time_used1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
    // printf("fun() took %f seconds to execute \n", cpu_time_used1);

    start1 = clock();
    run_a_test3(&p,x);
    end1 = clock();
    cpu_time_used1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
    printf("fun() took %f seconds to execute \n", cpu_time_used1);

    start1 = clock();
    run_a_test4(&p,x);
    end1 = clock();
    cpu_time_used1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
    printf("fun() took %f seconds to execute \n", cpu_time_used1);

    free_cqk_problem(&p);

    return 0;

}