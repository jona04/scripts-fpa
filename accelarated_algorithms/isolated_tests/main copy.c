#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// #include "../lib/polinomial.c"
#include "../lib/softmax.c"
#include "../lib/fixed_point.c"
#include "../lib/qwmca.c"
#include "../lib/qwmca2.c"
#include "../lib/cont_quad_knapsack.c"


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

double generatorNumber(int lower, int upper) 
{ 
    return (double)(rand() % (upper - lower + 1)) + lower; 
} 

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
void run_a_test(cqk_problem *restrict p, double *restrict x) {
    /* Run the test and measure its running time. */
    double status = -1;

    double delta_lambda = -1.0;

    // status = newton(p, NULL, x);

    status = fixed_point(p,x);
    printf("\nresult 1 = %f - %f\n",status,x[1]);
    double residual, dist;
    residual_dist(p, x, &residual, &dist);
    printf("\n res dist = %f -  %f \n", residual, dist);

    // status = qwmca(p,x);
    // status = qwmca2(p,x);
     status = softmax(p,x);
    printf("\nresult 1 = %f - %f\n",status,x[1]);
    residual_dist(p, x, &residual, &dist);
    printf("\n res dist = %f -  %f \n", residual, dist);

    // status = polinomial(p,x);
    // printf("\nresult 1 = %f - %f\n",status,x[1]);
    // residual_dist(p, x, &residual, &dist);
    // printf("\n res dist = %f -  %f \n", residual, dist);
} 



int main(int argc, char *argv[]) {

    clock_t start1, end1;
    double cpu_time_used1;

    problem_type ptype; /* Code for the problem type */
    unsigned n;         /* Problem dimension */ 
    n = 2e6;
    ptype = correlated;

    start1 = clock();
    // for (int i=0; i < n; i++){
    //     fasterlog(1);
    // }
    // double rr=0;
    // for (int i=0; i < n; i++){
    //     rr = fasterlog(1);
    // }
    end1 = clock();
    cpu_time_used1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
    printf("fun() took %f seconds to execute \n", cpu_time_used1);

    /* Problem data */
    cqk_problem p;
    allocate_cqk_problem(n, &p);

    /* Solution vector */
    double *x = (double *) malloc(p.n*sizeof(double));

    generate_cqk_problem(ptype, &p);


    start1 = clock();
    
    run_a_test(&p,x);

    end1 = clock();
    cpu_time_used1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
    printf("fun() took %f seconds to execute \n", cpu_time_used1);

    free_cqk_problem(&p);
    return 0;

}