/*
 * Module Name: synthetic_benchmark
 *
 * Description: Implements a set of sythentically generated random
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
#ifndef SYNT_BENCH_H
#define SYNT_BENCH_H
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
void generate_cqk_problem(problem_type type, cqk_problem *restrict p);

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
void generate_flow_problem(cqk_problem *restrict p);

#endif
