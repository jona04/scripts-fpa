#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/* include the user-visible definitions for BLG.h */
#include "BLG_user.h"

/* ========================================================================== */
/* === macros =============================================================== */
/* ========================================================================== */

#define PRIVATE static
#define BLGZERO ((BLGFLOAT) 0)
#define BLGONE ((BLGFLOAT) 1)
#define BLGTWO ((BLGFLOAT) 2)

#define BLGMIN(a,b) (((a) < (b)) ? (a) : (b))
#define BLGMAX(a,b) (((a) > (b)) ? (a) : (b))

typedef struct BLGcom_struct /* common variables */
{
    BLGobjective *user ; /* information passed to user when function or
                            gradient must be evaluated */
    BLGINT          nf ; /* total number of function evaluations */
    BLGINT          ng ; /* total number of gradient evaluations */
    BLGINT         nmu ; /* number of iterations in mu evaluation for AS */
    BLGFLOAT        *a ; /* a'x = b (linear constraint vector) */
    BLGFLOAT         b ; /* right side */
    BLGFLOAT       *lo ; /* lower bounds */
    BLGFLOAT       *hi ; /* upper bounds */
    BLGFLOAT       tol ; /* convergence tolerance */
    int              n ; /* problem dimension */
    int           nsub ; /* nominal dimension of subspace in subspace mode */
    BLGFLOAT   bndnorm ; /* absolute bound on the elements of lo and hi */
    BLGFLOAT     anorm ; /* infinity norm of a vector */
    BLGFLOAT     a2sum ; /* a'a in the subspace for nonzero d_i */
    int           *iw1 ; /* work array of size n */
    int           *iw2 ; /* work array of size n */
    int           *iw3 ; /* work array of size n (n+1 in subspace mode) */
    int         *ikeep ; /* array of size n initialized to 0 (subspace only) */
    BLGFLOAT       *w1 ; /* work array of size n */
    BLGFLOAT       *w2 ; /* work array of size n */
    BLGFLOAT       *w3 ; /* work array of size n */
    BLGparm      *Parm ; /* parameter structure */
    void    (*value) (BLGobjective *) ; /* evaluate objective function */
    void     (*grad) (BLGobjective *) ; /* evaluate objective gradient */
    void  (*valgrad) (BLGobjective *) ; /* function & gradient if given*/
    /* if objective function is quadratic (Parm->QuadCost = TRUE) */
    void       (*Ad) (BLGobjective *) ; /* compute A*d if cost quadratic */
    BLGFLOAT    deriv2 ; /* d'Ad */
} BLGcom ;

/* prototypes */

PRIVATE void BLGevaluate
(
    BLGcom *Com,  /* BLGcom structure */
    char  *what   /* f = function, g = gradient, fg = function and gradient */
) ;


PRIVATE int BLGbuildSubspace /* return:
                     0 (error tolerance satisfied)
                     1 (refine accuracy for prior subspace)
                     2 (subspace has been evaluated)
                    10 (null subspace, d_pert too big) */
(
    int            *i, /* indices of previous subspace (input)
                        indices of chosen subspace (returned) */
    int         *nsub, /* dimension of previous subspace (input)
                          size of subspace (returned) */
    BLGFLOAT     *Err, /* error estimate (returned) */
    BLGFLOAT      *nu, /* multiplier estimate (input and returned) */
    BLGFLOAT   SubErr, /* error for previous subspace */
    BLGFLOAT       *x, /* current estimate of solution to problem */
    BLGFLOAT       *g, /* gradient at current point */
    BLGFLOAT lambda_k, /* current BB parameter */
    BLGcom       *Com  /* SVM Com structure */
) ;

PRIVATE int BLGd /* return:
                   0 (search direction was computed)
                   4 (computation of AS direction fails when bracketing)
                   6 (too many iterations when computing AS direction) */
(
    BLGFLOAT       *d, /* search direction (returned) */
    BLGFLOAT      *mu, /* input guess for mu, return final mu */
    BLGFLOAT    *xnew, /* xnew = x + d (returned) */
    BLGFLOAT     *Gtd, /* g'd (returned) */
    BLGFLOAT     *Err, /* KKT error (returned) */
    BLGFLOAT       *x, /* current iterate */
    BLGFLOAT       *g, /* gradient at x */
    BLGFLOAT       *a, /* linear constraint is a'x = b */
    BLGFLOAT      *lo, /* lower bounds */
    BLGFLOAT      *hi, /* upper bounds */
    int             n, /* dimension */
    BLGFLOAT lambda_k, /* constant in denominator of d */
    int        use_gp, /* use_gp = TRUE means use gradient project */
    int        use_lp, /* use_lp = TRUE means must use Frank-Wolfe LP
                          overrides use_gp */
    BLGcom       *Com  /* common variables */
) ;

PRIVATE int BLGas /* return:
                   0 (affine scaling direction was computed)
                   4 (computation of AS direction fails when bracketing)
                   6 (too many iterations when computing AS direction) */
(
    BLGFLOAT       *c, /* input guess for C as parameter, return C */
    BLGFLOAT       *x, /* current iterate */
    BLGFLOAT       *g, /* gradient at x */
    BLGFLOAT lambda_k, /* constant in denominator of d */
    BLGFLOAT       *a, /* linear constraint is a'x = b */
    BLGFLOAT      *lo, /* lower bounds */
    BLGFLOAT      *hi, /* upper bounds */
    int             n, /* problem dimension */
    BLGcom       *Com  /* common variables */
) ;

PRIVATE void BLGval
(
    BLGFLOAT      *vC,
    BLGFLOAT      *dC,
    BLGFLOAT        C,
    BLGFLOAT lambda_k,
    BLGFLOAT       *x,
    BLGFLOAT       *g,
    BLGFLOAT       *a,
    BLGFLOAT      *lo,
    BLGFLOAT      *hi,
    int             n,
    BLGcom       *Com  /* common variables */
) ;

PRIVATE BLGFLOAT BLGbb
(
    BLGFLOAT alpha,
    int         mm, /* stop current cycle when Parm->MaxCycle - mm <= 0 */
    int          n, /* problem dimension */
    BLGcom    *Com  /* common variables */
) ;

PRIVATE void BLGgather
(
    BLGFLOAT *y, /* compressed vector */
    BLGFLOAT *x, /* current x */
    int      *i, /* indices associated with subspace */
    int       n  /* number of elements to extract from x */
) ;

PRIVATE void BLGscatter
(
    BLGFLOAT *x, /* scatter y into x */
    BLGFLOAT *y, /* compressed vector */
    int      *i, /* indices associated with subspace */
    int       n  /* number of elements to extract from x */
) ;

PRIVATE void BLGstep
(
    BLGFLOAT *xnew, /* updated x vector */
    BLGFLOAT    *x, /* current x */
    BLGFLOAT    *d, /* search direction */
    BLGFLOAT    st, /* stepsize */
    BLGINT       n  /* dimension */
) ;

PRIVATE void BLGstepi
(
    BLGFLOAT *xnew, /* updated x vector */
    int         *i, /* indices of xnew that change */
    BLGFLOAT    *x, /* current x */
    BLGFLOAT    *d, /* search direction */
    BLGFLOAT    st, /* stepsize */
    int         nz  /* number of components to be updated */
) ;

PRIVATE BLGFLOAT BLGdot
(
    BLGFLOAT *x, /* first vector */
    BLGFLOAT *y, /* second vector */
    BLGINT    n  /* length of vectors */
) ;

PRIVATE BLGFLOAT BLGdoti
(
    BLGFLOAT *x, /* first vector */
    int      *i, /* indices associated with elements of x */
    BLGFLOAT *y, /* second vector */
    int       n  /* length of x */
) ;

PRIVATE BLGFLOAT BLGgtd
(
    BLGFLOAT *g, /* gradient, dense */
    BLGFLOAT *a, /* linear constraint vector, dense */
    BLGFLOAT *d, /* search direction, sparse */
    BLGFLOAT  t, /* t chosen so that a'(g-ta) = 0 */
    int      *i, /* indices in g and a corresponding to indices in d */
    int      nz, /* length of d */
    int  a_is_1  /* TRUE if a is identically 1 */
) ;

PRIVATE BLGFLOAT BLGdphi
(
    BLGFLOAT *g, /* gradient, dense */
    BLGFLOAT *d, /* search direction, sparse */
    BLGcom *Com  /* common variables */
) ;

double BLGadd
(
    BLGFLOAT *x, /* vector */
    int       n  /* length of vector */
) ;

PRIVATE void BLGcopy
(
    BLGFLOAT *x, /* first vector */
    BLGFLOAT *y, /* second vector */
    BLGINT    n  /* length of vectors */
) ;

PRIVATE void BLGicopy
(
    int *x, /* copy of y */
    int *y, /* given vector */
    int  n  /* length of vectors */
) ;

BLGFLOAT BLGmax
(
    BLGFLOAT *x,
    BLGINT    n
) ;

void BLGint_init
(
    int *x,  /* array to be initialized */
    int  s,  /* scalar */
    int  n   /* length of x */
) ;

PRIVATE void BLGpartialMinSort
(
    int      *i, /* indices of m largest numbers */
    BLGFLOAT *x, /* numbers to sort */
    int       m, /* find m smallest numbers */
    int       n  /* length of x */
) ;

PRIVATE void BLGpartialMaxSort
(
    int      *i, /* indices of m largest numbers */
    BLGFLOAT *x, /* numbers to sort */
    int       m, /* find m smallest numbers */
    int       n  /* length of x */
) ;

PRIVATE void BLGminHeapSort
(
    int      *i, /* indices with decreasing order based on x */
    int   *heap, /* min heap */
    BLGFLOAT *x, /* numbers to sort */
    int       n  /* size of heap */
) ;

PRIVATE void BLGminheap_build
(
    int   *heap, /* on input, an unsorted set of indices */
    BLGFLOAT *x, /* numbers to sort */
    int       n  /* number of elements to build into the heap */
) ;

PRIVATE int BLGminheap_delete  /* return new size of heap */
(
    int   *heap, /* containing indices into x, 1..n on input */
    BLGFLOAT *x, /* not modified */
    int       n  /* number of items in heap */
) ;

PRIVATE int BLGminheap_add
(
    int    leaf, /* the new leaf */
    int   *heap, /* size n, containing indices into x */
    BLGFLOAT *x, /* not modified */
    int       n  /* number of elements in heap not counting new one */
) ;

PRIVATE void BLGminheapify
(
    int   *heap, /* size n, containing indices into x */
    BLGFLOAT *x, /* not modified */
    int       p, /* start at node p in the heap */
    int       n  /* heap [1 ... n] is in use */
) ;

PRIVATE void BLGmaxheap_build
(
    int   *heap, /* on input, an unsorted set of elements */
    BLGFLOAT *x,
    int       n  /* number of elements to build into the heap */
) ;

PRIVATE int BLGmaxheap_delete  /* return new size of heap */
(
    int   *heap, /* containing indices into x, 1..n on input */
    BLGFLOAT *x, /* not modified */
    int       n  /* number of items in heap */
) ;

PRIVATE int BLGmaxheap_add
(
    int   *heap, /* size n, containing indices into x */
    BLGFLOAT *x, /* not modified */
    int    leaf, /* the new leaf */
    int       n  /* number of elements in heap not counting new one */
) ;

PRIVATE void BLGmaxheapify
(
    int   *heap, /* size n, containing indices into x */
    BLGFLOAT *x, /* not modified */
    int       p, /* start at node p in the heap */
    int       n  /* heap [1 ... n] is in use */
) ;

int BLGnapsack /* return:
                  0 (found solution)
                  1 (feasible set empty) */
(
    BLGFLOAT      *x, /* holds y on input, and the solution x on output */
    BLGFLOAT *Lambda, /* input:  starting guess for multiplier
                         return: optimal multiplier */
    BLGFLOAT      *a, /* a'x = b (linear constraint vector) */
    BLGFLOAT       b, /* right side */
    BLGFLOAT     *lo, /* lower bounds */
    BLGFLOAT     *hi, /* upper bounds */
    int            n, /* problem dimension */
    BLGcom      *Com  /* common variables */
) ;

PRIVATE int BLGnapup /* return:
                                0 (found solution)
                                1 (feasible set empty) */
(
    BLGFLOAT        *y, /* input is target and output is solution */
    BLGFLOAT   *Lambda, /* input:  initial guess for multiplier
                           return: optimal multiplier */
    BLGFLOAT        *a, /* a'x = b (linear constraint vector) */
    BLGFLOAT         b, /* right side */
    BLGFLOAT       *lo, /* lower bounds */
    BLGFLOAT       *hi, /* upper bounds */
    int              n, /* problem dimension */
    BLGFLOAT       atx, /* a'x */
    BLGFLOAT     a2sum, /* sum_j a_j^2 */
    BLGFLOAT *breakpts, /* break points */
    int    *bound_heap, /* work array */
    int     *free_heap, /* work array */
    BLGcom        *Com  /* common variables */

) ;

PRIVATE int BLGnapdown /* return:
                                  0 (found solution)
                                  1 (feasible set empty) */
(
    BLGFLOAT        *y,
    BLGFLOAT   *Lambda, /* input:  initial guess for multiplier
                           return: optimal multiplier */
    BLGFLOAT        *a, /* a'x = b (linear constraint vector) */
    BLGFLOAT         b, /* right side */
    BLGFLOAT       *lo, /* lower bounds */
    BLGFLOAT       *hi, /* upper bounds */
    int              n, /* problem dimension */
    BLGFLOAT       atx, /* a'x */
    BLGFLOAT     a2sum, /* sum_j a_j^2 */
    BLGFLOAT *breakpts, /* break points */
    int    *bound_heap, /* work array */
    int     *free_heap, /* work array */
    BLGcom        *Com  /* common variables */
) ;

PRIVATE BLGFLOAT BLGlp1 /* return multiplier for linear constraint */
(
    BLGFLOAT  *x, /* solution */
    BLGFLOAT  mu, /* guess for multiplier */
    BLGFLOAT  *c, /* cost vector */
    BLGFLOAT   C, /* upper bound */
    BLGFLOAT   b, /* right side of linear constraint 1'x = b */
    int        n,
    int       *i, /* work array of size n */
    int   *index, /* work array of size n */
    BLGFLOAT  *y  /* work array of size n */
) ;

PRIVATE BLGFLOAT BLGlp /* return the final mu */
(
    BLGFLOAT     *x, /* holds y on input, and the solution x on output */
    BLGFLOAT     *c, /* vector c in the objective function */
    BLGFLOAT     mu, /* starting guess for mu */
    BLGFLOAT  x_bnd, /* -x_bnd <= x <= x_bnd */
    BLGFLOAT     *a, /* linear constraint is a'x = b */
    BLGFLOAT      b, /* right side */
    BLGFLOAT    *lo, /* lower bounds */
    BLGFLOAT    *hi, /* upper bounds */
    int           n, /* dimension */
    BLGcom     *Com  /* common variables */
) ;

PRIVATE void BLGindexMinSort
(
    int      *y,  /* n-by-1, indices of x giving decreasing order, output */
    int      *w,  /* n-by-1 work array, input modified */
    BLGFLOAT *x,  /* n-by-1, input not modified */
    int       n   /* number of elements to sort */
) ;

PRIVATE void BLGindexMaxSort
(
    int      *y,  /* n-by-1, indices of x giving decreasing order, output */
    int      *w,  /* n-by-1 work array, input modified */
    BLGFLOAT *x,  /* n-by-1, input not modified */
    int       n   /* number of elements to sort */
) ;

PRIVATE void BLGminSort
(
    int *x,  /* numbers to sort (input), sorted numbers (output) */
    int *y,  /* work array of length n */
    int *z,  /* work array of length n+1 */
    int  n   /* number of elements to sort */
) ;
