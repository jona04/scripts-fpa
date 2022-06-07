/* This example solves the problem

        minimize sum_{i=1}^n exp (x [i]) - sqrt(i+1)*x [i]
        subject to a'x = b, lo <= x <= hi
        where a_i = 1 for i <= n1 or i > n2
              a_i =-1 for n1 < i <= n2
              lo_i = 0 and hi_i = 2 for i <= n1
              lo_i = 1 and hi_i = 3 for n1 < i <= n2
              lo_i = 2 and hi_i = 4 for i > n2
   The output on a linux workstation was the following:

------------------------------------------------------
    Convergence tolerance for gradient satisfied
    
    KKT error:          7.321921e-13
    Feasibility error:  5.684342e-14
    function value:    -6.470662e+02
    
    Iterations:           42
    Function evaluations: 49
    Gradient evaluations: 42
    Mu computation:       0
    
    USE affine scaling search direction
    
    Convergence tolerance for gradient satisfied
    
    KKT error:          3.907985e-13
    Feasibility error:  4.263256e-14
    function value:    -6.470662e+02
    
    Iterations:           54
    Function evaluations: 70
    Gradient evaluations: 56
    Mu computation:       107
------------------------------------------------------

To solve the problem by optimizing over a series of subspace,
uncomment the two subspace statements below. Note that
the evaluations routines myvalue, mygrad, and myvalgrad currently
do not exploit the fact that x is only changing in a subspace */

#include <math.h>
#include "BLG_user.h"

void myvalue
(
    BLGobjective *user
) ;

void mygrad
(
    BLGobjective *user
) ;

void myvalgrad
(
    BLGobjective *user
) ;

int main (void)
{
    double *x, *lo, *hi, *a, b ;
    int i, n, n1, n2 ;
    BLGparm Parm ;

    /* allocate space for solution */
    n = 100 ;
    x  = (double *) malloc (n*sizeof (double)) ;
    lo = (double *) malloc (n*sizeof (double)) ;
    hi = (double *) malloc (n*sizeof (double)) ;
    a  = (double *) malloc (n*sizeof (double)) ;

    /* set up the constraints */
    n1 = n/3 ;
    n2 = 2*n/3 ;
    for (i = 0; i < n1; i++) { lo [i] = 0. ; hi [i] = 2. ; }
    for (; i < n2; i++) { lo [i] = 1. ; hi [i] = 3. ; }
    for (; i < n; i++) { lo [i] = 2. ; hi [i] = 4. ; }
    for (i = 0; i < n; i++) a [i] = 1. ;
    for (i = n1; i < n2; i++) a [i] = -1. ;

    /* set starting guess */
    for (i = 0; i < n; i++) x [i] = 2. ;

    /* for this test problem, we take b = a'x */
    b = 0. ;
    for (i = 0; i < n; i++) b += a [i]*x [i] ;

    BLGdefault(&Parm) ;
    Parm.PrintFinal = BLGTRUE ;
/*  Parm.Subspace = BLGTRUE ;
    Parm.nsub = 10 ;*/

    /* run the code */
    BLG (x, NULL, lo, hi, a, b, NULL, n, 1.e-12, NULL, NULL, NULL, &Parm,
         NULL, myvalue, mygrad, myvalgrad) ;

    /* with some loss of efficiency, you could omit the valgrad routine */
    Parm.PrintParms = 0 ;   /* do not print the parameters this time */
    printf ("\nUSE affine scaling search direction\n") ;
    Parm.GP = 0 ;           /* use affine scaling direction */
    for (i = 0; i < n; i++) x [i] = 2. ; /* starting guess */
    BLG (x, NULL, lo, hi, a, b, NULL, n, 1.e-12, NULL, NULL, NULL, &Parm,
         NULL, myvalue, mygrad, NULL) ;

    free (x) ;
    free (a) ;
    free (lo) ;
    free (hi) ;
    return (1) ;
}

void myvalue
(
    BLGobjective *user
)
{
    double f, t, *x ;
    int i, n ;
    x = user->x ;
    n = user->n ;
    f = 0. ;
    for (i = 0; i < n; i++)
    {
        t = i+1 ;
        t = sqrt (t) ;
        f += exp (x [i]) - t*x [i] ;
    }
    user->f = f ;
    return ;
}

void mygrad
(
    BLGobjective *user
)
{
    double t, *g, *x ;
    int i, n ;
    x = user->x ;
    g = user->g ;
    n = user->n ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        g [i] = exp (x [i]) -  t ;
    }
    return ;
}

void myvalgrad
(
    BLGobjective *user
)
{
    double ex, f, t, *g, *x ;
    int i, n ;
    f = (double) 0 ;
    x = user->x ;
    g = user->g ;
    n = user->n ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        ex = exp (x [i]) ;
        f += ex - t*x [i] ;
        g [i] = ex -  t ;
    }
    user->f = f ;
    return ;
}
