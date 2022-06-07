#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define BLGTRUE 1
#define BLGFALSE 0
#define BLGFLOAT double
#define BLGINT int
#define BLGINF DBL_MAX
#define BLGIINT INT_MAX
#ifndef NULL
#define NULL 0
#endif

/* ========================================================================== */
/* === structures =========================================================== */
/* ========================================================================== */

typedef struct BLGobjective_struct /* passed to user's evaluation routines */
{
    /* if objective function not quadratic: */
    BLGFLOAT     f ; /* value at x (provided by user) */
    BLGFLOAT  fold ; /* final value of f at prior iteration (from BLG) */
    BLGFLOAT    *x ; /* current iterate (from BLG) */
    BLGFLOAT *xold ; /* final x at prior iteration (from BLG) */
    BLGFLOAT    *g ; /* gradient at x (provided by user) */
    BLGFLOAT *gold ; /* previous gradient (from BLG) */
    BLGFLOAT    *d ; /* nonzero components of search direction (from BLG) */
    int         *i ; /* indices of nonzeros elements (from BLG) */
    int         nz ; /* number of nonzeros in d (from BLG) */
    BLGFLOAT  step ; /* evaluate gradient at x = xold + (*step)*d (from BLG) */
    BLGFLOAT deriv1; /* when BLG computes derivative in search direction at
                        start of iteration (info = 0), it stores it in deriv1 */
    int       info ; /* -1 = initial request for value/grad
                         0 = first evaluation in search direction
                         1 = later evaluation in search direction (from BLG) */
    int          n ; /* problem dimension (from BLG) */

    /* If objective function is quadratic (Parm->QuadCost = TRUE),
       then f(x) = .5x'Ax + c'x. The product A*d should be return in p
       during each call to the user Ad routine */
    BLGFLOAT    *p ; /* the product A*d (provided by user)*/
} BLGobjective ;

/* see README file and BLGdefault for more explanation concerning parameters */

typedef struct BLGparm_struc /* user parameters (default in BLGdefault) */
{
/* =================== PRINTING ============================================= */
    int     PrintParms ; /* TRUE means print parameter values */
    int     PrintLevel ; /* Level 0  = no print, ... , Level 2 = maximum print*/
    int     PrintFinal ; /* T => print final statistics, F => no */

/* =================== BLG LIMITS =========================================== */
    BLGINT      maxits ; /* maximum number of iterations */

/* =================== PROBLEM DESCRIPTION ================================== */
    int      QuadCost ; /* TRUE means the cost is quadratic */
    int      QuadLazy ; /* objective function is quadratic, let code compute
                           starting gradient Ax + c and starting cost
                           .5x'Ax + c'x */
    int        a_is_1 ; /* TRUE means the a vector is identically 1 */
    int bnd_is_scalar ; /* TRUE means scalar bounds */

/* =================== FEASIBILITY TESTING ================================== */
    int    check_feas ; /* TRUE means check the problem feasibility */
    BLGFLOAT feas_tol ; /* problem infeasible if |b-a'x| >= feas_tol*||a||_inf*/

/* =================== BB PARAMETER ========================================= */
    BLGFLOAT   lambda0 ; /* lower bound for lambda_k */
    BLGFLOAT lambda0Factor ;/* require that lambda_k >= lambda0Factor*||g|| */

/* =================== BB CYCLE LENGTH ====================================== */
    int            nm ; /* nominal CBB cycle length */
    int      MaxCycle ; /* absolute max number of iterations in a cycle */

/* =================== LINE SEARCH ========================================== */
    BLGFLOAT Armijo_delta ; /* f(st) <= fR+st*delta*f'(0), fR = mfax or fr */
    BLGFLOAT Armijo_decay ; /* backtracking factor in Armijo line search */
    int  max_backsteps ; /* max number of Armijo backsteps */
    BLGFLOAT     gamma ; /* required decay factor of width in bracket interval*/
    int              M ; /* memory for nonmonotone line search */
    /* T => use approximate nonmonotone Armijo line search
       F => use ordinary nonmonotone Armijo line search, switch to
            approximate Armijo when |f_r-f| < AArmijoFac*|min (f_r, f_{max})| */
    int        AArmijo ;
    BLGFLOAT AArmijoFac ;
    /* T => estimated error in function value = eps*|min (f_r, f_{max}) |
       F => estimated error in function value = eps */
    int       PertRule ;
    BLGFLOAT       eps ;
    /* search for non nan function value by shrinking search interval
       at most nshrink times */
    int       nshrink ;
    /* factor by which interval shrinks when searching for non nan value */
    BLGFLOAT  nan_fac ;
    int             L ; /* update fr if fmin was not improved after
                                L iterations */
    int             P ; /* update fr if P previous initial stepsize was
                                all accepted */
    BLGFLOAT   gamma1 ; /* criterion for updating reference value fr */
    BLGFLOAT   gamma2 ; /* criterion for updating reference value fr */
    BLGFLOAT  armijo0 ; /* criterion for quadratic interpolation in
                                cbb line search */
    BLGFLOAT  armijo1 ; /* criterion for quadratic interpolation in
                                cbb line search */

/* =================== SUBSPACE SELECTION =================================== */
    int      Subspace ; /* TRUE means use subspace implementation */
    int          nsub ; /* dimension of subspace */
    BLGFLOAT     fsub ; /* if nsub < 0, dimension = fsub*n */
    BLGFLOAT   d_pert ; /* only include d_i in subspace if >= d_pert */
    BLGFLOAT  MaxKept ; /* fraction of indices kept from prior subspace */
    BLGFLOAT   PreAmp ; /* factor >= 1 gives weight placed on prior
                           subspace indices */
    BLGFLOAT    serr1 ; /* err1 = serr1 times SVM tolerance */
    BLGFLOAT    serr2 ; /* err2 = serr2 times prior KKT error */
    BLGFLOAT bnd_pert ; /* perturbation used to decide free variables */

/* =================== SEARCH DIRECTION ===================================== */
    int            GP ; /* TRUE means use gradient projection direction.
                           Otherwise use affine scaling */
    int            FW ; /* TRUE means use Frank-Wolfe direction (has
                           precedence over GP) */
    BLGFLOAT lpfactor ; /* parameter in decision rule of QP versus LP */
    BLGFLOAT   lo_cut ; /* g_j > 0 and x_j - lo_j < lo_cut => d_j = 0 */
    BLGFLOAT   hi_cut ; /* g_j < 0 and hi_j - x_j < hi_cut => d_j = 0 */

/* =================== AFFINE SCALING CONTROL =============================== */
    BLGFLOAT    epsmu0 ; /* error tol for mu: B - A <= epsmu0*|g|/|a| */
    BLGFLOAT    epsmu1 ; /* error tol for mu: B - A <= epsmu1*(|A|+|B|) */
    int     max_mu_its ; /* max number of iterations when computing mu */
} BLGparm ;

typedef struct BLGstat_struct /* statistics returned to user */
{
    BLGFLOAT      f ; /* function value */
    BLGFLOAT lambda ; /* multiplier for constraint */
    BLGFLOAT    err ; /* max norm of g - mu d */
    BLGINT       it ; /* number of iterations */
    BLGINT       nf ; /* number of function evaluations */
    BLGINT       ng ; /* number of gradient evaluations*/
    BLGINT      nmu ; /* total number of iterations in computation of mu for
                         affine scaling algorithm */
} BLGstat ;

/* prototypes */

int BLG /* return:
                   0 (error tolerance satisfied)
                   1 (number of iterations exceeded limit)
                   2 (out of memory)
                   3 (line search failed)
                   4 (search direction not a descent direction)
                   5 (computation of search direction fails)
                   6 (search direction not a descent direction)
                   7 (computation of mu fails)
                   8 (function value is nan in line search)*/
(
    BLGFLOAT        *x, /* input: starting guess, output: the solution */
    BLGFLOAT *Lambda_k, /*either NULL or starting guess for lambda_k, if
                          not NULL, then final lambda_k is returned */
    BLGFLOAT       *lo, /* lower bounds (single bound if bndType = 1) */
    BLGFLOAT       *hi, /* upper bounds (single bound if bndType = 1) */
    BLGFLOAT        *a, /* linear constraint vector in a'x = b
                           NULL = no linear constraint or a identically 1
                           if Uparm->a_is_1 = 1 */
    BLGFLOAT         b, /* right side */
    BLGFLOAT       *g0, /* if cost quadratic, g0 = Ax + c = starting gradient */
    int              n, /* problem dimension */
    BLGFLOAT  grad_tol,/* bound on KKT error */
    /* any of the remaining arguments can be NULL */
    BLGFLOAT    *xWork, /* NULL = allocate memory internally */
    int         *iWork, /* NULL = allocate memory internally */
    BLGstat      *Stat, /* structure with statistics (can be NULL) */
    BLGparm     *Uparm, /* user parameters, NULL = use default parameters */
    void      (*Ad) (BLGobjective *), /* compute A*d if quadratic objective */
    void   (*value) (BLGobjective *), /* evaluate objective function */
    void    (*grad) (BLGobjective *), /* evaluate objective gradient */
    void (*valgrad) (BLGobjective *)  /* evaluate function and gradient
                                         NULL = use value & grad routines */
) ;

void BLGdefault
(
    BLGparm *Parm /* pointer to parameter structure */
) ;

void BLGprint_parms
(
    BLGparm     *Parm, /* SSMparm structure to be printed */
    BLGFLOAT grad_tol, /* gradient tolerance */
    int             n
) ;
