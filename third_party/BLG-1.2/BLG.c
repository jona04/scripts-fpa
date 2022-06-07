/* =========================================================================
   ================================ BLG ====================================
   =========================================================================
       ________________________________________________________________
      |Solve an optimization problem with linear and bound constraints |
      |                                                                |
      |  (BL)    min f(x)  subject to  lo <= x <= hi,  a'x = b         |
      |                                                                |
      |  The search direction is computed by either affine scaling(AS),|
      |  gradient projection (GP), or a Frank-Wolfe LP. The affine     |
      |  scaling algorithm of the code ASL is contained in BLG as a    |
      |  special case.                                                 |
      |                                                                |
      |                Version 1.0 (September 30, 2010)                |
      |                Version 1.1 (November 16, 2010)                 |
      |                Version 1.2 (February 21, 2013)                 |
      |                                                                |
      |       William W. Hager                 Hongchao Zhang          |
      |      hager@math.ufl.edu             hozhang@math.lsu.edu       |
      |   Department of Mathematics       Department of Mathematics    |
      |    University of Florida         Louisiana State University    |
      |  Gainesville, Florida 32611        Baton Rouge, LA 70803       |
      |      352-392-0281 x 244                 225-578-1982           |
      |http://www.math.ufl.edu/~hager     www.math.lsu.edu/~hozhang    |
      |                                                                |
      |      Copyright by William W. Hager and Hongchao Zhang          |
      |________________________________________________________________|

       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|

      NOTE: The code can exploit special structure:

            1. If the upper and lower bounds are scalars, not vectors,
               then set Uparm->bnd_is_scalar = TRUE. *lo and *hi each contain one
               entry, the scalar lower and upper bounds.
               Uparm->bnd_is_scalar = FALSE by default.

            2. If a = 1 identically, then set Uparm->a_is_1 = TRUE and put
               NULL for the a argument of BLG. Uparm->a_is_1 = FALSE
               by default.

            3. If f(x) = .5x'Ax + c'x (quadratic objective function),
               set Uparm->QuadCost = TRUE. When the routine Ad below is
               invoked, it should return the product A*d in Objective->p.
               Uparm->QuadCost = FALSE by default. For a quadratic
               objective, the user should provide the gradient at the
               starting x in BLG argument g0 (= Ax + c). The returned
               objective function value is the difference between the
               starting cost and the final (optimal) cost. If the user
               wants the code to compute g0 and the final cost (not the
               cost difference), then store c in the g0 argument and set
               Parm->QuadLazy = TRUE.

            4. For a general objective function, if the function value
               and gradient can be computed together at the same speed
               as the value alone, then replace the value argument of
               BLG by valgrad if valgrad is not NULL.

      The search direction is controlled by two parameters:

               Uparm->FW = TRUE means always use Frank-Wolfe search direction
               Uparm->FW = FALSE (default) means use either gradient project (GP)
                               or affine scaling (AS)
               Uparm->GP = TRUE (default) means use GP
               Uparm->GP = FALSE means use AS */
/* Contents:

    1. BLG         - solve optimization problem with linear + bound constraints
    2. BLGevaluate - evaluate objective function or gradient
    3. BLGd        - evaluate the search direction
    4. BLGval      - evaluate r(C) = a'd(C) and r'(C)
    5. BLGbb       - compute lambda_k using BB formula
    6. BLGgather   - gather and compress elements of an array
    7. BLGscatter  - scatter elements of a compressed array
    8. BLGstep     - evaluate xnew = x + st*d
    9. BLGstepi    - evaluate xnew = x + st*d when d is sparse
    10. BLGdot      - compute x dot y, x and y dense
    11. BLGdoti    - compute x dot y, x sparse and y dense
    12. BLGgtd     - compute (g - ta)'d = g'd, t = g'a/a'a given
    13. BLGdphi    - compute (g - ta)'d = g'd, t = g'a/a'a computed
    14. BLGadd     - add components of a vector
    15. BLGcopy    - copy real y to real x
    16. BLGicopy   - copy int y to int x
    17. BLGmax     - compute infinity norm of a vector
    18. BLGint_init- initialize an integer array
    19. BLGindexMinSort   - find indices that put real array in increasing order
    20. BLGindexMaxSort   - find indices that put real array in decreasing order
    21. BLGMinSort        - find indices that put int array in increasing order
    22. BLGPartialMinSort - extract some smallest elements of an array
    23. BLGpartialMaxSort - extract some largest elements of an array
    24. BLGminHeapSort    - sort a min heap
    25. BLGminheap_build  - build a min heap
    26. BLGminheap_delete - delete top element in a min heap
    27. BLGminheap_add    - add new leaf to a min heap
    28. BLGminheapify     - min heapify starting at given node
    29. BLGmaxheap_build  - build a max heap
    30. BLGmaxheap_delete - delete top element in a max heap
    31. BLGmaxheap_add    - add new leaf to a max heap
    32. BLGmaxheapify     - max heapify starting at given node
    33. BLGnapsack        - solve projection problem for napsack polytope
    34. BLGnapup          - used by napsack when multiplier guess too small
    35. BLGnapdown        - used by napsack when multiplier guess too big
    36. BLGlp1            - min c'x s.t. 1'x = b, 0 <= x <= C
    37. BLGlp             - min c'x s.t. a'x = b, max{lo, -x_bnd} <= x
                                                              <= min{hi, x_bnd}
    38. BLGdefault        - set default parameter values
    39. BLGprint_parms    - print the parameter array */

#include "BLG.h"

int BLG /* return:
                   0 (error tolerance satisfied)
                   1 (number of iterations exceeded limit)
                   2 (out of memory)
                   3 (line search failed)
                   4 (computation of AS search direction fails when bracketing)
                   5 (search direction not a descent direction )
                   6 (too many iterations when computing AS direction)
                   7 (function value is nan in line search)
                   8 (problem is infeasible)
                   9 (subproblem is infeasible)
                  10 (null subspace) */
(
    BLGFLOAT        *x, /* input: starting guess, output: the solution */
    BLGFLOAT *Lambda_k, /*either NULL or starting guess for lambda_k, if
                          not NULL, then final lambda_k is returned */
    BLGFLOAT       *lo, /* lower bounds (single bound if bnd_is_scalar = TRUE)*/
    BLGFLOAT       *hi, /* upper bounds (single bound if bnd_is_scalar = TRUE)*/
    BLGFLOAT        *a, /* linear constraint vector in a'x = b
                           NULL = no linear constraint or a identically 1
                           if Uparm->a_is_1 = TRUE */
    BLGFLOAT         b, /* right side */
    BLGFLOAT       *g0, /* if cost quadratic and Parm->QuadLazy = FALSE
                               g0 = Ax + c = starting gradient
                           if cost quadratic and Parm->QuadLazy = TRUE
                               g0 = c, linear term in cost */
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
)
{
    int fmemp, j, k, l, ia, M, max_backsteps, PrintLevel,
        status, QuadCost, mm, ll, nl, np, nz, nsub, no_a,
        AArmijo, Subspace, a_is_1, bnd_is_scalar, use_gp, use_lp, repeat,
        *i, *iwork, *iw1, *iw2, *iold ;
    long N ;
    BLGFLOAT B, fmax, fr, fdelete, fnew, f, fmin, fc, alpha, err, lambda,
          lambda0, lambda_k, atemp, normg, gtd, normx, s, t, nu, xj,
          Armijo_hi, AArmijo_hi, delta, armijo_decay, armijo0, armijo1,
          fcomp, fr_pert, fcomp_pert, dphia, lo_cut, hi_cut, Clo, Chi, SubErr,
         *xComp, *xnewComp, *gComp, *aComp, *loComp, *hiComp,
         *d, *g, *p, *xnew, *gnew, *work, *fmem, *w1, *w2 ;
    BLGINT it ;
    BLGobjective Objective ;
    BLGparm *Parm, ParmStruc ;
    BLGcom Com ;

    /* ---------------------------------------------------------------------- */
    /* set parameter values */
    /* ---------------------------------------------------------------------- */
    if ( Uparm == NULL )
    {
        Parm = &ParmStruc ;
        BLGdefault (Parm) ;
    }
    else
    {
        Parm = Uparm ;
    }

    if ( Parm->PrintParms ) BLGprint_parms (Parm, grad_tol, n) ;
    lo_cut = Parm->lo_cut ;
    hi_cut = Parm->hi_cut ;
    QuadCost = Parm->QuadCost ;
    PrintLevel = Parm->PrintLevel ;
    lambda0 = Parm->lambda0 ;
    Subspace = Parm->Subspace ;
    bnd_is_scalar = Parm->bnd_is_scalar ;
    a_is_1 = Parm->a_is_1 ;
    no_a = BLGTRUE ;
    if ( (a_is_1 == BLGTRUE) || (a != NULL) ) no_a = BLGFALSE ;
    if ( bnd_is_scalar == BLGTRUE )
    {
        Clo = *lo ;
        Chi = *hi ;
        Com.bndnorm = BLGMAX(fabs(Clo), fabs(Chi)) ;
    }
    else
    {
        s = BLGmax (lo,n) ;
        t = BLGmax (hi,n) ;
        Com.bndnorm = BLGMAX(s, t) ;
    }
    Com.a = a ;
    Com.b = b ;
    Com.lo = lo ;
    Com.hi = hi ;
    Com.tol = grad_tol ;
    Com.n = n ;

    Com.Parm = Parm ;
    N = (long) n ;
    /* compute max norm of linear constraint vector */
    if ( a_is_1 == BLGTRUE ) t = 1 ;
    else
    {
        t = BLGZERO ;
        if ( a != NULL ) t = BLGmax (a, n) ;
    }
    Com.anorm = t ;

    /* ---------------------------------------------------------------------- */
    /* allocate memory */
    /* ---------------------------------------------------------------------- */
    if ( Subspace == BLGTRUE )
    {
        /* subspace size */
        nsub = Parm->nsub ;
        if ( nsub <= 0 )
        {
            nsub = n*Parm->fsub ;
            nsub = BLGMAX (2, nsub) ;
        }
        Com.nsub = nsub ;

        /* d, g, p, xnew, gnew, w1, xnewComp, fmem, aComp, loComp, hiComp,
           xComp, gComp, w2, w3 */
        k = 14 ;
        /* no need to store a if a_is_1 */
        if ( a_is_1 == BLGTRUE ) k-- ;
        /* no need to store p if cost function is not quadratic */
        if ( QuadCost == BLGFALSE ) k-- ;
        j = 0 ;
        /* no need to store loComp and hiComp if bnd_is_scalar */
        if ( bnd_is_scalar == BLGTRUE ) { k -= 2 ; j = 2 ; }
        if ( xWork != NULL ) work = xWork ;
        else work = (BLGFLOAT *) malloc ((k*N + Parm->M + j)*sizeof(BLGFLOAT)) ;
        /* i (n), iw1, iw2, iw3, ikeep, iold (nsub)
             add 1 to iw2 due to minSort */
        if ( iWork != NULL ) iwork = iWork ;
        else iwork = (int *) malloc ((5*N+1+nsub)*sizeof (int)) ;
    }
    else
    {
        /* d, g, p, xnew, gnew, w1, xnewComp, fmem, aComp, loComp, hiComp */
        k = 10 ;
        /* no need to store a if a_is_1 */
        if ( a_is_1 == BLGTRUE ) k-- ;
        /* no need to store p if cost function is not quadratic */
        if ( QuadCost == BLGFALSE ) k-- ;
        j = 0 ;
        /* no need to store loComp and hiComp if bnd_is_scalar */
        if ( bnd_is_scalar == BLGTRUE ) { k -= 2 ; j = 2 ; }
        if ( xWork != NULL ) work = xWork ;
        else work = (BLGFLOAT *) malloc ((k*N + Parm->M + j)*sizeof(BLGFLOAT)) ;
        /* i, iw1, iw2, iw3 */
        if ( iWork != NULL ) iwork = iWork ;
        else iwork = (int *) malloc ((4*N)*sizeof (int)) ;
    }
    if ( (work == NULL) || (iwork == NULL) )
    {
        printf ("Insufficient memory for specified problem dimension %e\n",
                 (double) n) ;
        return (2) ;
    }
    /* ---------------------------------------------------------------------- */
    /* assign pointers to allocated memory */
    /* ---------------------------------------------------------------------- */
    i = iwork ;
    iw1 = i+n ;
    Com.iw1 = iw1 ;
    iw2 = iw1+n ;
    Com.iw2 = iw2 ;
    Com.iw3 = iw2+n ;
    d = work ;
    g = d+n ;
    if ( QuadCost == BLGTRUE )
    {
        p = g+n ;
        xnew = p+n ;
    }
    else
    {
        p = 0 ;
        xnew = g+n ;
    }
    gnew = xnew+n ;
    w1 = gnew+n ;
    Com.w1 = w1 ;
    xnewComp = Com.w1+n ;
    fmem = xnewComp+n ;
    if ( (a_is_1 == BLGTRUE) || (no_a == BLGTRUE) )
    {
        aComp = NULL ;
        loComp = fmem+Parm->M ;
    }
    else
    {
        aComp = fmem+Parm->M ;
        loComp = aComp+n ;
    }
    if ( bnd_is_scalar == BLGTRUE )
    {
        hiComp = loComp+1 ;
        *loComp = *lo ;
        *hiComp = *hi ;
        xComp = hiComp+1 ;
    }
    else
    {
        hiComp = loComp+n ;
        xComp = hiComp+n ;
    }
    if ( Subspace == BLGTRUE )
    {
        Com.ikeep = Com.iw3+(n+1) ; /* add 1 due to minSort */
        iold = Com.ikeep+n ;
        w2 = Com.w2 = xComp+n ;  /* for y in Subspace */
        Com.w3 = w2+n ; /* for z in Subspace */
        gComp = Com.w3+n ;
    }

    /* test for feasibility and make sure the starting guess is feasible
       by projecting on the feasible set */
    if ( Parm->check_feas == BLGTRUE )
    {
        /* check the bounds */
        if ( bnd_is_scalar == BLGTRUE )
        {
            if ( Clo > Chi ) return (8) ;
            if ( no_a == BLGTRUE )
            {
                for (j = 0; j < n; j++)
                {
                    xj = BLGMAX (x [j], Clo) ;
                    x [j] = BLGMIN (xj, Chi) ;
                }
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                if ( lo [j] > hi [j] ) return (8) ;
                if ( no_a == BLGTRUE )
                {
                    xj = BLGMAX (x [j], lo [j]) ;
                    x [j] = BLGMIN (xj, hi [j]) ;
                }
            }
        }
        if ( no_a == BLGFALSE )
        {
            t = BLGZERO ; /* starting guess for multiplier */
            if ( BLGnapsack (x, &t, a, b, lo, hi, n, &Com) ) return (8) ;
        }
    }

    /* Armijo line search parameters */
    delta = Parm->Armijo_delta ;
    armijo_decay = Parm->Armijo_decay ;
    max_backsteps = Parm->max_backsteps ;
    armijo0 = Parm->armijo0 ;
    armijo1 = Parm->armijo1 ;
    AArmijo = Parm->AArmijo ;

    /* starting function value and gradient */
    Objective.x = x ;
    Objective.g = g ;
    Objective.d = d ;
    Objective.i = i ;
    Objective.n = n ;
    Objective.info = -1 ; /* initial evaluation */
    Objective.p = p ;
    Com.user = &Objective ;
    Com.Ad = Ad ;
    Com.value = value ;
    Com.grad = grad ;
    Com.valgrad = valgrad ;
    /* initialize counters for statistics */
    it = 0 ;
    Com.nf = 0 ;
    Com.ng = 0 ;
    Com.nmu = 0 ;

    if ( QuadCost == BLGTRUE )
    {
        if ( Parm->QuadLazy == BLGTRUE )
        {
            nz = 0 ;
            for (j = 0; j < n; j++)
            {
                if ( x [j] != BLGZERO )
                {
                    d [nz] = x [j] ;
                    i [nz] = j ;
                    nz++ ;
                }
            }
            Objective.nz = nz ;
            Ad (&Objective) ;
            Com.ng = 1 ;
            /* NOTE: g0 stores linear term c in quadratic cost, g = Ax + c */
            BLGstep (g, g0, Objective.p, BLGONE, n) ;
            /* f = .5(x'g + c'x) = .5(x'(Ax+c) + c'x) = .5x'Ax + c'x */
            Objective.f = .5*(BLGdot (x, g, n) + BLGdot (x, g0, n)) ;
        }
        else
        {
            BLGcopy (g, g0, n) ;
            Objective.f = BLGZERO ;
            Com.ng = 0 ;
        }
    }
    else BLGevaluate (&Com, "fg") ; /* evaluate function and gradient */
    f = Objective.f ;

    /* henceforth, function value = fnew and gradient = gnew */
    Objective.x = xnew ;
    BLGcopy (xnew, x, n) ;
    Objective.xold = x ;
    Objective.g = gnew ;
    Objective.gold = g ;
    Objective.fold = f ;

    if ( PrintLevel >= 1 ) printf ("starting function value:  %e\n", f) ;

    /* initialize function value memory */
    M = Parm->M ;
    for (k = 1; k < M; k++) fmem [k] = -BLGINF ;
    fmem [0] = f ;
    fmemp = 1 ;    /* pointer to next available location in fmem */
    fmax = f ;
    err = BLGINF ;
    mm = 0 ; /* number of iterations in the current CBB cycle */
    ll = 0 ; /* zero as long as unit steps in line search */
    nl = 0 ; /* number of iterations since fmin decreased in value */
    np = 0 ; /* number of times initial stepsize was accepted in line search */
    fmin = f ;
    fr = f ;
    fc = f ;

    normg = BLGmax (g, n) ;
    normx = BLGmax (x, n) ;
    if ( PrintLevel >= 2 )
    {
        printf ("starting max norm of iterate: %e gradient:  %e\n",
                 normx, normg) ;
    }

    /* starting value for BB parameter */
    if ( Lambda_k == NULL )
    {
        if ( normx != BLGZERO)  lambda_k = BLGMAX (normg/normx, lambda0) ;
        else                    lambda_k = BLGMAX (normg, lambda0) ;
    }
    else lambda_k = *Lambda_k ;

    /* guess for linear constraint multiplier based on all free components */
    if      ( a_is_1 == BLGTRUE ) nu = -BLGadd (g, n)/((BLGFLOAT) n) ;
    else if ( no_a == BLGFALSE )  nu = -BLGdot (g, a, n)/BLGdot (a, a, n) ;

    /* if Parm->GP is FALSE and search direction not descent direction,
       switch use_gp to TRUE */
    while ( it < Parm->maxits )
    {
        it++ ;
        use_lp = Parm->FW ;
        use_gp = Parm->GP ;
        if ( Subspace == BLGTRUE )
        {
            if ( it == 1 ) { nsub = 0 ; SubErr = BLGZERO ; }
            else if ( nz != nsub ) BLGicopy (i, iold, nsub) ;
            repeat = 0 ;
            while ( 1 )
            {
                status = BLGbuildSubspace (i, &nsub, &err, &nu, SubErr,
                            x, g, lambda_k, &Com) ;
                if ( repeat || status || (no_a == BLGTRUE)
                            || (use_lp == BLGTRUE) ) break ;
                /* error tolerance reached, project on linear constraint */
                if ( PrintLevel >= 1 )
                {
                    printf ("Convergence tolerance satisfied (%e)\n", err) ;
                    printf ("Solve knapsack problem to improve feasibility\n") ;
                }
                repeat = 1 ;
                if ( nz != nsub )
                {
                    if ( (a_is_1 == BLGFALSE) && (a != NULL) )
                        BLGgather (aComp, a, i, nsub);
                    if ( bnd_is_scalar == BLGFALSE )
                    {
                        BLGgather (loComp, lo, i, nsub) ;
                        BLGgather (hiComp, hi, i, nsub) ;
                    }
                }
                BLGgather (xComp, x, i, nsub) ;
                BLGcopy (w2, xComp, nsub) ;

                /* evaluate right side for subspace problem, note that
                   indices in i are in increasing order */
                B = b ;
                j = 0 ;
                if ( a_is_1 == BLGTRUE )
                {
                    for (k = 0; k < nsub; k++)
                    {
                        l = i [k] ;
                        while ( j < l ) { B -= x [j] ; j++ ; }
                        j++ ;
                    }
                    while ( j < n ) { B -= x [j] ; j++ ; }
                }
                else
                {
                    for (k = 0; k < nsub; k++)
                    {
                        l = i [k] ;
                        while ( j < l ) { B -= x [j]*a [j] ; j++ ; }
                        j++ ;
                    }
                    while ( j < n ) { B -= x [j]*a [j] ; j++ ; }
                }

                /* project compressed x, stored in w2, onto feasible set,
                   napsack returns the projection in w2 */
                t = BLGZERO ;
                k = BLGnapsack (w2, &t, aComp, B, loComp, hiComp, nsub, &Com);
                if ( k ) { status = 9 ; break ; }
                /* copy the projection into xnew */
                BLGscatter (xnew, w2, i, nsub) ;
                /* compute direction from xComp (compressed x) to projection */
                BLGstep (d, w2, xComp, -BLGONE, nsub) ;
                Objective.nz = nsub ;
                Objective.fold = f ;
                Objective.step = BLGONE ;
                /* evaluate function and gradient at the new x */
                BLGevaluate (&Com, "fg") ;
                /* update the x, g, and f */
                BLGscatter (x, w2, i, nsub) ;
                BLGcopy (g, gnew, n) ;
                f = Objective.f ;
            }

            if ( (status == 0) || (status == 9) || (status == 10) ) break ;

            /* compute new search direction d, xnew = x + d, and gtd = g'd */
            if ( PrintLevel >= 1 )
            {
                printf ("subspace size before computing d: %i\n", nsub) ;
            }

            if ( status == 2 )
            {
                BLGicopy (iold, i, nsub) ; /* save new i */
                if ( (a_is_1 == BLGFALSE) && (a != NULL) )
                    BLGgather (aComp, a, i, nsub) ;
                if ( bnd_is_scalar == BLGFALSE )
                {
                    BLGgather (loComp, lo, i, nsub) ;
                    BLGgather (hiComp, hi, i, nsub) ;
                }
            }
            else /* status == 1 */
            {
                if ( nz != nsub )
                {
                   if ( (a_is_1 == BLGFALSE) && (a != NULL) )
                       BLGgather (aComp, a, i, nsub) ;
                   if ( bnd_is_scalar == BLGFALSE )
                   {
                        BLGgather (loComp, lo, i, nsub) ;
                        BLGgather (hiComp, hi, i, nsub) ;
                   }
                }
            }

            /* compress subspace elements of x */
            BLGgather (xComp, x, i, nsub) ;

            /* compress subspace elements of g */
            BLGgather (gComp, g, i, nsub) ;

            if ( it == 1 ) lambda = nu ;
            while ( 1 )
            {
                status = BLGd (d, &lambda, xnewComp, &gtd, &SubErr,
                               xComp, gComp, aComp, loComp, hiComp,
                               nsub, lambda_k, use_gp, use_lp, &Com);
                nz = Objective.nz ;
                if ( use_lp == BLGTRUE ) break ;

                if ( gtd > BLGZERO )
                {
                    if ( nz != nsub )
                    {
                       BLGicopy (i, iold, nsub) ; /* recover subspace index */
                       if ( (a_is_1 == BLGFALSE) && (a != NULL) )
                       {
                           BLGgather (aComp, a, i, nsub) ;
                       }
                       if ( bnd_is_scalar == BLGFALSE )
                       {
                            BLGgather (loComp, lo, i, nsub) ;
                            BLGgather (hiComp, hi, i, nsub) ;
                       }
                    }
                    if ( use_gp == BLGFALSE ) use_gp = BLGTRUE ;
                    else                      use_lp = BLGTRUE ;
                    if ( PrintLevel >= 1 )
                    {
                        printf ("No descent in search direction gtd: %e use_gp:"
                                " %i use_lp: %i\n", gtd, use_gp, use_lp) ;
                    }
                }
                else break ;
            }
            BLGscatter (xnew, xnewComp, i, nz) ;
        }
        else
        {
            if ( it == 1 )
            {
                lambda = nu ;
                err = grad_tol ;
            }

            /* compute new search direction d and xnew = x + d */
            repeat = 0 ;
            while ( 1 )
            {
                status = BLGd (d, &lambda, xnewComp, &gtd, &err, x, g, a,
                               lo, hi, n, lambda_k, use_gp, use_lp, &Com) ;
                if ( (no_a == BLGTRUE) || (repeat == 2) ) break ;
                /* error tolerance reached, project on linear constraint */
                if ( (err <= grad_tol) && !repeat )
                {
                    repeat = 1 ;
                    BLGcopy (xnew, x, n) ;
                    t = BLGZERO ;
                    BLGnapsack (xnew, &t, a, b, lo, hi, n, &Com) ;
                    nz = 0 ;
                    for (k = 0; k < n; k++)
                    {
                        t = xnew [k] - x [k] ;
                        if ( t != BLGZERO )
                        {
                           d [nz] = t ;
                           i [nz] = k ;
                           nz++ ;
                        }
                    }
                    Objective.nz = nz ;
                    Objective.fold = f ;
                    Objective.step = BLGONE ;

                    /* evaluate function and gradient at the new x */
                    BLGevaluate (&Com, "fg") ;
                    /* update the x, g, and f */
                    BLGcopy (x, xnew, n) ;
                    BLGcopy (g, gnew, n) ;
                    f = Objective.f ;
                    if ( PrintLevel >= 1 )
                    {
                        printf ("Convergence tolerance satisfied (%e)\n", err) ;
                        printf ("Solve knapsack prob to improve feasibility\n");
                    }
                }
                else if ( (gtd > BLGZERO) && (use_gp == BLGFALSE) )
                {
                    use_gp = BLGTRUE ;
                    repeat = 0 ;
                }
                else if ( (gtd > BLGZERO) && (use_lp == BLGFALSE) )
                {
                    use_lp = BLGTRUE ;
                    repeat = 2 ;
                }
                else break ;
                if ( (gtd > BLGZERO) && (PrintLevel >= 1) )
                {
                    printf ("No descent in search direction gtd: %e use_gp:"
                            " %i use_lp: %i\n", gtd, use_gp, use_lp) ;
                }
            }

            nz = Objective.nz ;
            if ( err < grad_tol )  break ;
            if (nz < n )   BLGscatter (xnew, xnewComp, i, nz) ;
            else           BLGcopy (xnew, xnewComp, n) ;
        }

        /* line_search: */
        if ( status ) goto Exit ;
        if ( PrintLevel >= 1 )
        {
            printf ("iter: %i f: %17.10e err: %e\n", it, f, err) ;
            printf ("nonzeros in d: %i\n", nz) ;
            if ( Subspace == BLGTRUE ) printf ("subspace error: %e\n", SubErr) ;
        }

        if ( gtd > BLGZERO ) /* search direction not a descent direction */
        {
            status = 5 ;
            break ;
        }

        /* start line search */
        alpha = BLGONE ;
        ia = 0 ;
        Objective.info = 0 ;
        Objective.fold = f ;
        Objective.deriv1 = gtd ;
        Objective.step = alpha ;
        BLGevaluate (&Com, "f") ;
        fnew = Objective.f ;
        Objective.info = 1 ;
        if ( PrintLevel >= 2 ) printf ("iter: %i iarmijo: %i gtd: %e f: "
                               "%15.6e\n", (int) it, ia, gtd, fnew) ;

        if ( nl == Parm->L )
        {
            fr = fmax ;
            t = (fr-fmin)/(fc-fmin) ;
            if ( t > Parm->gamma1 ) fr = fc ;
            nl = 0 ;
        }

        if ( (np > Parm->P) && (fmax > f) )
        {
           t = (fr-f)/(fmax-f) ;
           if ( t > Parm->gamma2 ) fr = fmax ;
        }

        /* if initial step is not a CBB step, use nonmonotone line search */
        /*  fcomp = MIN (fmax, fr) ;  */
        fcomp = f;

        /* simplified line search when objective function is quadratic */
        if ( QuadCost == BLGTRUE )
        {
            if ( AArmijo == BLGTRUE )
            {
                if ( Parm->PertRule == BLGTRUE ) t = Parm->eps*fabs(fcomp) ;
                else                             t = Parm->eps ;
                fr_pert = fr + t ;
                fr = fr_pert ;
                fcomp_pert = fcomp + t ;
                if ( Parm->PrintLevel >= 1 )
                {
                    printf ("Perform approximate Armijo line search\n") ;
                    if ( Parm->PrintLevel >= 2 )
                    {
                        printf ("fr_pert: %14.6e fcomp_pert: %14.6e\n",
                                 fr_pert, fcomp_pert) ;
                    }
                }
                AArmijo_hi = (BLGTWO*delta - BLGONE)*gtd ;
                dphia = Objective.deriv1 + Com.deriv2 ;
                if ( mm == 0 )
                {
                    if ( fnew <= fr_pert )
                    {
                        if (dphia <= BLGTWO*(fr_pert - f) + AArmijo_hi )
                        {
                            mm++ ;
                            goto exit_cbbls ;
                        }
                    }
                }
                else
                {
                    if ( fnew <= fcomp_pert )
                    {
                        if (dphia <= BLGTWO*(fcomp_pert - f) + AArmijo_hi )
                        {
                            mm++ ;
                            goto exit_cbbls ;
                        }
                    }
                }
            }
            else
            {
                if ( PrintLevel >= 1 )
                {
                    printf ("Perform ordinary Armijo line search\n") ;
                }
                Armijo_hi = delta*gtd ;
                if ( mm == 0 ) t = fr ;
                else           t = fcomp ;
                if ( fnew <= t+Armijo_hi )
                {
                    mm++ ;
                    goto exit_cbbls ;
                }
            }
            /* if unit step failed, compute exact minimizer */
            ia = 1 ;
            alpha = -gtd/Com.deriv2 ;
            fnew = f + alpha*(gtd + .5*alpha*Com.deriv2) ;
            goto exit_cbbls ;
        }

        /* Approximate nonmonotone Armijo line search, decrease alpha until:
           phi'(alpha) <= [2(phi_r - phi(0))/alpha] + (2 delta - 1) phi'(0) and
           phi(alpha) <= phi_r, where phi_r = fr_pert or fcomp_pert. */
        if ( AArmijo == BLGTRUE )
        {
            if ( Parm->PertRule == BLGTRUE ) t = Parm->eps*fabs(fcomp) ;
            else                             t = Parm->eps ;
            fr_pert = fr + t ;
            fr = fr_pert ;
            fcomp_pert = fcomp + t ;
            if ( Parm->PrintLevel >= 1 )
            {
                printf ("Perform approximate Armijo line search\n") ;
                if ( Parm->PrintLevel >= 2 )
                {
                    printf ("fr_pert: %14.6e fcomp_pert: %14.6e\n",
                             fr_pert, fcomp_pert) ;
                }
            }

            AArmijo_hi = (BLGTWO*delta - BLGONE)*gtd ;
            if ( fnew != fnew ) /* function value is nan, reduce stepsize */
            {
                for (ia = 0; ia < Parm->nshrink; ia++)
                {
                    alpha *= Parm->nan_fac ;
                    if ( nz ==n )  BLGstep (xnew, x, d, alpha, n) ;
                    else           BLGstepi (xnew, i, x, d, alpha, nz) ;
                    Objective.step = alpha ;
                    BLGevaluate (&Com, "f") ;
                    fnew = Objective.f ;
                    if ( fnew == fnew ) break ;
                }
                if ( (ia == Parm->nshrink) || (alpha == BLGZERO) )
                {
                    status = 7 ;
                    goto Exit ;
                }

                if ( fnew <= fcomp_pert)
                {
                    /* compute gradient gnew if not computed at xnew */
                    if ( value != valgrad ) BLGevaluate (&Com, "g") ;
                    dphia = BLGdphi (g, d, &Com) ;
                    if (dphia <= BLGTWO*(fcomp_pert - f)/alpha + AArmijo_hi )
                        goto exit_cbbls ; /* unit step is valid */
                }
            }
            else
            {
                if ( mm == 0 )
                {
                    if ( fnew <= fr_pert )
                    {
                        /* compute gradient gnew if not computed at xnew */
                        if ( value != valgrad ) BLGevaluate (&Com, "g") ;
                        dphia = BLGdphi (g, d, &Com) ;
                        if (dphia <= BLGTWO*(fr_pert - f) + AArmijo_hi )
                        {
                            mm++ ;
                            goto exit_cbbls ;
                        }
                    }
                }
                else
                {
                    if ( fnew <= fcomp_pert )
                    {
                        /* compute gradient gnew if not computed at xnew */
                        if ( value != valgrad ) BLGevaluate (&Com, "g") ;
                        dphia = BLGdphi (g, d, &Com) ;
                        if (dphia <= BLGTWO*(fcomp_pert - f) + AArmijo_hi )
                        {
                            mm++ ;
                            goto exit_cbbls ;
                        }
                    }
                }
            }

            while ( 1 )
            {
                ia++ ;
                /* Modified Raydan's quadratic interpolation line search */
                t = BLGTWO*(fnew-f-alpha*gtd) ;
                if ( t != BLGZERO )
                {
                    atemp = (-gtd*alpha*alpha)/t ;
                    if ( (atemp < armijo0*alpha) || (atemp > armijo1*alpha ) )
                    {
                        atemp = armijo_decay*alpha ;
                    }
                    alpha = atemp ;
                }
                else alpha *= armijo_decay ;

                if ( nz ==n )  BLGstep (xnew, x, d, alpha, n) ;
                else           BLGstepi (xnew, i, x, d, alpha, nz) ;
                Objective.step = alpha ;
                BLGevaluate (&Com, "f") ;
                fnew = Objective.f ;
                if ( PrintLevel >= 2 )
                {
                     printf ("iter: %i iarmijo: %i f: %15.6e\n",
                              (int) it, ia, fnew) ;
                }

                if ( fnew <= fcomp_pert )
                {
                    /* compute gradient gnew if not computed at xnew */
                    if ( value != valgrad ) BLGevaluate (&Com, "g") ;
                    dphia = BLGdphi (g, d, &Com) ;
                    if (dphia <= BLGTWO*(fcomp_pert - f)/alpha + AArmijo_hi )
                        goto exit_cbbls ;
                }

                if ( (alpha <= BLGZERO) || (ia >= Parm->max_backsteps) )
                {
                   status = 3 ;
                   goto Exit ;
                }
            }
            /* End of approximate Armijo line search */
        }
        /* Ordinary nonmonotone Armijo line search, decrease alpha until
           phi(alpha) <= phi_r + alpha * delta * phi'(0)
           where phi_r = fr or fcomp. */
        else
        {
            if ( PrintLevel >= 1 )
            {
                printf ("Perform ordinary Armijo line search\n") ;
            }
            Armijo_hi = delta*gtd ;
            if ( fnew != fnew ) /* function value is nan, reduce stepsize */
            {
                for (ia = 1; ia < Parm->nshrink; ia++)
                {
                    alpha *= Parm->nan_fac ;
                    if ( nz ==n )  BLGstep (xnew, x, d, alpha, n) ;
                    else           BLGstepi (xnew, i, x, d, alpha, nz) ;
                    Objective.step = alpha ;
                    BLGevaluate (&Com, "f") ;
                    fnew = Objective.f ;
                    if ( fnew == fnew ) break ;
                }
                if ( (ia == Parm->nshrink) || (alpha == BLGZERO) )
                {
                    status = 7 ;
                    goto Exit ;
                }
                if ( fnew <= fcomp+alpha*Armijo_hi )
                {
                    if ( value != valgrad ) BLGevaluate (&Com, "g") ;
                    goto exit_cbbls ;
                }
            }
            else
            {
                if ( mm == 0 ) t = fr ;
                else           t = fcomp ;
                if ( fnew <= t+Armijo_hi )
                {
                    mm++ ;
                    if ( value != valgrad ) BLGevaluate (&Com, "g") ;
                    goto exit_cbbls ;
                }
            }

            while ( 1 )
            {
                ia++ ;
                /* Modified Raydan's quadratic interpolation line search */
                t = BLGTWO*(fnew-f-alpha*gtd) ;
                if ( t != BLGZERO )
                {
                    atemp = (-gtd*alpha*alpha)/t ;
                    if ( (atemp < armijo0*alpha) || (atemp > armijo1*alpha ) )
                    {
                        atemp = armijo_decay*alpha ;
                    }
                    alpha = atemp ;
                }
                else alpha *= armijo_decay ;

                Objective.step = alpha ;
                if ( nz ==n )  BLGstep (xnew, x, d, alpha, n) ;
                else           BLGstepi (xnew, i, x, d, alpha, nz) ;
                BLGevaluate (&Com, "f") ;
                fnew = Objective.f ;
                if ( PrintLevel >= 2 )
                {
                     printf ("iter: %i iarmijo: %i f: %15.6e\n",
                              (int) it, ia, fnew) ;
                }

                if ( fnew <= fcomp+alpha*Armijo_hi )
                {
                    if ( value != valgrad ) BLGevaluate (&Com, "g") ;
                    break ;
                }

                if ( (alpha <= BLGZERO) || (ia >= Parm->max_backsteps) )
                {
                    /* try approximate Armijo line search  */
                    if ( Parm->AArmijoFac > BLGZERO )
                    {
                       fr = fcomp ;
                       printf ("Armijo linesearch fails, "
                               "try Approximate Armijo linesearch\n") ;
                    }
                    else  /* line search fails */
                    {
                        status = 3 ;
                        goto Exit ;
                    }
                }
            } /* End of Armijo line search */
        }

        exit_cbbls:
        if ( fnew < fmin )
        {
             fmin = fnew ;
             fc = fnew ;
             nl = 0 ;
        }
        else nl++ ;
        if ( fnew > fc ) fc = fnew ;

        ll += ia ;
        if ( ia == 0 ) np++ ;
        else           np = 0 ;

        /* end of cbbls */

        if ( status ) goto Exit ;

        /* if conditions hold, update BB step */
        if ( (ll >= 1) || (mm >= Parm->nm) || (it == 1)
                       || (QuadCost && Com.deriv2 < BLGZERO) )
        {
            ll = 0 ;
            t = BLGbb (alpha, mm, n, &Com) ;
            /* if t = -1, lambda_k does not change */
            if  ( t != -BLGONE )
            {
                lambda_k = t ;
                mm = 0 ;
            }
        }
        else /* set x = xnew, g = gnew */
        {
            /* update x and g */
            if ( QuadCost == BLGTRUE )
            {
                if ( alpha == BLGONE )
                {
                    for (k = 0; k < nz; k++)
                    {
                        j = i [k] ;
                        x [j] = xnew [j] ;
                    }
                }
                else
                {
                    for (k = 0; k < nz; k++) x [i [k]] += alpha*d [k] ;
                }
                BLGstep (g, g, Objective.p, alpha, n) ;
            }
            else
            {
                for (k = 0; k < nz; k++)
                {
                    j = i [k] ;
                    x [j] = xnew [j] ;
                }
                BLGcopy (g, gnew, n) ;
            }
        }

        f = fnew ;

        /* update fmax */
        fdelete = fmem [fmemp] ;
        fmem [fmemp] = fnew ;
        fmemp++ ;
        if ( fmemp == M ) fmemp = 0 ; /* memory pointer loops from M to 0 */
        if ( fmax == fdelete ) /* find the max over all the memory */
        {
             fmax = -BLGINF ;
             for (k = 0; k < M; k++) if ( fmax < fmem [k] ) fmax = fmem [k] ;
        }
        else if ( fmax < fnew ) fmax = fnew ;

        if ( AArmijo == BLGFALSE )
        {
            if ( fabs(fr - f) < Parm->AArmijoFac*fabs(fcomp) )
            {
                if (Parm->PrintLevel >= 1)
                    printf("Switch to approximate Armijo line search\n") ;
                AArmijo = BLGTRUE ;
            }
        }
    }
    if ( it == Parm->maxits ) status = 1 ; /* iteration limit exceeded */

Exit:
    if ( Lambda_k != NULL ) *Lambda_k = lambda_k ;
    if ( f > fmin) f = fmin ;
    if ( err < grad_tol )     status = 0 ; /* termination criterion satisfied */

    if ( (Parm->PrintFinal == BLGTRUE) || (PrintLevel >= 1) )
    {
        printf ("\n") ;
        if ( status == 0 )
        {
            printf ("Convergence tolerance for gradient satisfied\n") ;
        }
        else if ( status == 1 )
        {
            printf ("Number of iterations exceed specified limit\n") ;
            printf ("Iterations: %i maxit: %i\n", (int) it, Parm->maxits) ;
        }
        else if ( status == 2 )
        {
            printf ("Insufficient memory\n") ;
        }
        else if ( status == 3 )
        {
            printf ("Line search failed\n") ;
        }
        else if ( status == 4 )
        {
            printf ("Computation of search direction fails\n") ;
        }
        else if ( status == 5 )
        {
            printf ("Search direction not a descent direction (gtd = %e)\n",
                     gtd) ;
        }
        else if ( status == 6 )
        {
            printf ("Number of mu iterations exceeds limit %i\n",
                     Parm->max_mu_its) ;
        }
        else if ( status == 7 )
        {
            printf ("Function value is nan in line search\n") ;
        }
        else if ( status == 8 )
        {
            printf ("Problem is infeasible\n") ;
        }
        else if ( status == 9 )
        {
            printf ("Subproblem is infeasible\n") ;
        }
        else if ( status == 10 )
        {
            printf ("Null subspace encountered. This could happen if d_pert\n");
            printf ("is too big or the error tolerance is too small.\n") ;
        }

        printf ("\nKKT error:         %13.6e\n", err) ;
        if ( a_is_1 == BLGTRUE )
        {
            printf ("Feasibility error: %13.6e\n", fabs (BLGadd(x,n)-b)) ;
        }
        else if ( no_a == BLGFALSE )
        {
            printf ("Feasibility error: %13.6e\n", fabs (BLGdot(a,x,n)-b)) ;
        }
        printf ("function value:    %13.6e\n", f) ;
        printf ("\nIterations:           %i\n", it) ;;
        printf ("Function evaluations: %i\n", (int) Com.nf) ;
        printf ("Gradient evaluations: %i\n", (int) Com.ng) ;
        printf ("Mu computation:       %i\n", (int) Com.nmu) ;
    }

    if ( iWork == NULL ) free (iwork) ;
    if ( xWork == NULL ) free (work) ;
    if ( Stat != NULL )
    {
        Stat->f = f ;
        Stat->lambda = -lambda ;
        Stat->err = err ;
        Stat->it = it ;
        Stat->nf = Com.nf ;
        Stat->ng = Com.ng ;
        Stat->nmu = Com.nmu ;
    }
    return (status) ;
}
/* ==========================================================================
   === BLGbuildSubspace =====================================================
   ==========================================================================
    The first step in the subspace calculation is to estimate the
    multiplier nu associated with the linear constraint in the bound and
    linearly constrained problem

    (BL) min f(x)  subject to  lo <= x <= hi,  a'x = b.

    We use two different approaches depending on the size of the BB
    estimate lambda_k for the Hessian of (BL). One approach is based on
    the quadratic program

    (QP) min .5*lambda_k ||d||^2 + g'd subject to  lo <= x + d <= hi, a'd = 0,

    where g = grad f(x).  However, if lambda_k is tiny, then there
    could be rounding error problems when we try to solve this quadratic
    program. When lambda_k is sufficiently small, we drop the quadratic
    term and consider the linear program

    (LP) min g'd subject to lo <= x + d <= hi, a'd = 0.

    For the SVM problem, we use (LP) when

         max_i |g_i| >= lpfactor*C*lambda_k

    where the constraints for SVM are lo <= x <= hi.

    To solve (QP) or (LP), we make the change of variables
    z = x + d. This leads to the equivalent problems

    (QP') min ||z-y|| subject to lo <= z <= hi, a'z = b

    where y = x - g/lambda_k and b = a'x, and

    (LP') min g'z subject to lo <= z <= hi, a'z = b

    At a solution of (BL), the multiplier nu associated with the linear
    constraint is the same as the multiplier associated with linear
    constraints for (QP), (LP), and (LP'). However, for (QP') the multiplier
    mu for the linear constraint satisfies nu = lambda_k*mu.

    From the multiplier estimate, we obtain an estimate for the error
    in a approximate solution to (LB) by computing the smallest perturbation
    that allows us to satisfy the KKT conditions

         (g_i + nu) =  0 if lo_i < x_i < hi_i,
         (g_i + nu) <= 0 if x_i = hi_i,
         (g_i + nu) >= 0 if x_i = lo_i.

    We have the following ways to perturb the KKT conditions:

         (g_i + p_i + nu) =  0 if lo_i < x_i < hi_i
         (g_i + p_i + nu) <= 0 if x_i + q_i = hi_i
         (g_i + p_i + nu) >= 0 if x_i + q_i = lo_i

    If the error estimate satisfies the stopping conditions, then we stop.
    If the error in the solution to the previous subspace problem is too
    big relative to the current error, then we return to the previous
    subspace problem and further refine the error.

    The subspace is chosen in the following way. We consider the
    individual terms in the Lagrangian for either (QP) or (LP):

    t_i = .5*lambda_k d_i^2 + (g_i + nu)d_i or t_i = (g_i + nu)d_i

    Let J and K be indices corresponding to d_i > 0 and d_i < 0 respectively.
    Arrange the indices in J and K so that the t_i for i \in J or
    t_i for i \in K are in increasing order. Start with the first
    element in J or K and put it in the subspace. Then go to the
    opposite list and add indices until the d_i sum changes sign.
    Repeat this process until obtaining the desired number of indices <=
    alpha*bsub, where alpha denotes the fraction of the subspace indices
    chosen in this first phase (alpha = 2/3 default). In the second
    phase, add indices from the previous subspace until obtaining bsub
    indices if possible.  Only add indices corresponding to free variables,
    with preference to indices for which the absolute gradient is largest.
   ========================================================================== */

PRIVATE int BLGbuildSubspace /* return:
                     0 (error tolerance satisfied)
                     1 (refine accuracy for prior subspace)
                     2 (subspace has been evaluated)
                    10 (null subspace) */
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
)
{
    short qp, take1 ;
    int I1, I2, j, k, kj, n, nn, np, npre, MaxKept, Nsub, Nsub1, Nsub2, Nsort,
        PrintLevel, top, top1, top2, a_is_1, bnd_is_scalar, nfree,
       *iy, *iw1, *iw2, *ikeep, *i1, *i2, *p1, *p2 ;
    BLGFLOAT b, d_pert, dsum, err, gnorm, lpfactor, mu, Nu, PreAmp, r, s, t, u,
             bnd_pert, lo_pert, hi_pert, xj, x_bnd, X, Y, Clo, Chi, at, B,
            *a, *lo, *hi, *d, *d1, *d2, *y, *z, *z1, *z2 ;
    BLGparm *Parm ;

    /* ======== extract variables ======== */
    a_is_1 = Com->Parm->a_is_1 ;   /* TRUE means the a vector is 1 identically*/
    n = Com->n ;
    a = Com->a ;
    b = Com->b ;
    lo = Com->lo ;
    hi = Com->hi ;
    d = Com->w1 ;
    y = Com->w2 ;
    z = Com->w3 ;
    iw1 = Com->iw1 ;
    iw2 = Com->iw2 ;
    iy = Com->iw3 ;
    ikeep = Com->ikeep ; /* array initialized to 0 */
    bnd_is_scalar = Com->Parm->bnd_is_scalar ; /* TRUE means scalar bounds */
    /* ======== extract variables ======== */

    BLGint_init (ikeep, 0, n) ;
    if ( bnd_is_scalar == BLGTRUE )
    {
        Chi = *hi ;
        Clo = *lo ;
    }

    Parm = Com->Parm ;
    PrintLevel = Parm->PrintLevel ;
    d_pert = Com->tol*Parm->d_pert ;
    lpfactor = Parm->lpfactor ;
    Nsub = Com->nsub ;          /* nominal subspace dimension */
    MaxKept = Nsub*Parm->MaxKept ;
    PreAmp = Parm->PreAmp ;
    bnd_pert = Parm->bnd_pert ;
    /* ======== end of extract ======== */

    gnorm = sqrt (BLGdot (g, g, n)) ;

    if ( (a_is_1 == BLGTRUE) || (a != NULL) ) /* linear constraint present */
    {
        if ( a_is_1 == BLGTRUE ) B = BLGadd(x, n) ;
        else                     B = BLGdot(a, x, n) ;
        qp = 0 ;
        /* if the box that you are projecting onto is tiny compared to
           the vector that you are projecting, then it is better to solve
           an LP gotten by neglecting the quadratic term in (QP) */
        if ( gnorm >= lpfactor*lambda_k*Com->bndnorm ) /* solve LP */
        {
            if ( PrintLevel >= 1 )
            {
                printf ("Use LP in buildSubspace\n") ;
            }
            /* obtain bound on infinity norm of solution based on the
               following bound for the solution to the napsack problem:
               ||x||_inf <= ||x0||_inf + 2||g||_2/lambda_k
               here x0 is any feasible point */
            x_bnd = BLGmax (x, n) + 2.*gnorm/lambda_k ;
            Nu = BLGlp (y, g, *nu, x_bnd, a, B, lo, hi, n, Com) ;
        }
        else                        /* solve QP */
        {
            if ( PrintLevel >= 1 )
            {
                printf ("Use QP in buildSubspace\n") ;
            }
            BLGstep (y, x, g, -BLGONE/lambda_k, n) ;
            {
                mu = *nu/lambda_k ;
                BLGnapsack (y, &mu, a, B, lo, hi, n, Com) ;
                Nu = lambda_k*mu ;
                qp = 1 ;
            }
        }
        *nu = Nu ;
    }
    else                       /* no linear constraint */
    {
        qp = -1 ;
        BLGstep (y, x, g, -BLGONE/lambda_k, n) ;
        if ( bnd_is_scalar == BLGTRUE )
        {
             for (j = 0; j < n; j++)
             {
                 t = y [j] ;
                 if      ( t < Clo ) y [j] = Clo ;
                 else if ( t > Chi ) y [j] = Chi ;
             }
        }
        else
        {
             for (j = 0; j < n; j++)
             {
                 t = y [j] ;
                 if      ( t < lo [j] ) y [j] = lo [j] ;
                 else if ( t > hi [j] ) y [j] = hi [j] ;
             }
        }
    }

    nfree = 0 ;
    /* flag prior free indices */
    if  ( bnd_is_scalar == BLGTRUE )
    {
        lo_pert = Clo + bnd_pert ;
        hi_pert = Chi - bnd_pert ;
        for (j = 0; j < *nsub; j++)
        {
            k = i [j] ;
            if ( (x [k] > lo_pert) && (x [k] < hi_pert) )
            {
                ikeep [k] = -1 ;
                nfree++ ;
            }
       }
    }
    else
    {
        for (j = 0; j < *nsub; j++)
        {
            k = i [j] ;
            if ( (x [k] > lo [k] + bnd_pert) && (x [k] < hi [k] - bnd_pert) )
            {
                ikeep [k] = -1 ;
                nfree++ ;
            }
        }
    }
    /* y now stores z, compute the following:
       1. error estimate based on the smallest change
          in g_j or xj for which optimality conditions hold
       2. d = y - x
       3. the subspace is based on the smallest elements of
          .5 lambda_k d_j^2 + (g_j + nu a_j)d_j or
                              (g_j + nu a_j)d_j */

    err = BLGZERO ;
    s = .5*lambda_k ;
    nn = 0 ;
    np = n ;
    if ( a_is_1 == BLGTRUE )
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
           for (j = 0; j < n; j++)
           {
              xj = x [j] ;
              t = y [j] - xj ;
              r = g [j] + Nu ;
              if ( fabs (t) >= d_pert ) /* include in d */
              {
                  if ( qp ) u = t*(s*t + r) ;
                  else      u = t*r ;
                  if ( t > BLGZERO )
                  {
                      np-- ;
                      if ( ikeep [j] ) z [np] = PreAmp*u ;
                      else             z [np] = u ;
                      iy [np] = j ;
                      d [np] = t ;
                  }
                  else
                  {
                      if ( ikeep [j] ) z [nn] = PreAmp*u ;
                      else             z [nn] = u ;
                      iy [nn] = j ;
                      d [nn] = t ;
                      nn++ ;
                  }
              }
              if ( r > BLGZERO )
              {
                  X = xj - Clo ;
                  if ( r > X ) err = BLGMAX (err, X) ;
                  else         err = BLGMAX (err, r) ;
              }
              else
              {
                  Y = Chi - xj ;
                  if (-r > Y ) err = BLGMAX (err, Y) ;
                  else         err = BLGMAX (err,-r) ;
              }
           }
        }
        else
        {
           for (j = 0; j < n; j++)
           {
              xj = x [j] ;
              t = y [j] - xj ;
              r = g [j] + Nu ;
              if ( fabs (t) >= d_pert ) /* include in d */
              {
                  if ( qp ) u = t*(s*t + r) ;
                  else      u = t*r ;
                  if ( t > BLGZERO )
                  {
                      np-- ;
                      if ( ikeep [j] ) z [np] = PreAmp*u ;
                      else             z [np] = u ;
                      iy [np] = j ;
                      d [np] = t ;
                  }
                  else
                  {
                      if ( ikeep [j] ) z [nn] = PreAmp*u ;
                      else             z [nn] = u ;
                      iy [nn] = j ;
                      d [nn] = t ;
                      nn++ ;
                  }
              }
              if ( r > BLGZERO )
              {
                  X = xj - lo [j] ;
                  if ( r > X ) err = BLGMAX (err, X) ;
                  else         err = BLGMAX (err, r) ;
              }
              else
              {
                  Y = hi [j] - xj ;
                  if (-r > Y ) err = BLGMAX (err, Y) ;
                  else         err = BLGMAX (err,-r) ;
              }
           }
        }
    }
    else if ( a == NULL ) /* no linear constraint */
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
           for (j = 0; j < n; j++)
           {
              xj = x [j] ;
              t = y [j] - xj ;
              r = g [j] ;
              if ( fabs (t) >= d_pert ) /* include in d */
              {
                  u = t*(s*t + r) ;
                  if ( ikeep [j] ) z [nn] = PreAmp*u ;
                  else             z [nn] = u ;
                  iy [nn] = j ;
                  nn++ ;
              }
              if ( r > BLGZERO )
              {
                  X = xj - Clo ;
                  if ( r > X ) err = BLGMAX (err, X) ;
                  else         err = BLGMAX (err, r) ;
              }
              else
              {
                  Y = Chi - xj ;
                  if (-r > Y ) err = BLGMAX (err, Y) ;
                  else         err = BLGMAX (err,-r) ;
              }
           }
        }
        else
        {
           for (j = 0; j < n; j++)
           {
              xj = x [j] ;
              t = y [j] - xj ;
              r = g [j] ;
              if ( fabs (t) >= d_pert ) /* include in d */
              {
                  u = t*(s*t + r) ;
                  if ( ikeep [j] ) z [nn] = PreAmp*u ;
                  else             z [nn] = u ;
                  iy [nn] = j ;
                  nn++ ;
              }
              if ( r > BLGZERO )
              {
                  X = xj - lo [j] ;
                  if ( r > X ) err = BLGMAX (err, X) ;
                  else         err = BLGMAX (err, r) ;
              }
              else
              {
                  Y = hi [j] - xj ;
                  if (-r > Y ) err = BLGMAX (err, Y) ;
                  else         err = BLGMAX (err,-r) ;
              }
           }
        }
    }
    else                  /* linear constraint and a != 1 */
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
           for (j = 0; j < n; j++)
           {
              xj = x [j] ;
              t = y [j] - xj ;
              r = g [j] + Nu*a[j] ;
              if ( fabs (t) >= d_pert ) /* include in d */
              {
                  if ( qp ) u = t*(s*t + r) ;
                  else      u = t*r ;
                  at = a[j]*t ;
                  if ( at > BLGZERO )
                  {
                      np-- ;
                      if ( ikeep [j] ) z [np] = PreAmp*u ;
                      else             z [np] = u ;
                      iy [np] = j ;
                      d [np] = at ;
                  }
                  else
                  {
                      if ( ikeep [j] ) z [nn] = PreAmp*u ;
                      else             z [nn] = u ;
                      iy [nn] = j ;
                      d [nn] = at ;
                      nn++ ;
                  }
              }
              if ( r > BLGZERO )
              {
                  X = xj - Clo ;
                  if ( r > X ) err = BLGMAX (err, X) ;
                  else         err = BLGMAX (err, r) ;
              }
              else
              {
                  Y = Chi - xj ;
                  if (-r > Y ) err = BLGMAX (err, Y) ;
                  else         err = BLGMAX (err,-r) ;
              }
           }
        }
        else
        {
           for (j = 0; j < n; j++)
           {
              xj = x [j] ;
              t = y [j] - xj ;
              r = g [j] + Nu*a[j] ;
              if ( fabs (t) >= d_pert ) /* include in d */
              {
                  if ( qp ) u = t*(s*t + r) ;
                  else      u = t*r ;
                  at = a[j]*t ;
                  if ( at > BLGZERO )
                  {
                      np-- ;
                      if ( ikeep [j] ) z [np] = PreAmp*u ;
                      else             z [np] = u ;
                      iy [np] = j ;
                      d [np] = at ;
                  }
                  else
                  {
                      if ( ikeep [j] ) z [nn] = PreAmp*u ;
                      else             z [nn] = u ;
                      iy [nn] = j ;
                      d [nn] = at ;
                      nn++ ;
                  }
              }
              if ( r > BLGZERO )
              {
                  X = xj - lo [j] ;
                  if ( r > X ) err = BLGMAX (err, X) ;
                  else         err = BLGMAX (err, r) ;
              }
              else
              {
                  Y = hi [j] - xj ;
                  if (-r > Y ) err = BLGMAX (err, Y) ;
                  else         err = BLGMAX (err,-r) ;
              }
           }
        }
    }

    *Err = err ;
     if ( PrintLevel > 0 ) printf ("    error:                  %e\n", err) ;

    /* check if termination condition is satisfied */
    if ( err <= Com->tol ) return (0) ;

    /* check whether the prior subspace should be repeated */
    if ( SubErr > BLGMAX (Parm->serr1*Com->tol, Parm->serr2*err ) )
    {
        if ( PrintLevel >= 1 ) printf ("    repeat subspace\n") ;
        return (1) ;
    }

    if (nn == 0 || ( (np == n) && ((a_is_1 == BLGTRUE) || (a != NULL)) ) )
        return (10) ;

    if ( Nsub >= nn + (n - np) ) /* put all the indices in the subspace */
    {
        npre = nfree ;
        for (k = 0; k < nn; k++) i [k] = iy [k] ;
        Nsub = nn ;
        for (k = np; k < n; k++)
        {
            j = iy [k] ;
            i [Nsub] = j ;
            Nsub++ ;
            ikeep [j] = 1 ;
        }
    }
    /* need to choose indices for subspace */
    else if ( (a_is_1 == BLGTRUE) || (a != NULL) )
    {
        /* find Nsub smallest elements in the bottom part of z */
        Nsub1 = BLGMIN(Nsub, nn) ;
        BLGpartialMinSort (iw2, z, Nsub1, nn) ;
        BLGminHeapSort (iw1, iw2, z, Nsub1) ;

        /* find Nsub smallest elements in the top part of z */
        Nsub2 = BLGMIN(Nsub, n-np) ;
        BLGpartialMinSort (iw2+np, z+np, Nsub2, n-np) ;
        BLGminHeapSort (iw1+np, iw2+np, z+np, Nsub2) ;

        z1 = z ;     /* cost improvement when d < 0 */
        i1 = iw1 ;   /* indices giving increasing order in z1 */
        p1 = iy ;    /* original indices */
        d1 = d ;     /* direction y - x */

        z2 = z+np ;  /* cost improvement when d > 0 */
        i2 = iw1+np ;/* indices giving increasing order in z2 */
        p2 = iy+np ; /* original indices */
        d2 = d+np ;  /* direction y - x */

        npre = 0 ;   /* # of previous free subspace indices in new subspace */
        k = 0 ;
        top1 = 0 ;
        top2 = 0 ;
        I1 = i1 [0] ;
        I2 = i2 [0] ;
        dsum = BLGZERO ;

        take1 = 1 ;
        if ( z1 [I1] > z2 [I2] ) take1 = 0 ;
        while ( k < Nsub )
        {
            if ( take1 )
            {
                I1 = i1 [top1] ;
                top1++ ;
                j = p1 [I1] ;
                kj = ikeep [j] ;
                if ( !kj || (kj && (npre < MaxKept)) )
                {
                    if ( kj ) npre++ ;
                    dsum += d1 [I1] ;
                    i [k] = j ;
                    ikeep [j] = 1 ;
                    k++ ;
                }
            }
            else
            {
                I2 = i2 [top2] ;
                top2++ ;
                j = p2 [I2] ;
                kj = ikeep [j] ;
                if ( !kj || (kj && (npre < MaxKept)) )
                {
                    if ( kj ) npre++ ;
                    dsum += d2 [I2] ;
                    i [k] = j ;
                    ikeep [j] = 1 ;
                    k++ ;
                }
            }
            if ( (top1 == Nsub1) || (top2 == Nsub2) ) break ;
            if      ( dsum > BLGZERO ) take1 = 1 ;
            else if ( dsum < BLGZERO ) take1 = 0 ;
            else
            {
                if ( z1 [I1] < z2 [I2] ) take1 = 1 ;
                else                     take1 = 0 ;
            }
        }
        if ( k < Nsub )
        {
            if ( top1 == Nsub1 )
            {
                while ( (k < Nsub) && (top2 < Nsub2) )
                {
                    j = p2 [i2 [top2]] ;
                    top2++ ;
                    kj = ikeep [j] ;
                    if ( !kj || (kj && (npre < MaxKept)) )
                    {
                        if ( kj ) npre++ ;
                        i [k] = j ;
                        ikeep [j] = 1 ;
                        k++ ;
                    }
                }
            }
            else
            {
                while ( (k < Nsub) && (top1 < Nsub1) )
                {
                    j = p1 [i1 [top1]] ;
                    top1++ ;
                    kj = ikeep [j] ;
                    if ( !kj || (kj && (npre < MaxKept)) )
                    {
                        if ( kj ) npre++ ;
                        i [k] = j ;
                        ikeep [j] = 1 ;
                        k++ ;
                    }
                }
            }
        }
        Nsub = k ;
    }
    else    /* a not present, need to choose indices for subspace */
    {
        /* find Nsub+MaxKept smallest elements in the bottom part of z */
        Nsort = BLGMIN(Nsub+MaxKept, nn) ;
        BLGpartialMinSort (iw2, z, Nsort, nn) ;
        BLGminHeapSort (iw1, iw2, z, Nsort) ;
        npre = 0 ;   /* # of previous free subspace indices in new subspace */
        k = 0 ;
        top = 0 ;
        while ( k < Nsub )
        {
            if ( top == Nsort ) break ;
            j = iy [iw1[top]] ;
            top++ ;
            kj = ikeep [j] ;
            if ( !kj || (kj && (npre < MaxKept)) )
            {
                if ( kj ) npre++ ;
                i [k] = j ;
                k++ ;
            }
        }
        Nsub = k ;
    }

    if (PrintLevel > 1)
    {
        printf ("    prior subspace indices: %i\n", npre) ;
        printf ("    total subspace size:    %i\n", Nsub) ;
    }

    /* At this point, the subspace is complete.
       Nsub = nsub = dimension of current subspace */
    *nsub = Nsub ;

    /* sort subspace indices in increasing order */
    BLGminSort (i, iw1, Com->iw3, Nsub) ;

    return (2) ;
}

/* ==========================================================================
   === BLGd =================================================================
   ==========================================================================
   This routine constructs the search direction. For a bound constrained
   problem without linear constraints, the search direction is simply
   the gradient projection search direction. When the linear constraint
   is present, there are 3 different ways to construct the search direction:

   1. The Frank-Wolfe (FW) search direction which corresponds to solving the LP

      (LP) min g'd subject to lo <= x + d <= hi, a'd = 0.

      To ensure that the LP has a solution, we impose a bound chosen
      so that ||x + d||_inf <= an upper bound for the solution to the QP.

   2. The gradient projection (GP) search direction which corresponds
      to solving a QP

      (QP) min .5*lambda_k ||d||^2 + g'd subject to  lo <= x + d <= hi, a'd = 0,

      Here lambda_k is the BB parameter.

   3. The affine scaling (AS) search direction d(mu) with the property that
      a'd(mu) = 0 where

                             g_j - mu a_j
         d_j (mu) = - ----------------------------
                       lambda_k + |g_j - mu a_j|/X

         X = hi_j - x_j if g_j - mu a_j <= 0
         X = x_j - lo_j if g_j - mu a_j  > 0

    The choice of the search direction is based on the logical parameters
    use_lp and use_gp. If use_lp is TRUE, then FW is used. If use_lp
    is FALSE and use_gp is TRUE, the GP is nominally used. If use_lp and
    use_gp are both FALSE, then AS is nominally used. When lambda_k is
    tiny, it may be difficult to compute the GP or AS directions accurately.
    In this case, we use FW since the solution to the LP may be a
    more accurate approximation to the GP or the AS directions
    than what either the GP or AS routines would give.

    For more information concerning the solution of either (LP) or (GP),
    see the comments that precede BLGbuildSubspace. For more information
    concerning the AS direction, see the comments that precede BLGas.

    After generating one of the 3 directions, we then convert the search
    direction to a sparse vector and compute gtd = g'd.
   ========================================================================== */
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
)
{
    int j, nz, a_is_1, bnd_is_scalar, PrintLevel, QuadCost,
        Subspace, status, *i ;
    BLGFLOAT a2sum, aj, B, C, gtd, err, r, t, dj, gj, pj, xj, loj, hij,
             Clo, Chi, rk, X, Y, gnorm, lo_cut, hi_cut, psi, x_bnd, xnewj ;
    BLGparm *Parm ;
    i = Com->user->i ;
    Parm = Com->Parm ;
    a_is_1 = Parm->a_is_1 ;
    bnd_is_scalar = Parm->bnd_is_scalar ;
    QuadCost = Parm->QuadCost ;
    Subspace = Parm->Subspace ;
    PrintLevel = Parm->PrintLevel ;

    lo_cut = Parm->lo_cut ;
    hi_cut = Parm->hi_cut ;
    if ( bnd_is_scalar == BLGTRUE )
    {
        Clo = *lo ;
        Chi = *hi ;
    }
    gtd = BLGZERO ;
    err = BLGZERO ;
    nz = 0 ;

    /* if no linear constraint, the search direction is obtained
       by gradient projection onto the feasible set */
    if ( (a == NULL) && (a_is_1 == BLGFALSE) )
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
            for (j = 0; j < n; j++)
            {
                gj = g [j] ;
                xj = x [j] ;
                X = xj - Clo ;
                Y = Chi - xj ;
                if ( gj > BLGZERO )
                {
                    if ( gj > X ) err = BLGMAX (err, X) ;
                    else          err = BLGMAX (err, gj) ;
                    if  ( X > lo_cut )
                    {
                        if ( use_gp == BLGTRUE ) t = -BLGMIN(gj/lambda_k, X) ;
                        else       { r = gj/X ;  t = -gj/(lambda_k + r) ; }
                        if ( Y-t <= hi_cut ) dj = BLGZERO;
                        else                 dj = t ;
                    }
                    else                     dj = BLGZERO ;
                }
                else
                {
                    t = -gj ;
                    if ( t > Y ) err = BLGMAX (err, Y) ;
                    else         err = BLGMAX (err, t) ;
                    if ( Y > hi_cut )
                    {
                        if ( use_gp == BLGTRUE ) t = BLGMIN(t/lambda_k, Y) ;
                        else        { r = t/Y ;  t = t/(lambda_k + r) ; }
                        if ( X+t <= lo_cut ) dj = BLGZERO;
                        else                 dj = t ;
                    }
                    else                     dj = BLGZERO ;
                }
                /* only store nonzeros in d */
                if ( dj != BLGZERO )
                {
                    d [nz] = dj ;
                    if ( QuadCost == BLGFALSE ) xnew [nz] = xj + dj ;
                    if ( Subspace == BLGTRUE ) i [nz] = i [j] ;
                    else                       i [nz] = j ;
                    gtd += dj*gj ;
                    nz++ ;
                }
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                gj = g [j] ;
                xj = x [j] ;
                X = xj - lo [j] ;
                Y = hi [j] - xj ;
                if ( gj > BLGZERO )
                {
                    if ( gj > X ) err = BLGMAX (err, X) ;
                    else          err = BLGMAX (err, gj) ;
                    if  ( X > lo_cut )
                    {
                        if ( use_gp == BLGTRUE )   t = -BLGMIN(gj/lambda_k, X) ;
                        else         { r = gj/X ;  t = -gj/(lambda_k + r) ; }
                        if ( Y-t <= hi_cut ) dj = BLGZERO;
                        else                 dj = t ;
                    }
                    else                     dj = BLGZERO ;
                }
                else
                {
                    t = -gj ;
                    if ( t > Y ) err = BLGMAX (err, Y) ;
                    else         err = BLGMAX (err, t) ;
                    if ( Y > hi_cut )
                    {
                        if ( use_gp == BLGTRUE ) t = BLGMIN(t/lambda_k, Y) ;
                        else        { r = t/Y ;  t = t/(lambda_k + r) ; }
                        if ( X+t <= lo_cut ) dj = BLGZERO;
                        else                 dj = t ;
                    }
                    else                     dj = BLGZERO ;
                }
                /* only store nonzeros in d */
                if ( dj != BLGZERO )
                {
                    d [nz] = dj ;
                    if ( QuadCost == BLGFALSE ) xnew [nz] = xj + dj ;
                    if ( Subspace == BLGTRUE )  i [nz] = i [j] ;
                    else                        i [nz] = j ;
                    gtd += dj*gj ;
                    nz++ ;
                }
            }
        }
        *Err = err ;
        *Gtd = gtd ;
        Com->user->nz = nz ;
        return (0) ;
    }

    gnorm = sqrt (BLGdot(g, g, n)) ;
    if ( gnorm == BLGZERO )
    {
        *Err = BLGZERO ;
        return (0) ;
    }
    if ( gnorm >= Parm->lpfactor*Com->bndnorm*lambda_k)  use_lp = BLGTRUE ;

    /* Compute right side of linear constraint to ensure descent.
       The error in the linear constraint is corrected at the end */
    if ( a_is_1 == BLGTRUE ) B = BLGadd(x, n) ;
    else                     B = BLGdot(a, x, n) ;

    if ( use_lp == BLGTRUE )  /* use lp to find d */
    {
        if ( PrintLevel >= 1 )
        {
            printf ("Use LP to compute search direction\n") ;
        }
        /* obtain bound on infinity norm of solution based on the
           following bound for the solution to the napsack problem:
           ||x||_inf <= ||x0||_inf + 2||g||_2/lambda_k
           here x0 is any feasible point */
        x_bnd = BLGmax (x, n) + 2.*gnorm/lambda_k ;
        C = BLGlp (d, g, *mu, x_bnd, a, B, lo, hi, n, Com) ;
    }
    else if ( use_gp == BLGTRUE ) /* compute the gradient projection direction*/
    {
        if ( PrintLevel >= 1 )
        {
            printf ("Use gradient projection search direction\n") ;
        }
        rk = -BLGONE/lambda_k ;
        BLGstep (d, x, g, rk, n) ;
        psi = *mu/lambda_k ;
        BLGnapsack (d, &psi, a, B, lo, hi, n, Com) ;
        C = lambda_k*psi ;
    }
    else           /* compute affine scaling direction */
    {
        if ( PrintLevel >= 1 )
        {
            printf ("Use affine scaling search direction\n") ;
        }
        C = *mu ;
        status = BLGas (&C, x, g, lambda_k, a, lo, hi, n, Com) ;
        if ( status ) return (status) ;
    }

    /* compute the search direction and the error at x */
    gtd = BLGZERO ;
    a2sum = BLGZERO ;
    if ( bnd_is_scalar == BLGTRUE )
    {
        for (j = 0; j < n; j++)
        {
            gj = g [j] ;
            if ( a_is_1 == BLGTRUE )       t = gj + C ;
            else            { aj = a [j] ; t = gj + C*aj ; }
            pj = t ;
            xj = x [j] ;
            dj = d [j] ;
            xnewj = dj ;
            X = xj - Clo ;
            Y = Chi - xj ;
            if ( use_lp == BLGTRUE )
            {
               if ( t > BLGZERO ) err = BLGMAX (err, BLGMIN(t, X )) ; /* KKT */
               else               err = BLGMAX (err, BLGMIN(-t,Y )) ;
               dj -= xj ;
            }
            else
            {
               if ( t > BLGZERO )
               {
                   if ( t > X ) err = BLGMAX (err, X) ; /* KKT error */
                   else         err = BLGMAX (err, t) ;
                   if  ( X > lo_cut )
                   {
                       if ( use_gp == BLGTRUE ) /* GP */
                       {
                           if      ( dj == Clo ) t = dj - xj ;
                           else if ( dj == Chi ) t = dj - xj ;
                           else
                           {
                               if ( a_is_1 == BLGTRUE ) t = rk*gj - psi ;
                               else                     t = rk*gj - psi*aj ;
                           }
                       }
                       else /* AF */
                       {
                           r = t/X ;
                           t = -t/(lambda_k + r) ;
                           xnewj = xj + t ;
                       }
                       if ( Y-t <= hi_cut )
                       {
                           dj = BLGZERO ;
                           xnewj = Chi ;
                       }
                       else dj = t ;
                   }
                   else
                   {
                       xnewj = Clo ;
                       dj = BLGZERO ;
                   }
               }
               else
               {
                   t = -t ;
                   if ( t > Y ) err = BLGMAX (err, Y) ; /* KKT error */
                   else         err = BLGMAX (err, t) ;
                   if ( Y > hi_cut )
                   {
                       if ( use_gp == BLGTRUE ) /* GP */
                       {
                           if      ( dj == Clo ) t = dj - xj ;
                           else if ( dj == Chi ) t = dj - xj ;
                           else
                           {
                               if ( a_is_1 == BLGTRUE ) t = rk*gj - psi ;
                               else                     t = rk*gj - psi*aj ;
                           }
                       }
                       else /* AF */
                       {
                           r = t/Y ;
                           t = t/(lambda_k + r) ;
                           xnewj = xj + t ;
                       }
                       if ( X+t <= lo_cut )
                       {
                           dj = BLGZERO ;
                           xnewj = Clo ;
                       }
                       else dj = t ;
                   }
                   else
                   {
                       xnewj = Chi ;
                       dj = BLGZERO ;
                   }
               }
            }
            /* only store nonzeros in d */
            if ( dj != BLGZERO )
            {
                gtd += pj*dj ;
                if ( a_is_1 == BLGFALSE ) a2sum += aj*aj ;
                d [nz] = dj ;
                xnew [nz] = xnewj ;
                if ( Subspace == BLGTRUE ) i [nz] = i [j] ;
                else                       i [nz] = j ;
                nz++ ;
            }
        }
    }
    else
    {
        for (j = 0; j < n; j++)
        {
            gj = g [j] ;
            if ( a_is_1 == BLGTRUE ) t = gj + C ;
            else      { aj = a [j] ; t = gj + C*aj ; }
            pj = t ;
            xj = x [j] ;
            dj = d [j] ;
            xnewj = dj ;
            loj = lo [j] ;
            X = xj - loj ;
            hij = hi [j] ;
            Y = hij - xj ;
            if ( use_lp == BLGTRUE )
            {
               if ( t > BLGZERO ) err = BLGMAX (err, BLGMIN(t, X )) ; /* KKT */
               else               err = BLGMAX (err, BLGMIN(-t,Y )) ;
               dj -= xj ;
            }
            else
            {
               if ( t > BLGZERO )
               {
                   if ( t > X ) err = BLGMAX (err, X) ; /* KKT error */
                   else         err = BLGMAX (err, t) ;
                   if  ( X > lo_cut )
                   {
                       if ( use_gp == BLGTRUE ) /* GP */
                       {
                           if      ( dj == loj ) t = dj - xj ;
                           else if ( dj == hij ) t = dj - xj ;
                           else
                           {
                               if ( a_is_1 == BLGTRUE ) t = rk*gj - psi ;
                               else                     t = rk*gj - psi*a [j] ;
                           }
                       }
                       else /* AF */
                       {
                           r = t/X ;
                           t = -t/(lambda_k + r) ;
                           xnewj = xj + t ;
                       }
                       if ( Y-t <= hi_cut )
                       {
                           dj = BLGZERO ;
                           xnewj = hij ;
                       }
                       else dj = t ;
                   }
                   else
                   {
                       xnewj = loj ;
                       dj = BLGZERO ;
                   }
               }
               else
               {
                   t = -t ;
                   if ( t > Y ) err = BLGMAX (err, Y) ; /* KKT error */
                   else         err = BLGMAX (err, t) ;
                   if ( Y > hi_cut )
                   {
                       if ( use_gp == BLGTRUE ) /* GP */
                       {
                           if      ( dj == loj ) t = dj - xj ;
                           else if ( dj == hij ) t = dj - xj ;
                           else
                           {
                               if ( a_is_1 == BLGTRUE ) t = rk*gj - psi ;
                               else                     t = rk*gj - psi*a [j] ;
                           }
                       }
                       else /* AF */
                       {
                           r = t/Y ;
                           t = t/(lambda_k + r) ;
                           xnewj = xj + t ;
                       }
                       if ( X+t <= lo_cut )
                       {
                           dj = BLGZERO ;
                           xnewj = loj ;
                       }
                       else dj = t ;
                   }
                   else
                   {
                       xnewj = hij ;
                       dj = BLGZERO ;
                   }
               }
            }
            /* only store nonzeros in d */
            if ( dj != BLGZERO )
            {
                gtd += pj*dj ;
                if ( a_is_1 == BLGFALSE ) a2sum += aj*aj ;
                d [nz] = dj ;
                xnew [nz] = xnewj ;
                if ( Subspace == BLGTRUE ) i [nz] = i [j] ;
                else                       i [nz] = j ;
                nz++ ;
            }
        }
    }
    *Err = err ;
    *mu = C ;
    *Gtd = gtd ;
    Com->user->nz = nz ;
    if ( a_is_1 == BLGTRUE ) Com->a2sum = (BLGFLOAT) nz ;
    else                     Com->a2sum = a2sum ;
    return (0) ;
}

/* ==========================================================================
   === BLGas ================================================================
   ==========================================================================
   Evaluate the affine scaling search direction d(mu) with the property
   that a'd(mu) = 0 where

                      g_j - mu a_j
   d_j (mu) = - ----------------------------
                 lambda_k + |g_j - mu a_j|/X

   X = hi_j - x_j if g_j - mu a_j <= 0
   X = x_j - lo_j if g_j - mu a_j  > 0

   If a is NULL, take mu = 0 to evaluate d.
   We first find an interval [a, b] containing mu:

           A = min_j g_j/a_j
           B = max_j g_j/a_j

    A trial mu is obtained by treating the denominator as constant and
    solving the linear equation a'd(mu) = 0 for mu = a'g/a'a.
    If mu is outside the interval [A, B], then we replace mu by the
    interval midpoint. Next we alternate between a Newton step and
    a secant step. If the pair of steps do not reduce the interval
    width by gamma*(B-A), then we perform a bisection step.

                                   a_j^2 lambda_k
    Note that a_j d_j'(mu) = --------------------------- .
                              (lambda_k - |g - mu a|/X)^2

    In the code which follows, the sign of d_j has been reversed.
   ========================================================================== */

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
)
{
    int iu, j, PrintLevel, method, max_mu_its, prior_side,
        a_is_1, bnd_is_scalar, ok ;
    BLGFLOAT A, B, C, D, r, t, u, alk, aj, da, db, dc, dd, adsum,
             va, vb, vc, vd, gj, X, Y, Clo, Chi, gamma, tol0, tol1, width ;
    BLGparm *Parm ;

    Parm = Com->Parm ;
    a_is_1 = Parm->a_is_1 ;
    bnd_is_scalar = Parm->bnd_is_scalar ;

    if ( bnd_is_scalar == BLGTRUE )
    {
        Clo = *lo ;
        Chi = *hi ;
    }
    /* compute affine scaling d */
    gamma = Parm->gamma ;
    max_mu_its = Parm->max_mu_its ;
    PrintLevel = Parm->PrintLevel ;

    /* bracket the root */
    A =  BLGINF ; /* stores lower bound for bracketing interval */
    B = -BLGINF ; /* stores upper bound for bracketing interval */
    if ( a_is_1 == BLGTRUE )
    {
        for (j = 0; j < n; j++)
        {
            gj = g [j] ;
            if ( gj < A ) A = gj ;
            if ( gj > B ) B = gj ;
        }
    }
    else
    {
        for (j = 0; j < n; j++)
        {
            gj = g [j] ;
            aj = a [j] ;
            if ( aj != BLGZERO )
            {
                t = gj/aj ;
                if ( t < A ) A = t ;
                if ( t > B ) B = t ;
            }
        }
    }

    C = -*c ;
    /* if the estimate lies outside bracketing interval, use the midpoint */
    if ( (C <= A) || (C >= B) ) C = .5*(A+B) ;

    if ( PrintLevel >= 2 )
    {
        printf ("A: %e C: %e B: %e\n", A, C, B) ;
    }

    /* evaluate function and derivative at a, b, and c */
    va = BLGZERO ;
    vb = BLGZERO ;
    vc = BLGZERO ;
    da = BLGZERO ;
    db = BLGZERO ;
    dc = BLGZERO ;
    /* The goal is to obtain a'd = 0. Based on the starting guess for mu,
       we compute sum |a_j d_j| and use this in the stopping condition */
    adsum = BLGZERO ;

    if ( a_is_1 == BLGTRUE )
    {
        for (j = 0; j < n; j++)
        {
            gj = g [j] ;
            t = gj - A ;
            if ( bnd_is_scalar == BLGTRUE )
            {
                X = x [j] - Clo ;
                Y = Chi - x [j] ;
            }
            else
            {
                X = x [j] - lo [j] ;
                Y = hi [j] - x [j] ;
            }
            if ( t > BLGZERO )
            {
                if ( X > BLGZERO )
                {
                    r = t/X ;
                    u = lambda_k + r ;
                    va += t/u ;
                    da -= lambda_k/(u*u);
                }
            }
            else
            {
                if ( Y > BLGZERO )
                {
                    r = t/Y ;
                    u = lambda_k - r ;
                    va += t/u ;
                    da -= lambda_k/(u*u) ;
                }
            }
            t = gj - B ;
            if ( t > BLGZERO )
            {
                if ( X > BLGZERO )
                {
                    r = t/X ;
                    u = lambda_k + r ;
                    vb += t/u ;
                    db -= lambda_k/(u*u) ;
                }
            }
            else
            {
                if ( Y > BLGZERO )
                {
                    r = t/Y ;
                    u = lambda_k - r ;
                    vb += t/u ;
                    db -= lambda_k/(u*u) ;
                }
            }
            t = gj - C ;
            if ( t > BLGZERO )
            {
                if ( X > BLGZERO )
                {
                    r = t/X ;
                    u = lambda_k + r ;
                    r = t/u ;
                    vc += r ;
                    adsum += fabs(r) ;
                    dc -= lambda_k/(u*u) ;
                }
            }
            else
            {
                if ( Y > BLGZERO )
                {
                    r = t/Y ;
                    u = lambda_k - r ;
                    r = t/u ;
                    vc += r ;
                    adsum += fabs(r) ;
                    dc -= lambda_k/(u*u) ;
                }
            }
        }
    }
    else
    {
        for (j = 0; j < n; j++)
        {
            aj = a [j] ;
            gj = g [j] ;
            t = gj - A*aj ;
            alk = aj*aj*lambda_k ;
            if ( bnd_is_scalar == BLGTRUE )
            {
                X = x [j] - Clo ;
                Y = Chi - x [j] ;
            }
            else
            {
                X = x [j] - lo [j] ;
                Y = hi [j] - x [j] ;
            }
            if ( t > BLGZERO )
            {
                if ( X > BLGZERO )
                {
                    r = t/X ;
                    u = lambda_k + r ;
                    va += aj*t/u ;
                    da -= alk/(u*u);
                }
            }
            else
            {
                if ( Y > BLGZERO )
                {
                    r = t/Y ;
                    u = lambda_k - r ;
                    va += aj*t/u ;
                    da -= alk/(u*u) ;
                }
            }
            t = gj - B*aj ;
            if ( t > BLGZERO )
            {
                if ( X > BLGZERO )
                {
                    r = t/X ;
                    u = lambda_k + r ;
                    vb += aj*t/u ;
                    db -= alk/(u*u) ;
                }
            }
            else
            {
                if ( Y > BLGZERO )
                {
                    r = t/Y ;
                    u = lambda_k - r ;
                    vb += aj*t/u ;
                    db -= alk/(u*u) ;
                }
            }
            t = gj - C*aj ;
            if ( t > BLGZERO )
            {
                if ( X > BLGZERO )
                {
                    r = t/X ;
                    u = lambda_k + r ;
                    r = aj*t/u ;
                    vc += r ;
                    adsum += fabs(r) ;
                    dc -= alk/(u*u) ;
                }
            }
            else
            {
                if ( Y > BLGZERO )
                {
                    r = t/Y ;
                    u = lambda_k - r ;
                    r = aj*t/u ;
                    vc += r ;
                    adsum += fabs(r) ;
                    dc -= alk/(u*u) ;
                }
            }
        }
    }

    if ( (va < BLGZERO) || (vb > BLGZERO) )
    {
        printf ("computation of muk: either va = %e < 0 or vb = %e > 0\n",
                 va, vb) ;
        return (4) ;
    }

    if ( vc > BLGZERO ) /* if vc > 0 move right */
    {
        A = C ;
        va = vc ;
        da = dc ;
    }
    else             /* move left */
    {
        B = C ;
        vb = vc ;
        db = dc ;
    }

    method = 0 ;                             /* start with Newton's method */
    width = BLGINF ;
    iu = 0 ;
    if ( va <= BLGZERO )
    {
        C = A ;
        goto Exit ;
    }
    else if ( vb >= BLGZERO )
    {
        C = B ;
        goto Exit ;
    }
    prior_side = 0 ;
    ok = 1 ;
    tol0 = adsum*Parm->epsmu0 ;
    for (iu = 0; iu < max_mu_its; iu++)
    {
        /* test for convergence */
        tol1 = Parm->epsmu1*BLGMAX (fabs (A), fabs (B)) ;
        t = B - A ;
        if ( ((t <= tol1) && ok) || (fabs(va) <= tol0) || (fabs(vb) <= tol0) )
        {
            if ( iu > 0 )
            {
                if ( va > -vb )
                {
                    C = B ;
                    vc = vb ;
                }
                else
                {
                    C = A ;
                    vc = va ;
                }
                goto Exit ;
            }
        }

        ok = 1 ;
        /* if the previous method was secant and b-a > width, use bisection */
        if ( (method == 0) && (B-A > width) )
        {
            C = .5*(A+B) ;
            width = BLGINF ;                  /* ensures Newton is used next */
            prior_side = 0 ;
        }
        else if ( va > -vb )
        {
            if ( method == 0 )               /* Newton */
            {
                width = gamma*(B-A) ;        /* Secant must give b-a <= width */
                C = B - vb/db ;              /* Newton from b */
                /* bisect if Newton outside [A,B]*/
                if ( (C >= B) || (C <= A) ) C = .5*(A+B) ;
                method = 1 ;                 /* switch to Secant */
            }
            else                             /* Secant */
            {
                if ( iu == 1 )
                {
                    if ( va/vb <= -1.e3 )
                    {
                        C = B - 2*vb/db ;
                        if ( C <= A ) C = B - (A-B)*(vb/(va-vb)) ;
                    }
                    else C = B - (A-B)*(vb/(va-vb)) ; /*Secant from b*/
                }
                else
                {
                    C = B - (A-B)*(vb/(va-vb)) ; /* Secant from b */
                }
                method = 0 ;                 /* switch to Newton */
            }
        }
        else
        {
            if ( method == 0 )               /* Newton */
            {
                width = gamma*(B-A) ;        /* Secant must give b-a <= width */
                C = A - va/da ;              /* Newton from a */
                /* bisect if Newton outside [A,B]*/
                if ( (C >= B) || (C <= A) ) C = .5*(A+B) ;
                method = 1 ;                 /* switch to Secant */
            }
            else                             /* Secant */
            {
                if ( iu == 1 )
                {
                    if ( vb/va <= -1.e3 )
                    {
                        C = A - 2*va/da ;
                        if ( C >= B ) C = A - (A-B)*(va/(va-vb)) ;
                    }
                    else C = A - (A-B)*(va/(va-vb)) ; /*Secant from a*/
                }
                else
                {
                    C = A - (A-B)*(va/(va-vb)) ; /* Secant from a */
                }
                method = 0 ;                 /* switch to Newton */
            }
        }

        if ( PrintLevel == 2 )
        {
            if ( width == BLGINF )
            {
                printf ("Bisect A: %12.5e va: %12.5e C: %12.5e B: "
                        "%12.5e vb: %12.5e\n", A, va, C, B, vb) ;
            }
            else if ( method == 1 )
            {
                printf ("Newton A: %12.5e va: %12.5e C: %12.5e B: "
                        "%12.5e vb: %12.5e\n", A, va, C, B, vb) ;
            }
            else
            {
                printf ("Secant A: %12.5e va: %12.5e C: %12.5e B: "
                        "%12.5e vb: %12.5e\n", A, va, C, B, vb) ;
            }
        }
        /* evaluate function and derivative */
        BLGval (&vc, &dc, C, lambda_k, x, g, a, lo, hi, n, Com) ;
        if ( vc > BLGZERO ) /* move right */
        {
            if ( C <= A ) goto Exit ;
            if ( prior_side == -1 ) /* left side updated twice */
            {
                D = BLGMIN (A + tol0, .5*(A+B)) ;
                if ( D > C )
                {
                    BLGval (&vd, &dd, D, lambda_k, x, g, a, lo, hi, n, Com);
                    if ( vd >= BLGZERO ) /* update left side */
                    {
                        A = D ;
                        va = vd ;
                        da = dd ;
                    }
                    else /* update right side and flag for Newton step */
                    {
                        prior_side = 1 ;
                        B = D ;
                        vb = vd ;
                        db = dd ;
                        A = C ;
                        va = vc ;
                        da = dc ;
                        width = BLGINF ;
                        method = 0 ;
                        ok = 0 ;
                    }
                }
                else
                {
                    A = C ;
                    va = vc ;
                    da = dc ;
                }
            }
            else
            {
                prior_side = -1 ;
                A = C ;
                va = vc ;
                da = dc ;
            }
        }
        else             /* move left */
        {
            if ( C >= B ) goto Exit ;
            if ( prior_side == 1 ) /* right side updated twice */
            {
                D = BLGMAX (B - tol0, .5*(A+B)) ;
                if ( D < C )
                {
                    BLGval (&vd, &dd, D, lambda_k, x, g, a, lo, hi, n, Com);
                    if ( vd <= BLGZERO ) /* update right side */
                    {
                        B = D ;
                        vb = vd ;
                        db = dd ;
                    }
                    else /* update left side and flag for Newton step */
                    {
                        prior_side = -1 ;
                        A = D ;
                        va = vd ;
                        da = dd ;
                        B = C ;
                        vb = vc ;
                        db = dc ;
                        width = BLGINF ;
                        method = 0 ;
                        ok = 0 ;
                    }
                }
                else
                {
                    B = C ;
                    vb = vc ;
                    db = dc ;
                }
            }
            else
            {
                prior_side = 1 ;
                B = C ;
                vb = vc ;
                db = dc ;
            }
        }
        if ( va <= BLGZERO)
        {
            C = A ;
            goto Exit ;
        }
        else if ( vb >= BLGZERO )
        {
            C = B ;
            goto Exit ;
        }
    }
    return (6) ;

    Exit:
    if ( PrintLevel >= 2 )
    {
        printf ("Mu iterations: %6i A: %12.5e C: %12.5e B: %12.5e vc: %12.5e\n",
                 iu, A, C, B, vc) ;
    }
    *c = -C ;
    Com->nmu += iu ;
    return (0) ;
}
/* ==========================================================================
   === BLGval ===============================================================
   ==========================================================================
       Evaluate r(C)  and r'(C) where r (C) = a transpose d (C)
   ========================================================================== */
PRIVATE void BLGval
(
    BLGFLOAT      *vC, /* a transpose d(C)  (returned) */
    BLGFLOAT      *dC, /* a transpose d'(C) (returned) */
    BLGFLOAT        C,
    BLGFLOAT lambda_k,
    BLGFLOAT       *x,
    BLGFLOAT       *g,
    BLGFLOAT       *a,
    BLGFLOAT      *lo,
    BLGFLOAT      *hi,
    int             n,
    BLGcom       *Com  /* common variables */
)
{
    int j, bnd_is_scalar, a_is_1 ;
    BLGFLOAT aj, alk, r, t, u, vc, dc, X, Clo, Chi, lo_cut, hi_cut ;
    BLGparm *Parm ;
    Parm = Com->Parm ;
    bnd_is_scalar = Parm->bnd_is_scalar ;
    a_is_1 = Parm->a_is_1 ;
    vc = BLGZERO ;
    dc = BLGZERO ;
    if ( bnd_is_scalar == BLGTRUE )
    {
        Clo = *lo ;
        Chi = *hi ;
    }
    lo_cut = Parm->lo_cut ;
    hi_cut = Parm->hi_cut ;
    if ( a_is_1 == BLGTRUE )
    {
        for (j = 0; j < n; j++)
        {
            t = g [j] - C ;
            if ( t > BLGZERO )
            {
                if ( bnd_is_scalar == BLGTRUE ) X = x [j] - Clo ;
                else                            X = x [j] - lo [j] ;
                if ( X > lo_cut )
                {
                    r = t/X ;
                    u = lambda_k + r ;
                    vc += t/u ;
                    dc -= lambda_k/(u*u) ;
                }
            }
            else
            {
                if ( bnd_is_scalar == BLGTRUE ) X = Chi - x [j] ;
                else                            X = hi [j] - x [j] ;
                if ( X > hi_cut )
                {
                    r = t/X ;
                    u = lambda_k - r ;
                    vc += t/u ;
                    dc -= lambda_k/(u*u) ;
                }
            }
        }
    }
    else
    {
        for (j = 0; j < n; j++)
        {
            aj = a [j] ;
            t = g [j] - C*aj ;
            alk = aj*aj*lambda_k ;
            if ( t > BLGZERO )
            {
                if ( bnd_is_scalar == BLGTRUE ) X = x [j] - Clo ;
                else                            X = x [j] - lo [j] ;
                if ( X > lo_cut )
                {
                    r = t/X ;
                    u = lambda_k + r ;
                    vc += aj*t/u ;
                    dc -= alk/(u*u) ;
                }
            }
            else
            {
                if ( bnd_is_scalar == BLGTRUE ) X = Chi - x [j] ;
                else                            X = hi [j] - x [j] ;
                if ( X > hi_cut )
                {
                    r = t/X ;
                    u = lambda_k - r ;
                    vc += aj*t/u ;
                    dc -= alk/(u*u) ;
                }
            }
        }
    }
    *vC = vc ;
    *dC = dc ;
}

/* ==========================================================================
   === BLGevaluate ==========================================================
   ==========================================================================
   Evaluate function, gradient, or both
   ========================================================================== */
PRIVATE void BLGevaluate
(
    BLGcom *Com,  /* BLGcom structure */
    char  *what   /* f = function, g = gradient, fg = function and gradient */
)
{
    BLGobjective *user ;
    user = Com->user ;
    if ( !strcmp (what, "fg") )     /* compute function and gradient */
    {
        Com->ng++ ;
        if ( Com->Parm->QuadCost == BLGTRUE )
        {
            Com->Ad (user) ;

            /* update the function value */
            if ( user->nz == user->n )
                   user->f = user->fold
                           + BLGdot (user->d, user->gold, user->n)
                           + .5*BLGdot (user->d, user->p, user->n) ;
            else   user->f = user->fold
                           + BLGdoti (user->d, user->i, user->gold, user->nz)
                           + .5*BLGdoti (user->d, user->i, user->p, user->nz) ;
            /* update the gradient */
            BLGstep (user->g, user->gold, user->p, BLGONE, Com->n) ;
        }
        else
        {
            Com->nf++ ;
            user->nz = user->n ;
            if ( Com->valgrad != NULL ) Com->valgrad (user) ;
            else
            {
                Com->value (user) ;
                Com->grad (user) ;
            }
        }
    }
    else if ( !strcmp (what, "f") ) /* compute function */
    {
        if ( Com->Parm->QuadCost == BLGTRUE )
        {
            Com->ng++ ;
            Com->Ad (user) ;
            if ( user->nz == user->n )
                 Com->deriv2 = BLGdot  (user->d, user->p, user->n) ; /* d'Ad */
            else Com->deriv2 = BLGdoti (user->d, user->i, user->p, user->nz) ;
            user->f = user->fold + user->deriv1 + .5*Com->deriv2 ;
        }
        else
        {
            Com->nf++ ;
            Com->value (user) ;
            if ( Com->valgrad == Com->value ) Com->ng++ ;
        }
    }
    else                                 /* compute gradient */
    {
        Com->ng++ ;
        Com->grad (user) ;
    }
}

/* ==========================================================================
   === BLGbb ================================================================
   ==========================================================================
       Compute lambda_k using BB formula.  If the returned value of t is
       positive, then BB the parameter was updated, otherwise set t = -1
   ========================================================================== */
PRIVATE BLGFLOAT BLGbb /* return lambda_k */
(
    BLGFLOAT alpha,
    int         mm, /* stop current cycle when Parm->MaxCycle - mm <= 0 */
    int          n, /* problem dimension */
    BLGcom    *Com  /* common variables */
)
{
    int j, k, nz, *i ;
    BLGFLOAT lambda0, dtd, dk, sy, ss, normg, normx,
            t, *x, *d, *g, *xold, *gold ;
    BLGobjective *Objective ;
    BLGparm *Parm ;
    Parm = Com->Parm ;
    lambda0 = Parm->lambda0 ;
    Objective = Com->user ;
    d = Objective->d ;
    nz= Objective->nz ;
    i = Objective->i ;
    x = Objective->x ;
    g = Objective->g ;
    xold = Objective->xold ;
    gold = Objective->gold ;
    sy = BLGZERO ;
    ss = BLGZERO ;
    dtd= BLGZERO ;
    if ( Parm->QuadCost == BLGTRUE )
    {
        sy = Com->deriv2 ; /* d'Ad (the alpha's cancel) */
        if ( alpha == BLGONE )
        {
            for (k = 0; k < nz; k++)
            {
                dk = d [k] ;
                ss += dk*dk ;
                j = i [k] ;
                xold [j] = x [j] ;
            }
        }
        else
        {
            for (k = 0; k < nz; k++)
            {
                dk = d [k] ;
                ss += dk*dk ;
                j = i [k] ;
                xold [j] += alpha*dk ;
            }
        }
        /* update g */
        BLGstep (gold, gold, Objective->p, alpha, n) ;
    }
    else
    {
        for (k = 0; k < nz; k++)
        {
            dk = d [k] ;
            ss += dk*dk ;
            j = i [k] ;
            xold [j] = x[j] ;
            sy += dk*(g [j] - gold [j]) ;
        }
        ss *= alpha ;
        BLGcopy (gold, g, n) ;
    }
    t = -BLGONE ;
    if ( sy <= BLGZERO || ss == BLGZERO )
    {
       if ( Parm->QuadCost == BLGTRUE )
       {
          normg = BLGmax (gold, n) ;
          normx = BLGmax (xold, n) ;
          if ( normx > BLGZERO )
          {
             t= BLGMAX (Parm->lambda0Factor*normg/normx, lambda0) ;
          }
          else t = lambda0 ;
       }
       else if ( Parm->MaxCycle - mm <= 0 )
       {
          if ( Parm->QuadCost == BLGTRUE )
          {
              normg = BLGmax (gold, n) ;
              normx = BLGmax (xold, n) ;
          }
          else
          {
              normg = BLGmax (g, n) ;
              normx = BLGmax (x, n) ;
          }
          if ( normx > BLGZERO )
          {
              t = BLGMAX (Parm->lambda0Factor*normg/normx, lambda0) ;
          }
          else t = lambda0 ;
       }
    }
    else
    {
       t = BLGMAX (sy/ss, lambda0) ;
    }
    return (t) ;
}

/* =========================================================================
   === BLGgather ===========================================================
   =========================================================================
   compress x by extracting the subspace elements and copying them to
   the compressed array y
   ========================================================================= */
PRIVATE void BLGgather
(
    BLGFLOAT *y, /* compressed vector */
    BLGFLOAT *x, /* current x */
    int      *i, /* indices associated with subspace */
    int       n  /* number of elements to extract from x */
)
{
    int k, n5, *ik ;
    BLGFLOAT *yk ;
    n5 = n % 5 ;
    for (k = 0; k < n5; k++) y [k] = x [i [k]] ;
    yk = y+k ;
    ik = i+k ;
    for (; k < n; k += 5)
    {
        *(yk++) = x [*(ik++)] ;
        *(yk++) = x [*(ik++)] ;
        *(yk++) = x [*(ik++)] ;
        *(yk++) = x [*(ik++)] ;
        *(yk++) = x [*(ik++)] ;
    }
}

/* =========================================================================
   === BLGscatter ==========================================================
   =========================================================================
   scatter elements from the compressed array y into x
   ========================================================================= */
PRIVATE void BLGscatter
(
    BLGFLOAT *x, /* scatter y into x */
    BLGFLOAT *y, /* compressed vector */
    int      *i, /* indices associated with subspace */
    int       n  /* number of elements to extract from x */
)
{
    int k, n5, *ik ;
    BLGFLOAT *yk ;
    n5 = n % 5 ;
    for (k = 0; k < n5; k++) x [i [k]] = y [k] ;
    yk = y+k ;
    ik = i+k ;
    for (; k < n; k += 5)
    {
        x [*(ik++)] = *(yk++) ;
        x [*(ik++)] = *(yk++) ;
        x [*(ik++)] = *(yk++) ;
        x [*(ik++)] = *(yk++) ;
        x [*(ik++)] = *(yk++) ;
    }
}

/* =========================================================================
   === BLGstep =============================================================
   =========================================================================
   set xnew = x + st*d
   ========================================================================= */
PRIVATE void BLGstep
(
    BLGFLOAT *xnew, /* updated x vector */
    BLGFLOAT    *x, /* current x */
    BLGFLOAT    *d, /* search direction */
    BLGFLOAT    st, /* stepsize */
    int          n  /* dimension */
)
{
    int j, n5 ;
    n5 = n % 5 ;
    if ( st == BLGONE )
    {
        for (j = 0; j < n5; j++) xnew [j] = x [j] + d [j] ;
        for (; j < n; j += 5)
        {
            xnew [j  ] = x [j  ] + d [j  ] ;
            xnew [j+1] = x [j+1] + d [j+1] ;
            xnew [j+2] = x [j+2] + d [j+2] ;
            xnew [j+3] = x [j+3] + d [j+3] ;
            xnew [j+4] = x [j+4] + d [j+4] ;
        }
    }
    else if ( st == -BLGONE )
    {
        for (j = 0; j < n5; j++) xnew [j] = x [j] - d [j] ;
        for (; j < n; j += 5)
        {
            xnew [j  ] = x [j  ] - d [j  ] ;
            xnew [j+1] = x [j+1] - d [j+1] ;
            xnew [j+2] = x [j+2] - d [j+2] ;
            xnew [j+3] = x [j+3] - d [j+3] ;
            xnew [j+4] = x [j+4] - d [j+4] ;
        }
    }
    else
    {
        for (j = 0; j < n5; j++)
        {
            xnew [j] = x [j] + st*d [j] ;
        }
        for (; j < n; j += 5)
        {
            xnew [j  ] = x [j  ] + st*d [j  ] ;
            xnew [j+1] = x [j+1] + st*d [j+1] ;
            xnew [j+2] = x [j+2] + st*d [j+2] ;
            xnew [j+3] = x [j+3] + st*d [j+3] ;
            xnew [j+4] = x [j+4] + st*d [j+4] ;
        }
    }
}

/* =========================================================================
   === BLGstepi ============================================================
   =========================================================================
   set xnew_j = x_j + st*d_j, j = i [k], 0 <= k < nz
   ========================================================================= */
PRIVATE void BLGstepi
(
    BLGFLOAT *xnew, /* updated x vector */
    int         *i, /* indices of xnew that change */
    BLGFLOAT    *x, /* current x */
    BLGFLOAT    *d, /* search direction */
    BLGFLOAT    st, /* stepsize */
    int         nz  /* number of components to be updated */
)
{
    int j, k, n5, *jk ;
    n5 = nz % 5 ;
    jk = i ;
    if ( st == BLGONE )
    {
        for (k = 0; k < n5; k++)
        {
            j = *(jk++) ;
            xnew [j] = x [j] + d [k] ;
        }
        for (; k < nz; )
        {
            j = *(jk++) ;
            xnew [j] = x [j] + d [k++] ;
            j = *(jk++) ;
            xnew [j] = x [j] + d [k++] ;
            j = *(jk++) ;
            xnew [j] = x [j] + d [k++] ;
            j = *(jk++) ;
            xnew [j] = x [j] + d [k++] ;
            j = *(jk++) ;
            xnew [j] = x [j] + d [k++] ;
        }
    }
    else
    {
        for (k = 0; k < n5; )
        {
            j = *(jk++) ;
            xnew [j] = x [j] + st*d [k++] ;
        }
        for (; k < nz; )
        {
            j = *(jk++) ;
            xnew [j] = x [j] + st*d [k++] ;
            j = *(jk++) ;
            xnew [j] = x [j] + st*d [k++] ;
            j = *(jk++) ;
            xnew [j] = x [j] + st*d [k++] ;
            j = *(jk++) ;
            xnew [j] = x [j] + st*d [k++] ;
            j = *(jk++) ;
            xnew [j] = x [j] + st*d [k++] ;
        }
    }
}

/* =========================================================================
   === BLGdot ==============================================================
   =========================================================================
   Compute dot product of x and y
   ========================================================================= */
PRIVATE BLGFLOAT BLGdot
(
    BLGFLOAT *x, /* first vector */
    BLGFLOAT *y, /* second vector */
    int       n  /* length of vectors */
)
{
    int j, n5 ;
    BLGFLOAT t ;
    t = BLGZERO ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) t += x [j]*y [j] ;
    for (; j < n; j += 5)
    {
        t += x [j]*y[j] + x [j+1]*y [j+1] + x [j+2]*y [j+2]
                        + x [j+3]*y [j+3] + x [j+4]*y [j+4] ;
    }
    return (t) ;
}

/* =========================================================================
   === BLGdoti =============================================================
   =========================================================================
   Compute dot product of x and y, x is sparse (only store nonzeros)
   and y is dense
   ========================================================================= */
PRIVATE BLGFLOAT BLGdoti
(
    BLGFLOAT *x, /* first vector (sparse, nonzeros only) */
    int      *i, /* indices in y associated with elements of x */
    BLGFLOAT *y, /* second vector (dense) */
    int       n  /* length of x */
)
{
    int j, n5, *ij ;
    BLGFLOAT t, *xj ;
    t = BLGZERO ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) t += x [j]*y [i [j]] ;
    for (; j < n; j += 5)
    {
        ij = i+j ;
        xj = x+j ;
        t += xj[0]*y[ij[0]] + xj[1]*y[ij[1]] + xj[2]*y[ij[2]]
                            + xj[3]*y[ij[3]] + xj[4]*y[ij[4]] ;
    }
    return (t) ;
}

/* =========================================================================
   === BLGgtd ==============================================================
   =========================================================================
   Compute gtd = (g-t*a)'d = g'd
   ========================================================================= */
PRIVATE BLGFLOAT BLGgtd
(
    BLGFLOAT *g, /* gradient, dense */
    BLGFLOAT *a, /* linear constraint vector, dense */
    BLGFLOAT *d, /* search direction, sparse */
    BLGFLOAT  t, /* t chosen so that a'(g-ta) = 0 */
    int      *i, /* indices in g and a associated with indices in d */
    int      nz, /* length of d */
    int  a_is_1  /* TRUE if a is identically 1 */
)
{
    int j, n5 ;
    BLGFLOAT s ;
    s = BLGZERO ;
    n5 = nz % 5 ;
    if ( a_is_1 == BLGTRUE )
    {
        for (j = 0; j < n5; j++) s +=  (g [i [j]] - t)*d [j] ;
        for (; j < nz; j += 5)
        {
            s += (g [i [j  ]] - t)*d [j]
              +  (g [i [j+1]] - t)*d [j+1]
              +  (g [i [j+2]] - t)*d [j+2]
              +  (g [i [j+3]] - t)*d [j+3]
              +  (g [i [j+4]] - t)*d [j+4] ;
        }
    }
    else
    {
        for (j = 0; j < n5; j++) s +=  (g [i [j]] - a [i [j]]*t)*d [j] ;
        for (; j < nz; j += 5)
        {
            s += (g [i [j  ]] - a [i [j]]  *t)*d [j]
              +  (g [i [j+1]] - a [i [j+1]]*t)*d [j+1]
              +  (g [i [j+2]] - a [i [j+2]]*t)*d [j+2]
              +  (g [i [j+3]] - a [i [j+3]]*t)*d [j+3]
              +  (g [i [j+4]] - a [i [j+4]]*t)*d [j+4] ;
        }
    }
    return (s) ;
}

/* =========================================================================
   === BLGdphi =============================================================
   =========================================================================
   Compute gtd = (g-t*a)'d = g'd where t = a'g/a'a
   ========================================================================= */
PRIVATE BLGFLOAT BLGdphi
(
    BLGFLOAT *g, /* gradient, dense */
    BLGFLOAT *d, /* search direction, sparse */
    BLGcom *Com  /* common variables */
)
{
    int a_is_1, j, n5, nz, *i ;
    BLGFLOAT gtd, t, *a ;
    a = Com->a ;
    a_is_1 = Com->Parm->a_is_1 ;;
    nz = Com->user->nz ;
    i = Com->user->i ;

    /* if there is no linear constraint, then t = 0 */
    if ( (a == NULL) && (a_is_1 == BLGFALSE) ) return (BLGdoti (d, i, g, nz)) ;

    t = BLGZERO ;
    n5 = nz % 5 ;
    if ( a_is_1 == BLGTRUE )
    {
        for (j = 0; j < n5; j++) t += g [i [j]] ;
        for (; j < nz; j += 5)
        {
            t += g [i[j]] + g [i[j+1]] + g [i[j+2]] + g [i[j+3]] + g [i[j+4]] ;
        }

    }
    else
    {
        for (j = 0; j < n5; j++) t += g [i [j]]*a [i [j]] ;
        for (; j < nz; j += 5)
        {
            t += g [i[j]]  *a [i[j]]   + g [i[j+1]]*a [i[j+1]]
              +  g [i[j+2]]*a [i[j+2]] + g [i[j+3]]*a [i[j+3]]
              +  g [i[j+4]]*a [i[j+4]] ;
        }
    }
    if ( Com->a2sum != BLGZERO ) t = t/Com->a2sum ;
    gtd = BLGgtd (g, a, d, t, i, nz, a_is_1) ;
    return (gtd) ;
}

/* =========================================================================
   === BLGadd ==============================================================
   =========================================================================
   add the components of a vector
   ========================================================================= */
double BLGadd
(
    BLGFLOAT *x, /* vector */
    int       n  /* length of vector */
)
{
    int j, n5 ;
    BLGFLOAT t ;

    t = BLGZERO ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) t += x [j] ;
    for (; j < n; j += 5)
    {
        t += x [j] + x [j+1] + x [j+2] + x [j+3] + x [j+4] ;
    }
    return (t) ;
}

/* =========================================================================
   === BLGcopy ==============================================================
   =========================================================================
   Copy the second vector to the first x = y
   ========================================================================= */
PRIVATE void BLGcopy
(
    BLGFLOAT *x, /* copy of y */
    BLGFLOAT *y, /* given vector */
    int       n  /* length of vectors */
)
{
    int j, n5 ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = y [j] ;
    for (; j < n; j += 5)
    {
        x [j]   =  y[j] ;
        x [j+1] =  y[j+1] ;
        x [j+2] =  y[j+2] ;
        x [j+3] =  y[j+3] ;
        x [j+4] =  y[j+4] ;
    }
}

/* =========================================================================
   === BLGicopy ============================================================
   =========================================================================
   Copy the second vector to the first x = y
   ========================================================================= */
PRIVATE void BLGicopy
(
    int *x, /* copy of y */
    int *y, /* given vector */
    int  n  /* length of vectors */
)
{
    int j, n5 ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = y [j] ;
    for (; j < n; j += 5)
    {
        x [j]   =  y[j] ;
        x [j+1] =  y[j+1] ;
        x [j+2] =  y[j+2] ;
        x [j+3] =  y[j+3] ;
        x [j+4] =  y[j+4] ;
    }
}

/* =========================================================================
   === BLGmax ==============================================================
   =========================================================================
   Return max {fabs (x [j]) : 1 <= j < n}
   ========================================================================= */
BLGFLOAT BLGmax
(
    BLGFLOAT *x,
    int       n
)
{
    BLGFLOAT xnorm ;
    int j, n5 ;
    n5 = n % 5 ;
    xnorm = BLGZERO ;
    for (j = 0; j < n5; j++) if ( xnorm < fabs (x [j]) ) xnorm = fabs (x [j]) ;
    for (; j < n; j += 5)
    {
        if ( xnorm < fabs (x [j]  ) ) xnorm = fabs (x [j]) ;
        if ( xnorm < fabs (x [j+1]) ) xnorm = fabs (x [j+1]) ;
        if ( xnorm < fabs (x [j+2]) ) xnorm = fabs (x [j+2]) ;
        if ( xnorm < fabs (x [j+3]) ) xnorm = fabs (x [j+3]) ;
        if ( xnorm < fabs (x [j+4]) ) xnorm = fabs (x [j+4]) ;
    }
    return (xnorm) ;
}

/* =========================================================================
   === BLGint_init =========================================================
   =========================================================================
   Initialize an array
   ========================================================================= */
void BLGint_init
(
    int *x,  /* array to be initialized */
    int  s,  /* scalar */
    int  n   /* length of x */
)
{
    int j, n5 ;
    int *xj ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

/* ==========================================================================
   === BLGindexMinSort ======================================================
   ==========================================================================
   Sort an array x in increasing order (from LPDASA).
   The indices giving increasing order are stored in y
   while x is not changed
   ========================================================================== */

PRIVATE void BLGindexMinSort
(
    int      *y,  /* n-by-1, indices of x giving decreasing order, output */
    int      *w,  /* n-by-1 work array, input modified */
    BLGFLOAT *x,  /* n-by-1, input not modified */
    int       n   /* number of elements to sort */
)
{
    int i, j, k, l, m, p, q, *yi, *wi ;
    BLGFLOAT s, t ;

    if ( n == 0 ) return ;
    y [0] = 0 ;
    if ( n < 2 ) return ;
    if ( n < 3 )
    {
        if ( x [0] > x [1] )
        {
            y [0] = 1 ;
            y [1] = 0 ;
        }
        else y [1] = 1 ;
        return ;
    }

    j = k = 0 ;
    for (i = 1; i < n; i++)
    {
        if ( x [i] < x [j] )
        {
            w [k] = i ;
            k = i ;
        }
        y [i] = j = i ;
    }

    w [k] = n ;
    while ( k > 0 )
    {
        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = y [i] ;
            s = x [p] ;
            j = w [i] ;
            k = j ;
            if ( j == n )
            {
                y [i] = j ;
                l = j ;
                w [m] = p ;
                k += m - i ;
                yi = y+(i-m) ;
                for (m++; m < k; m++) w [m] = yi [m] ;
            }
            else
            {
                q = y [j] ;
                t = x [q] ;
                l = w [j] ;
                y [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        w [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            w [m] = p ;
                            k += m - i ;
                            yi = y+(i-m) ;
                            for (m++; m < k; m++) w [m] = yi [m] ;
                            break ;
                        }
                        q = y [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        w [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            w [m] = q ;
                            k = m + l - j ;
                            yi = y+(j-m) ;
                            for (m++; m < k; m++) w [m] = yi [m] ;
                            break ;
                        }
                        p = y [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n )
        {
            for (i = 0; i < n; i++) y [i] = w [i] ;
            return ;
        }

        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = w [i] ;
            s = x [p] ;
            j = y [i] ;
            k = j ;
            if ( j == n )
            {
                w [i] = j ;
                l = j ;
                y [m] = p ;
                k += m - i ;
                wi = w+(i-m) ;
                for (m++; m < k; m++) y [m] = wi [m] ;
            }
            else
            {
                q = w [j] ;
                t = x [q] ;
                l = y [j] ;
                w [i] = l ;
                while ( 1 )
                {
                    if ( s > t )
                    {
                        y [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            y [m] = p ;
                            k += m - i ;
                            wi = w+(i-m) ;
                            for (m++; m < k; m++) y [m] = wi [m] ;
                            break ;
                        }
                        q = w [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        y [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            y [m] = q ;
                            k = m + l - j ;
                            wi = w+(j-m) ;
                            for (m++; m < k; m++) y [m] = wi [m] ;
                            break ;
                        }
                        p = w [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n ) return ;
    }
}

/* Sort a BLGFLOAT array x in decreasing order (from LPDASA).
   The indices giving decreasing order are stored in y
   while x is not changed */

PRIVATE void BLGindexMaxSort
(
    int      *y,  /* n-by-1, indices of x giving decreasing order, output */
    int      *w,  /* n-by-1 work array, input modified */
    BLGFLOAT *x,  /* n-by-1, input not modified */
    int       n   /* number of elements to sort */
)
{
    int i, j, k, l, m, p, q, *yi, *wi ;
    BLGFLOAT s, t ;

    y [0] = 0 ;
    if ( n < 2 ) return ;
    if ( n < 3 )
    {
        if ( x [0] < x [1] )
        {
            y [0] = 1 ;
            y [1] = 0 ;
        }
        else y [1] = 1 ;
        return ;
    }

    j = k = 0 ;
    for (i = 1; i < n; i++)
    {
        if ( x [i] > x [j] )
        {
            w [k] = i ;
            k = i ;
        }
        y [i] = j = i ;
    }

    w [k] = n ;
    while ( k > 0 )
    {
        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = y [i] ;
            s = x [p] ;
            j = w [i] ;
            k = j ;
            if ( j == n )
            {
                y [i] = j ;
                l = j ;
                w [m] = p ;
                k += m - i ;
                yi = y+(i-m) ;
                for (m++; m < k; m++) w [m] = yi [m] ;
            }
            else
            {
                q = y [j] ;
                t = x [q] ;
                l = w [j] ;
                y [i] = l ;
                while ( 1 )
                {
                    if ( s < t )
                    {
                        w [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            w [m] = p ;
                            k += m - i ;
                            yi = y+(i-m) ;
                            for (m++; m < k; m++) w [m] = yi [m] ;
                            break ;
                        }
                        q = y [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        w [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            w [m] = q ;
                            k = m + l - j ;
                            yi = y+(j-m) ;
                            for (m++; m < k; m++) w [m] = yi [m] ;
                            break ;
                        }
                        p = y [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n )
        {
            for (i = 0; i < n; i++) y [i] = w [i] ;
            return ;
        }

        l = m = 0 ;
        while ( l < n )
        {
            i = l ;
            p = w [i] ;
            s = x [p] ;
            j = y [i] ;
            k = j ;
            if ( j == n )
            {
                w [i] = j ;
                l = j ;
                y [m] = p ;
                k += m - i ;
                wi = w+(i-m) ;
                for (m++; m < k; m++) y [m] = wi [m] ;
            }
            else
            {
                q = w [j] ;
                t = x [q] ;
                l = y [j] ;
                w [i] = l ;
                while ( 1 )
                {
                    if ( s < t )
                    {
                        y [m] = q ;
                        m++ ;
                        j++ ;
                        if ( j == l )
                        {
                            y [m] = p ;
                            k += m - i ;
                            wi = w+(i-m) ;
                            for (m++; m < k; m++) y [m] = wi [m] ;
                            break ;
                        }
                        q = w [j] ;
                        t = x [q] ;
                    }
                    else
                    {
                        y [m] = p ;
                        m++ ;
                        i++ ;
                        if ( i == k )
                        {
                            y [m] = q ;
                            k = m + l - j ;
                            wi = w+(j-m) ;
                            for (m++; m < k; m++) y [m] = wi [m] ;
                            break ;
                        }
                        p = w [i] ;
                        s = x [p] ;
                    }
                }
            }
        }
        if ( y [0] == n ) return ;
    }
}

/* ==========================================================================
   === BLGminSort ======================================================
   ==========================================================================
   Sort an array x in increasing order.
   The indices giving increasing order are stored in y
   while x is not changed
   ========================================================================== */
PRIVATE void BLGminSort
(
    int *x,  /* numbers to sort (input), sorted numbers (output) */
    int *y,  /* work array of length n */
    int *z,  /* work array of length n+1 */
    int  n   /* number of elements to sort */
)
{
    int i, i1, i2, i1stop, i2stop, xstart, ystart, ix, iy, l,
        *xm, *xp, *yp, t ;
    /* find increasing sequences in x and store start points in z */
    if ( n <= 1 ) return ;
    if ( n == 2 )
    {
        if ( x [1] < x [0] )
        {
            t = x [0] ;
            x [0] = x [1] ;
            x [1] = t ;
            return ;
        }
    }

    xm = x-1 ;
    z [0] = 0 ;
    xstart = 1 ;
    for (i = 1; i < n; i++)
    {
        if ( xm [i] > x [i] )
        {
            z [xstart] = i ;
            xstart++ ;
        }
    }
    z [xstart] = n ;

    while ( 1 )
    {
        /* shuffle x into y, save break points */
        if ( xstart == 1 ) return ;      /* done */
        ystart = 0 ;                     /* new start */
        iy = 0 ;
        for (l = 0; l < xstart; )
        {
            z [ystart] = iy ;
            ystart++ ;
            i1 = z [l] ;
            l++ ;
            if ( l == xstart )           /* copy rest of x into y */
            {
                for (; i1 < n; i1++) y [i1] = x [i1] ;
                break ;
            }
            i2 = z [l] ;
            i1stop = i2 ;
            l++ ;
            i2stop = z [l] ;
            while ( 1 )
            {
                if ( x [i1] <= x [i2] )  /* pop x [i1] */
                {
                    y [iy] = x [i1] ;
                    iy++ ;
                    i1++ ;
                    if ( i1 == i1stop )  /* dump x [i2:i2stop] to y [iy: ] */
                    {
                        yp = y+(iy-i2) ;
                        iy += i2stop - i2 ;
                        for (; i2 < i2stop; i2++) yp [i2] = x [i2] ;
                        break ;
                    }
                }
                else                     /* pop x [i2] */
                {
                    y [iy] = x [i2] ;
                    iy++ ;
                    i2++ ;
                    if ( i2 == i2stop )  /* dump x [i1:i1stop] to y [iy: ] */
                    {
                        yp = y+(iy-i1) ;
                        iy += i1stop - i1 ;
                        for (; i1 < i1stop; i1++) yp [i1] = x [i1] ;
                        break ;
                    }
                }
            }
        }
        z [ystart] = n ;
        if ( ystart == 1 )               /* done, sorted numbers in y */
        {
            for (i = 0; i < n; i++) x [i] = y [i] ;
            return ;
        }

        /* shuffle y into x, save break points */
        xstart = 0 ;                     /* new start */
        ix = 0 ;
        for (l = 0; l < ystart; )
        {
            z [xstart] = ix ;
            xstart++ ;
            i1 = z [l] ;
            l++ ;
            if ( l == ystart )           /* copy rest of y into x */
            {
                for (; i1 < n; i1++) x [i1] = y [i1] ;
                break ;
            }
            i2 = z [l] ;
            i1stop = i2 ;
            l++ ;
            i2stop = z [l] ;
            while ( 1 )
            {
                if ( y [i1] <= y [i2] )  /* pop y [i1] */
                {
                    x [ix] = y [i1] ;
                    ix++ ;
                    i1++ ;
                    if ( i1 == i1stop )  /* dump y [i2:i2stop] to x [ix: ] */
                    {
                        xp = x+(ix-i2) ;
                        ix += i2stop - i2 ;
                        for (; i2 < i2stop; i2++) xp [i2] = y [i2] ;
                        break ;
                    }
                }
                else                     /* pop y [i2] */
                {
                    x [ix] = y [i2] ;
                    ix++ ;
                    i2++ ;
                    if ( i2 == i2stop )    /* dump y [i1:i1stop] queue */
                    {
                        xp = x+(ix-i1) ;
                        ix += i1stop - i1 ;
                        for (; i1 < i1stop; i1++) xp [i1] = y [i1] ;
                        break ;
                    }
                }
            }
        }
        z [xstart] = n ;
    }
}

/* ==========================================================================
   === BLGpartialMinSort ====================================================
   ==========================================================================
   find indices of m smallest numbers in a list and arrange in a max heap
   ========================================================================== */
PRIVATE void BLGpartialMinSort
(
    int       *i, /* indices of m smallest numbers */
    BLGFLOAT  *x, /* numbers to sort */
    int        m, /* length of i */
    int        n  /* length of x */
)
{
    int j ;
    BLGFLOAT top ;
    if ( (m <= 0) || (n <= 0) ) return ;
    if ( m >= n )
    {
        for (j = 0; j < n; j++) i [j] = j ;
        BLGmaxheap_build (i-1, x, n) ;
        return ;
    }
    for (j = 0; j < m; j++) i [j] = j ;
    BLGmaxheap_build (i-1, x, m) ;
    top = x [i [0]] ;
    for (; j < n; j++)
    {
        if ( top > x [j] )
        {
            i [0] = j ;
            BLGmaxheapify (i-1, x, 1, m) ;
            top = x [i [0]] ;
        }
    }
}

/* ==========================================================================
   === BLGpartialMaxSort ====================================================
   ==========================================================================
   find indices of m largest numbers in a list
   ========================================================================== */
PRIVATE void BLGpartialMaxSort
(
    int      *i, /* indices of m largest numbers */
    BLGFLOAT *x, /* numbers to sort */
    int       m, /* length of i */
    int       n  /* length of x */
)
{
    int j ;
    BLGFLOAT bot ;
    if ( m >= n )
    {
        for (j = 0; j < n; j++) i [j] = j ;
        return ;
    }
    if ( (m <= 0) || (n <= 0) ) return ;
    for (j = 0; j < m; j++) i [j] = j ;
    BLGminheap_build (i-1, x, m) ;
    bot = x [i [0]] ;
    for (; j < n; j++)
    {
        if ( bot < x [j] )
        {
            i [0] = j ;
            BLGminheapify (i-1, x, 1, m) ;
            bot = x [i [0]] ;
        }
    }
}

/* ==========================================================================
   === BLGminHeapSort =======================================================
   ==========================================================================
   given a max heap, generate the list of indices corresponding to increasing
   order
   ========================================================================== */
PRIVATE void BLGminHeapSort
(
    int       *i, /* indices with increasing order based on x */
    int    *heap, /* min heap */
    BLGFLOAT  *x, /* numbers to sort */
    int        n  /* size of heap */
)
{
    int j, k, l ;
    l = n ;
    for (j = 0; j < l; j++)
    {
        k = heap [0] ;
        n = BLGmaxheap_delete (heap-1, x, n) ;
        i [n] = k ;
    }
}

/* ==========================================================================
   === minheap_build ========================================================
   ==========================================================================
   build a min heap in heap [1..nheap]
   heap [i] is an index of an element of array x
   ========================================================================== */

PRIVATE void BLGminheap_build
(
    int    *heap, /* on input, an unsorted set of indices */
    BLGFLOAT  *x, /* numbers to sort */
    int        n  /* number of elements to build into the heap */
)
{
    int p ;

    for (p = n/2 ; p >= 1 ; p--)
    {
        BLGminheapify (heap, x, p, n) ;
    }
}

/* ==========================================================================
   === minheap_delete =======================================================
   ==========================================================================
   delete the top element in a min heap
   ========================================================================== */
PRIVATE int BLGminheap_delete  /* return new size of heap */
(
    int    *heap, /* containing indices into x, 1..n on input */
    BLGFLOAT  *x, /* not modified */
    int        n  /* number of items in heap */
)
{
    if (n <= 1)
    {
        return (0) ;
    }

    /* move element from the end of the heap to the top */
    heap [1] = heap [n] ;
    n-- ;
    BLGminheapify (heap, x, 1, n) ;
    return (n) ;
}

/* ========================================================================== */
/* === minheap_add ========================================================== */
/* ========================================================================== */
/* ========================================================================== */

/* add a new leaf to a min heap */

PRIVATE int BLGminheap_add
(
    int     leaf, /* the new leaf */
    int    *heap, /* size n, containing indices into x */
    BLGFLOAT  *x, /* not modified */
    int        n  /* number of elements in heap not counting new one */
)
{
    int l, new, old ;
    BLGFLOAT xold, xnew ;

    n++ ;
    old = n ;
    heap [old] = leaf ;
    xold = x [leaf] ;
    while ( old > 1 )
    {
        new = old/2 ;
        l = heap [new] ;
        xnew = x [l] ;
        if ( xnew > xold )
        { /* swap new and old */
            heap [new] = leaf ;
            heap [old] = l ;
        }
        else return (n) ;
        old = new ;
    }
    return (n) ;
}

/* ========================================================================== */
/* === minheapify =========================================================== */
/* ========================================================================== */

/* heapify starting at node p.  On input, the heap at node p satisfies the */
/* heap property, except for heap [p] itself.  On output, the whole heap */
/* satisfies the heap property. */

PRIVATE void BLGminheapify
(
    int   *heap, /* size n, containing indices into x */
    BLGFLOAT *x, /* not modified */
    int       p, /* start at node p in the heap */
    int       n  /* heap [1 ... n] is in use */
)
{
    int left, right, e, hleft, hright ;
    BLGFLOAT xe, xleft, xright ;

    e = heap [p] ;
    xe = x [e] ;

    while ( 1 )
    {
        left = p * 2 ;
        right = left + 1 ;

        if ( right <= n )
        {
            hleft  = heap [left] ;
            hright = heap [right] ;
            xleft  = x [hleft] ;
            xright = x [hright] ;
            if (xleft < xright)
            {
                if (xe > xleft)
                {
                    heap [p] = hleft ;
                    p = left ;
                }
                else
                {
                    heap [p] = e ;
                    return ;
                }
            }
            else
            {
                if (xe > xright)
                {
                    heap [p] = hright ;
                    p = right ;
                }
                else
                {
                    heap [p] = e ;
                    return ;
                }
            }
        }
        else
        {
            if ( left <= n )
            {
                hleft = heap [left] ;
                xleft = x [hleft] ;
                if (xe > xleft)
                {
                    heap [p] = hleft ;
                    p = left ;
                }
            }
            heap [p] = e ;
            return ;
        }
    }
}

/* ========================================================================== */
/* === maxheap_build ======================================================== */
/* ========================================================================== */

/* build a max heap in heap [1..n] */

PRIVATE void BLGmaxheap_build
(
    int    *heap, /* on input, an unsorted set of elements */
    BLGFLOAT  *x,
    int        n  /* number of elements to build into the heap */
)
{
    int p ;

    for (p = n/2 ; p >= 1 ; p--)
    {
        BLGmaxheapify (heap, x, p, n) ;
    }
}

/* ========================================================================== */
/* === maxheap_delete ======================================================= */
/* ========================================================================== */

/* delete the top element in a max heap */

PRIVATE int BLGmaxheap_delete  /* return new size of heap */
(
    int    *heap, /* containing indices into x, 1..n on input */
    BLGFLOAT  *x, /* not modified */
    int        n  /* number of items in heap */
)
{
    if (n <= 1)
    {
        return (0) ;
    }

    /* move element from the end of the heap to the top */
    heap [1] = heap [n] ;
    n-- ;
    BLGmaxheapify (heap, x, 1, n) ;
    return (n) ;
}

/* ========================================================================== */
/* === maxheap_add ========================================================== */
/* ========================================================================== */

/* add a new leaf to a max heap */

PRIVATE int BLGmaxheap_add
(
    int     *heap, /* size n, containing indices into x */
    BLGFLOAT   *x, /* not modified */
    int      leaf, /* the new leaf */
    int         n  /* number of elements in heap not counting new one */
)
{
    int l, new, old ;
    BLGFLOAT xold, xnew ;

    n++ ;
    old = n ;
    heap [old] = leaf ;
    xold = x [leaf] ;
    while ( old > 1 )
    {
        new = old/2 ;
        l = heap [new] ;
        xnew = x [l] ;
        if ( xnew < xold )
        { /* swap new and old */
            heap [new] = leaf ;
            heap [old] = l ;
        }
        else return (n) ;
        old = new ;
    }
    return (n) ;
}

/* ========================================================================== */
/* === maxheapify =========================================================== */
/* ========================================================================== */

/* heapify starting at node p.  On input, the heap at node p satisfies the */
/* heap property, except for heap [p] itself.  On output, the whole heap */
/* satisfies the heap property. */

PRIVATE void BLGmaxheapify
(
    int   *heap, /* size n, containing indices into x */
    BLGFLOAT *x, /* not modified */
    int       p, /* start at node p in the heap */
    int       n  /* heap [1 ... n] is in use */
)
{
    int left, right, e, hleft, hright ;
    BLGFLOAT xe, xleft, xright ;

    e = heap [p] ;
    xe = x [e] ;

    while ( 1 )
    {
        left = p * 2 ;
        right = left + 1 ;

        if (right <= n)
        {
            hleft  = heap [left] ;
            hright = heap [right] ;
            xleft  = x [hleft] ;
            xright = x [hright] ;
            if (xleft > xright)
            {
                if (xe < xleft)
                {
                    heap [p] = hleft ;
                    p = left ;
                }
                else
                {
                    heap [p] = e ;
                    return ;
                }
            }
            else
            {
                if (xe < xright)
                {
                    heap [p] = hright ;
                    p = right ;
                }
                else
                {
                    heap [p] = e ;
                    return ;
                }
            }
        }
        else
        {
            if (left <= n)
            {
                hleft = heap [left] ;
                xleft = x [hleft] ;
                if (xe < xleft)
                {
                    heap [p] = hleft ;
                    p = left ;
                }
            }
            heap [p] = e ;
            return ;
        }
    }
}

/* ==========================================================================
   === BLGnapsack ===========================================================
   ==========================================================================
    Find x that minimizes ||x-y|| while satisfying lo <= x <= hi, a'x = b.
    The approach is to solve the dual problem obtained by introducing
    a multiplier lambda for the constraint a'x = b.  The dual function is

    L (lambda) = min { .5||x-y||^2 + lambda (a'x - b): lo <= x <= hi }

    The dual function is concave. It is continuously differentiable.
    If mu denotes the maximizer of the dual function, then the solution
    of the primal problem is

    x = proj (y - mu*a) ,

    where proj (z) is the projection of z onto the set { x : lo <= x <= hi }.
    Thus we have

       proj (z)_i = lo_i if z_i <= lo_i,
                    hi_i if z_i >= hi_i,
                    z_i otherwise  .

    Note that for any lambda, the minimizing x in the dual function is

       x (lambda) = proj (y - lambda*a).

    The slope of the dual function is

      L'(lambda) = a'x(lambda) - b

    The solution technique is to start with an initial guess lambda for
    mu and search for a zero of L'. We have the following cases:

  1. L'(lambda) >= 0: mu >= lambda. If L' = 0, then done. Otherwise,
                      increase lambda using napup until slope vanishes

  2. L'(lambda) <= 0: mu <= lambda. If L' = 0, then done. Otherwise,
                      decrease lambda using napdown until slope vanishes

    The free set is those i for which lo_i < x_i (lambda) < hi_i.
    The active set is those i for which x_i = lo_i or x_i = hi_i.
    If we have a guess for the optimal active and free indices,
    then we can obtain a good guess for the starting lambda by
    setting the slope of the dual function to zero and solving for lambda.
    The guess for the active and free set is specified using the array B
    which has values:

        +1 => x_i = hi_i
        -1 => x_i = lo_i
         0 => lo_i < x_i < hi_i

    If there is no guess, then B is NULL.  Note that B is an INPUT array,
    it is not computed within the routine.

    The total time taken by napsack is O (n + h log n), where n is the size of
    y, h is the number of times an element of x (lambda) moves off the boundary
    into the free set plus the number of times elements move from the
    free set to the opposite boundary.  A heap is used to hold the entries
    in the boundary and in the free set.  If the slope vanishes at either
    the starting lambda, then no heap is constructed, and the time is
    just O (n).
   ========================================================================== */

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
)
{

    int a_is_1, bnd_is_scalar, j, status ;
    BLGFLOAT Clo, Chi, xj, atx, a2sum, aj, lambda, slope ;

    lambda = *Lambda ;
    bnd_is_scalar = Com->Parm->bnd_is_scalar ; /* TRUE means scalar bounds */
    if ( bnd_is_scalar == BLGTRUE )
    {
        Chi = *hi ;
        Clo = *lo ;
    }
    a_is_1 = Com->Parm->a_is_1 ;   /* TRUE means the a vector is 1 identically*/

    /* ---------------------------------------------------------------------- */
    /* compute the initial slope */
    /* ---------------------------------------------------------------------- */

    atx = BLGZERO ;
    a2sum = BLGZERO ;
    if ( a_is_1 == BLGTRUE )
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
            for (j = 0; j < n; j++)
            {
                xj = x [j] - lambda ;
                if      ( xj >= Chi ) atx += Chi ;
                else if ( xj <= Clo ) atx += Clo ;
                else                { atx += x [j] ; a2sum++ ; }
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                xj = x [j] - lambda ;
                if      ( xj >= hi [j] ) atx += hi [j] ;
                else if ( xj <= lo [j] ) atx += lo [j] ;
                else                   { atx += x [j] ; a2sum++ ; }
            }
        }
    }
    else
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
            for (j = 0; j < n; j++)
            {
                aj = a [j] ;
                xj = x [j] - lambda*aj ;
                if      ( xj >= Chi ) atx += aj*Chi ;
                else if ( xj <= Clo ) atx += aj*Clo ;
                else                { atx += aj*x [j] ; a2sum += aj*aj ; }
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                aj = a [j] ;
                xj = x [j] - lambda*aj ;
                if      ( xj >= hi [j] ) atx += aj*hi [j] ;
                else if ( xj <= lo [j] ) atx += aj*lo [j] ;
                else                   { atx += aj*x [j] ; a2sum += aj*aj ; }
            }
        }
    }
    slope = atx - lambda*a2sum - b ; /* slope of the dual function */
    if      ( slope > BLGZERO ) /* case 1, increase lambda */
    {
        status = BLGnapup (x, &lambda, a, b, lo, hi, n, atx, a2sum, Com->w1,
                           Com->iw1-1, Com->iw2-1, Com);
    }
    else if ( slope < BLGZERO ) /* case 2, decrease lambda */
    {
        status = BLGnapdown (x, &lambda, a, b, lo, hi, n, atx, a2sum, Com->w1,
                             Com->iw1-1, Com->iw2-1, Com) ;
    }
    else status = 0 ;           /* slope = 0, exact starting guess */
    if ( status ) return (status) ;

    /* ---------------------------------------------------------------------- */
    /* replace x by x (lambda) */
    /* ---------------------------------------------------------------------- */
    if ( a_is_1 == BLGTRUE )
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
            for (j = 0; j < n; j++)
            {
                xj = x [j] - lambda ;
                if      ( xj < Clo ) xj = Clo ;
                else if ( xj > Chi ) xj = Chi ;
                x [j] = xj ;
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                xj = x [j] - lambda ;
                if      ( xj < lo [j] ) xj = lo [j] ;
                else if ( xj > hi [j] ) xj = hi [j] ;
                x [j] = xj ;
            }
        }
    }
    else
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
            for (j = 0; j < n; j++)
            {
                xj = x [j] - lambda*a [j] ;
                if      ( xj < Clo ) xj = Clo ;
                else if ( xj > Chi ) xj = Chi ;
                x [j] = xj ;
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                xj = x [j] - lambda*a [j] ;
                if      ( xj < lo [j] ) xj = lo [j] ;
                else if ( xj > hi [j] ) xj = hi [j] ;
                x [j] = xj ;
            }
        }
    }

    *Lambda = lambda ;
    return (status) ;
}

/* ========================================================================== */
/* === BLGnapup ============================================================= */
/* ========================================================================== */

/* Find x that minimizes ||x-y|| while satisfying the constraints
   lo <= x <= hi, a'x = b.  The algorithm is described in the napsack comments.
   It is assumed that the starting guess lambda for the dual multiplier is <=
   the correct multiplier. Hence, lambda will be increased.  The slope of
   the dual function is a'x - b. We increase lambda until a'x reaches b.
   If a_i > 0, then as lambda increases, x_i decreases. If a_i < 0, then
   as lambda increases, x_i increases. */

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

)
{
    int a_is_1, bnd_is_scalar, j, k, maxsteps, n_bound, n_free ;
    BLGFLOAT aj, Clo, Chi, lambda, minbound, minfree, new_break, slope, t, xj ;

    lambda = *Lambda ;
    a_is_1 = Com->Parm->a_is_1 ;   /* TRUE means the a vector is 1 identically*/
    bnd_is_scalar = Com->Parm->bnd_is_scalar ; /* TRUE means scalar bounds */
    minbound = BLGINF ;
    minfree  = BLGINF ;

    /* -------------------------------------------------------------- */
    /* construct the heaps */
    /* -------------------------------------------------------------- */

    n_bound = 0 ;
    n_free = 0 ;
    if ( a_is_1 == BLGTRUE )
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
            Clo = *lo ;
            Chi = *hi ;
            for (j = 0; j < n; j++)
            {
                xj = y [j] - lambda ;
                if ( xj >= Chi )
                {
                    n_bound++ ;
                    bound_heap [n_bound] = j ;
                    t = y [j] - Chi ;
                    if ( minbound > t ) minbound = t ;
                    breakpts [j] = t ;
                }
                else if ( xj > Clo )
                {
                    n_free++ ;
                    free_heap [n_free] = j ;
                    t = y [j] - Clo ;
                    if ( minfree > t ) minfree = t ;
                    breakpts [j] = t ;
                }
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                xj = y [j] - lambda ;
                if ( xj >= hi [j] )
                {
                    n_bound++ ;
                    bound_heap [n_bound] = j ;
                    t = y [j] - hi [j] ;
                    if ( minbound > t ) minbound = t ;
                    breakpts [j] = t ;
                }
                else if ( xj > lo [j] )
                {
                    n_free++ ;
                    free_heap [n_free] = j ;
                    t = y [j] - lo [j] ;
                    if ( minfree > t ) minfree = t ;
                    breakpts [j] = t ;
                }
            }
        }
    }
    else
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
            Clo = *lo ;
            Chi = *hi ;
            for (j = 0 ; j < n ; j++)
            {
                aj = a [j] ;
                xj = y [j] - aj*lambda ;
                if ( aj > BLGZERO )
                {
                    if ( xj >= Chi )
                    {
                        n_bound++ ;
                        bound_heap [n_bound] = j ;
                        t = (y [j]-Chi)/aj ;
                        if ( minbound > t ) minbound = t ;
                        breakpts [j] = t ;
                    }
                    else if ( xj > Clo )
                    {
                        n_free++ ;
                        free_heap [n_free] = j ;
                        t = (y [j]-Clo)/aj ;
                        if ( minfree > t ) minfree = t ;
                        breakpts [j] = t ;
                    }
                }
                else if ( aj < BLGZERO )
                {
                    if ( xj <= Clo )
                    {
                        n_bound++ ;
                        bound_heap [n_bound] = j ;
                        t = (y [j]-Clo)/aj ;
                        if ( minbound > t ) minbound = t ;
                        breakpts [j] = t ;
                    }
                    else if ( xj < Chi )
                    {
                        n_free++ ;
                        free_heap [n_free] = j ;
                        t = (y [j]-Chi)/aj ;
                        if ( minfree > t ) minfree = t ;
                        breakpts [j] = t ;
                    }
                }
            }
        }
        else
        {
            for (j = 0 ; j < n ; j++)
            {
                aj = a [j] ;
                xj = y [j] - aj*lambda ;
                if ( aj > BLGZERO )
                {
                    if ( xj >= hi [j] )
                    {
                        n_bound++ ;
                        bound_heap [n_bound] = j ;
                        t = (y [j]-hi [j])/aj ;
                        if ( minbound > t ) minbound = t ;
                        breakpts [j] = t ;
                    }
                    else if ( xj > lo [j] )
                    {
                        n_free++ ;
                        free_heap [n_free] = j ;
                        t = (y [j]-lo [j])/aj ;
                        if ( minfree > t ) minfree = t ;
                        breakpts [j] = t ;
                    }
                }
                else if ( aj < BLGZERO )
                {
                    if ( xj <= lo [j] )
                    {
                        n_bound++ ;
                        bound_heap [n_bound] = j ;
                        t = (y [j]-lo [j])/aj ;
                        if ( minbound > t ) minbound = t ;
                        breakpts [j] = t ;
                    }
                    else if ( xj < hi [j] )
                    {
                        n_free++ ;
                        free_heap [n_free] = j ;
                        t = (y [j]-hi [j])/aj ;
                        if ( minfree > t ) minfree = t ;
                        breakpts [j] = t ;
                    }
                }
            }
        }
    }

    maxsteps = 2*n + 1;
    for (k = 1 ; k <= maxsteps ; k++)
    {
        new_break = BLGMIN (minfree, minbound) ;
        slope = atx - new_break*a2sum - b ;
        if ( (slope <= BLGZERO) || (new_break == BLGINF) ) /* done */
        {
            if ( new_break == BLGINF ) /* problem is likely infeasible */
            {
                if ( fabs(slope) >= Com->anorm*Com->Parm->feas_tol ) return (1);
            }
            if ( a2sum != BLGZERO ) lambda = (atx-b)/a2sum ;
            *Lambda = lambda ;
            return (0) ;
        }
        lambda = new_break ;

        if ( k == 1 )
        {
            BLGminheap_build (free_heap, breakpts, n_free) ;
            BLGminheap_build (bound_heap, breakpts, n_bound) ;
        }

        /* -------------------------------------------------------------- */
        /* update the heaps */
        /* -------------------------------------------------------------- */

        if ( n_free > 0 )
        {
            if ( a_is_1 == BLGTRUE )
            {
                while (breakpts [j = free_heap [1]] <= lambda)
                {
                    a2sum-- ;
                    if ( bnd_is_scalar == BLGTRUE ) atx +=  Clo - y [j] ;
                    else                            atx +=  lo [j] - y [j] ;
                    n_free = BLGminheap_delete (free_heap, breakpts, n_free) ;
                    if ( n_free == 0 ) break ;
                }
            }
            else
            {
                while (breakpts [j = free_heap [1]] <= lambda)
                {
                    aj = a [j] ;
                    a2sum -= aj*aj ;
                    if ( aj > BLGZERO )
                    {
                        if (bnd_is_scalar == BLGTRUE) atx += aj*(Clo-y [j]) ;
                        else                          atx += aj*(lo [j]-y [j]) ;
                    }
                    else
                    {
                        if (bnd_is_scalar == BLGTRUE) atx += aj*(Chi-y [j]) ;
                        else                          atx += aj*(hi [j]-y [j]) ;
                    }
                    n_free = BLGminheap_delete (free_heap, breakpts, n_free) ;
                    if ( n_free == 0 )
                    {
                        a2sum = BLGZERO ;
                        break ;
                    }
                }
            }
        }

        if ( n_bound > 0 )
        {
            if ( a_is_1 == BLGTRUE )
            {
                while (breakpts [j = bound_heap [1]] <= lambda)
                {
                    n_bound  = BLGminheap_delete (bound_heap, breakpts,n_bound);
                    a2sum++ ;
                    if ( bnd_is_scalar == BLGTRUE )
                    {
                        atx += y [j] - Chi ;
                        breakpts [j] = y [j] - Clo ;
                    }
                    else
                    {
                        atx += y [j] - hi [j] ;
                        breakpts [j] = y [j] - lo [j] ;
                    }
                    n_free = BLGminheap_add (j, free_heap, breakpts, n_free) ;
                    if ( n_bound == 0 ) break ;
                }
            }
            else
            {
                while (breakpts [j = bound_heap [1]] <= lambda)
                {
                    n_bound  = BLGminheap_delete(bound_heap, breakpts, n_bound);
                    aj = a [j] ;
                    a2sum += aj*aj ;
                    if ( aj > BLGZERO )
                    {
                        if ( bnd_is_scalar == BLGTRUE )
                        {
                            atx += aj*(y [j] - Chi) ;
                            breakpts [j] = (y [j]-Clo)/aj ;
                        }
                        else
                        {
                            atx += aj*(y [j] - hi [j]) ;
                            breakpts [j] = (y [j]-lo [j])/aj ;
                        }
                    }
                    else
                    {
                        if ( bnd_is_scalar == BLGTRUE )
                        {
                            atx += aj*(y [j] - Clo) ;
                            breakpts [j] = (y [j]-Chi)/aj ;
                        }
                        else
                        {
                            atx += aj*(y [j] - lo [j]) ;
                            breakpts [j] = (y [j]-hi [j])/aj ;
                        }
                    }
                    n_free = BLGminheap_add (j, free_heap, breakpts, n_free) ;
                    if ( n_bound == 0 ) break ;
                }
            }
        }

        /*------------------------------------------------------------------- */
        /* get the smallest entry in each heap */
        /*------------------------------------------------------------------- */

        if ( n_free > 0 )  minfree = breakpts [free_heap [1]] ;
        else               minfree = BLGINF ;
        if ( n_bound > 0 ) minbound = breakpts [bound_heap [1]] ;
        else               minbound = BLGINF ;
    }
    return (1) ;
}

/* ========================================================================== */
/* === BLGnapdown =========================================================== */
/* ========================================================================== */

/* Find x that minimizes ||x-y|| while satisfying the constraints
   lo <= x <= hi, a'x = b.  The algorithm is described in the napsack comments.
   It is assumed that the starting guess lambda for the dual multiplier is >=
   the correct multiplier. Hence, lambda will be decreased.  The slope of
   the dual function is a'x - b. We decrease lambda until a'x reaches b.
   If a_i > 0, then as lambda decreases, x_i increases. If a_i < 0, then
   as lambda decreases, x_i decreases. */

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
)
{
    int a_is_1, bnd_is_scalar, j, k, maxsteps, n_bound, n_free ;
    BLGFLOAT aj, Clo, Chi, lambda, maxbound, maxfree, new_break, slope, t, xj ;

    lambda = *Lambda ;
    a_is_1 = Com->Parm->a_is_1 ;   /* TRUE means the a vector is 1 identically*/
    bnd_is_scalar = Com->Parm->bnd_is_scalar ; /* TRUE means scalar bounds */
    maxbound = -BLGINF ;
    maxfree  = -BLGINF ;

    /* -------------------------------------------------------------- */
    /* construct the heaps */
    /* -------------------------------------------------------------- */

    n_bound = 0 ;
    n_free = 0 ;
    if ( a_is_1 == BLGTRUE )
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
            Clo = *lo ;
            Chi = *hi ;
            for (j = 0; j < n; j++)
            {
                xj = y [j] - lambda ;
                if ( xj <= Clo )
                {
                    n_bound++ ;
                    bound_heap [n_bound] = j ;
                    t = y [j] - Clo ;
                    if ( maxbound < t ) maxbound = t ;
                    breakpts [j] = t ;
                }
                else if ( xj < Chi )
                {
                    n_free++ ;
                    free_heap [n_free] = j ;
                    t = y [j] - Chi ;
                    if ( maxfree < t ) maxfree = t ;
                    breakpts [j] = t ;
                }
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                xj = y [j] - lambda ;
                if ( xj <= lo [j] )
                {
                    n_bound++ ;
                    bound_heap [n_bound] = j ;
                    t = y [j] - lo [j] ;
                    if ( maxbound < t ) maxbound = t ;
                    breakpts [j] = t ;
                }
                else if ( xj < hi [j] )
                {
                    n_free++ ;
                    free_heap [n_free] = j ;
                    t = y [j] - hi [j] ;
                    if ( maxfree < t ) maxfree = t ;
                    breakpts [j] = t ;
                }
            }
        }
    }
    else
    {
        if ( bnd_is_scalar == BLGTRUE )
        {
            Clo = *lo ;
            Chi = *hi ;
            for (j = 0 ; j < n ; j++)
            {
                aj = a [j] ;
                xj = y [j] - aj*lambda ;
                if ( aj > BLGZERO )
                {
                    if ( xj <= Clo )
                    {
                        n_bound++ ;
                        bound_heap [n_bound] = j ;
                        t = (y [j]-Clo)/aj ;
                        if ( maxbound < t ) maxbound = t ;
                        breakpts [j] = t ;
                    }
                    else if ( xj < Chi )
                    {
                        n_free++ ;
                        free_heap [n_free] = j ;
                        t = (y [j]-Chi)/aj ;
                        if ( maxfree < t ) maxfree = t ;
                        breakpts [j] = t ;
                    }
                }
                else if ( aj < BLGZERO )
                {
                    if ( xj >= Chi )
                    {
                        n_bound++ ;
                        bound_heap [n_bound] = j ;
                        t = (y [j]-Chi)/aj ;
                        if ( maxbound < t ) maxbound = t ;
                        breakpts [j] = t ;
                    }
                    else if ( xj > Clo )
                    {
                        n_free++ ;
                        free_heap [n_free] = j ;
                        t = (y [j]-Clo)/aj ;
                        if ( maxfree < t ) maxfree = t ;
                        breakpts [j] = t ;
                    }
                }
            }
        }
        else
        {
            for (j = 0 ; j < n ; j++)
            {
                aj = a [j] ;
                if ( aj > BLGZERO )
                {
                    xj = y [j] - aj*lambda ;
                    if ( xj <= lo [j] )
                    {
                        n_bound++ ;
                        bound_heap [n_bound] = j ;
                        t = (y [j]-lo [j])/aj ;
                        if ( maxbound < t ) maxbound = t ;
                        breakpts [j] = t ;
                    }
                    else if ( xj < hi [j] )
                    {
                        n_free++ ;
                        free_heap [n_free] = j ;
                        t = (y [j]-hi [j])/aj ;
                        if ( maxfree < t ) maxfree = t ;
                        breakpts [j] = t ;
                    }
                }
                else if ( aj < BLGZERO )
                {
                    xj = y [j] - aj*lambda ;
                    if ( xj >= hi [j] )
                    {
                        n_bound++ ;
                        bound_heap [n_bound] = j ;
                        t = (y [j]-hi [j])/aj ;
                        if ( maxbound < t ) maxbound = t ;
                        breakpts [j] = t ;
                    }
                    else if ( xj > lo [j] )
                    {
                        n_free++ ;
                        free_heap [n_free] = j ;
                        t = (y [j]-lo [j])/aj ;
                        if ( maxfree < t ) maxfree = t ;
                        breakpts [j] = t ;
                    }
                }
            }
        }
    }

    maxsteps = 2*n + 1;
    for (k = 1; k <= maxsteps; k++)
    {
        new_break = BLGMAX (maxfree, maxbound) ;
        slope = atx - new_break*a2sum - b ;
        if ( (slope >= BLGZERO) || (new_break == -BLGINF) ) /* done */
        {
            if ( new_break == BLGINF ) /* problem is likely infeasible */
            {
                if ( fabs(slope) >= Com->anorm*Com->Parm->feas_tol ) return (1);
            }
            if ( a2sum != BLGZERO ) lambda = (atx-b)/a2sum ;
            *Lambda = lambda ;
            return (0) ;
        }
        lambda = new_break ;

        if ( k == 1 )
        {
            BLGmaxheap_build (free_heap, breakpts, n_free) ;
            BLGmaxheap_build (bound_heap, breakpts, n_bound) ;
        }

        /* -------------------------------------------------------------- */
        /* update the heaps */
        /* -------------------------------------------------------------- */

        if ( n_free > 0 )
        {
            if ( a_is_1 == BLGTRUE )
            {
                while (breakpts [j = free_heap [1]] >= lambda)
                {
                    a2sum-- ;
                    if ( bnd_is_scalar == BLGTRUE ) atx +=  Chi - y [j] ;
                    else                            atx +=  hi [j] - y [j] ;
                    n_free = BLGmaxheap_delete (free_heap, breakpts, n_free) ;
                    if ( n_free == 0 ) break ;
                }
            }
            else
            {
                while (breakpts [j = free_heap [1]] >= lambda)
                {
                    aj = a [j] ;
                    a2sum -= aj*aj ;
                    if ( aj > BLGZERO )
                    {
                        if (bnd_is_scalar == BLGTRUE) atx += aj*(Chi-y [j]) ;
                        else                          atx += aj*(hi [j]-y [j]) ;
                    }
                    else
                    {
                        if (bnd_is_scalar == BLGTRUE) atx += aj*(Clo-y [j]) ;
                        else                            atx += aj*(lo [j]-y [j]) ;
                    }
                    n_free = BLGmaxheap_delete (free_heap, breakpts, n_free) ;
                    if ( n_free == 0 )
                    {
                        a2sum = BLGZERO ;
                        break ;
                    }
                }
            }
        }

        if ( n_bound > 0 )
        {
            if ( a_is_1 == BLGTRUE )
            {
                while (breakpts [j = bound_heap [1]] >= lambda)
                {
                    n_bound  = BLGmaxheap_delete (bound_heap, breakpts,n_bound);
                    a2sum++ ;
                    if ( bnd_is_scalar == BLGTRUE )
                    {
                        atx += y [j] - Clo ;
                        breakpts [j] = y [j] - Chi ;
                    }
                    else
                    {
                        atx += y [j] - lo [j] ;
                        breakpts [j] = y [j] - hi [j] ;
                    }
                    n_free = BLGmaxheap_add (free_heap, breakpts, j, n_free) ;
                    if ( n_bound == 0 ) break ;
                }
            }
            else
            {
                while (breakpts [j = bound_heap [1]] >= lambda)
                {
                    n_bound  = BLGmaxheap_delete(bound_heap, breakpts, n_bound);
                    aj = a [j] ;
                    a2sum += aj*aj ;
                    if ( aj > BLGZERO )
                    {
                        if ( bnd_is_scalar == BLGTRUE )
                        {
                            atx += aj*(y [j] - Clo) ;
                            breakpts [j] = (y [j]-Chi)/aj ;
                        }
                        else
                        {
                            atx += aj*(y [j] - lo [j]) ;
                            breakpts [j] = (y [j]-hi [j])/aj ;
                        }
                    }
                    else
                    {
                        if ( bnd_is_scalar == BLGTRUE )
                        {
                            atx += aj*(y [j] - Chi) ;
                            breakpts [j] = (y [j]-Clo)/aj ;
                        }
                        else
                        {
                            atx += aj*(y [j] - hi [j]) ;
                            breakpts [j] = (y [j]-lo [j])/aj ;
                        }
                    }
                    n_free = BLGmaxheap_add (free_heap, breakpts, j, n_free) ;
                    if ( n_bound == 0 ) break ;
                }
            }
        }

        /*------------------------------------------------------------------- */
        /* get the biggest entry in each heap */
        /*------------------------------------------------------------------- */

        if ( n_free > 0 )  maxfree = breakpts [free_heap [1]] ;
        else               maxfree = -BLGINF ;
        if ( n_bound > 0 ) maxbound = breakpts [bound_heap [1]] ;
        else               maxbound = -BLGINF ;
    }
    return (1) ;
}
/* ========================================================================= */
/* === lp ==== min c'x subject to 1'x = b, 0 <= x <= C ===================== */
/* ========================================================================= */

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
)
{
    int j, k, l, m, nn, np, *indexP ;
    BLGFLOAT slope, t ;

    /* Find the negative and positive components of c_i + mu.
       Store the indices of smallest components of c in bottom of
       array index, store indices of largest components of c in top
       of array index */
    nn = 0 ;
    np = n ;
    for (j = 0; j < n; j++)
    {
        if ( c [j] + mu < BLGZERO )
        {
            x [j] = C ;       /* tentatively set x_j = C */
            index [nn] = j ;
            nn++ ; /* count number of Cs in x */
        }
        else
        {
            np-- ;
            x [j] = BLGZERO ; /* tentatively set x_j = 0 */
            index [np] = j ;
        }
    }

    /* dual function L(mu) = min { c'x + mu(1'x-b) : 0 <= x <= C }
       compute the slope of the dual function */
    slope = (C*(double) nn) - b ;

    /* If slope is positive, then mu is increased and more components
       of x go to 0. If slope is negative, then mu is decreased and
       more components of x go to C.  As mu crosses break points,
       the slope changes by C.  Count the number of break points that
       must be crossed to make the slope of the dual function vanish. */
    m = ceil (fabs(slope)/C) ;

    if ( slope > BLGZERO ) /* increase mu, use nn indices, Cs switch to  0 */
    {
        /* extract components of c for sorting */
        for (j = 0; j < nn; j++) y [j] = c [index [j]] ;

        /* extract the m indices in lower part of index array corresponding to
           the largest components of cost vector c, the values of the
           solution x are all changed from C to 0 except for the index
           associated with the smallest c_i. The corresponding x_i is set
           to a value between 0 and C so that the linear constraint is
           satisfied */

        BLGpartialMaxSort (i, y, m, nn) ;
        t = BLGINF ;
        for (j = 0; j < m; j++)
        {
            k = index [i [j]] ;
            x [k] = 0 ;
            if ( c [k] < t )
            {
                t = c [k] ;
                l = k ;
            }
        }
        x [l] = b - (nn-m)*C ;
        /* c_l + mu = 0, return mu = -c_l */
        return (-c [l]) ;
    }
    else if ( slope < BLGZERO ) /* decrease mu, use np indices, 0s switch to C*/
    {
        /* extract remaining components of c for sorting */
        for (j = np; j < n; j++) y [j] = c [index [j]] ;

        /* extract the m indices in upper part of index array corresponding
           to the smallest components of cost vector c, the values of the
           solution x are all changed from 0 to C except for the index
           associated with the largest c_i. The corresponding x_i is set
           to a value between 0 and C so that the linear constraint is
           satisfied */

        BLGpartialMinSort (i, y+np, m, n-np) ;
        t = -BLGINF ;
        indexP = index+np ;
        for (j = 0; j < m; j++)
        {
            k = indexP [i [j]] ;
            x [k] = C ;
            if ( c [k] > t )
            {
                t = c [k] ;
                l = k ;
            }
        }
        x [l] = b - (nn+m-1)*C ;
        /* c_l + mu = 0, return mu = -c_l */
        return (-c [l]) ;
    }
    else return (mu) ; /* slope = 0 implies the guess for mu was correct */
}

/* =========================================================================
   = BLGlp = min c'x subject to a'x = b, max{lo, -x_bnd} <= x <= min{hi, x_bnd}
   ========================================================================= */

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
)
{
    int a_is_1, bnd_is_scalar, nbrk, j, k, *Br_order, *Br_index, *iw ;
    BLGFLOAT aj, Clo, Chi, loj, hij, fd, s, t, *Br_value ;

    bnd_is_scalar = Com->Parm->bnd_is_scalar ; /* TRUE means scalar bounds */
    a_is_1 = Com->Parm->a_is_1 ; /* TRUE means the a vector is 1 */
    /* in initial LP, compute bound for x based on napsack problem
       ||x||_inf <= ||x-x0||_inf + ||x0||_inf
                 <= ||x-x0||_2   + ||x0||_inf
                 <= 2||g||_2/lambda_k + ||x0||_inf
       x0 = any feasible point (current iterate) */
    fd = BLGZERO ;
    if ( bnd_is_scalar == BLGTRUE )
    {
        Chi = *hi ;
        Clo = *lo ;
        Clo = BLGMAX (-x_bnd, Clo) ;
        Chi = BLGMIN ( x_bnd, Chi) ;
        /* if a = 1, then there is a faster routine to find solution */
        if ( a_is_1 == BLGTRUE )
        {
            b -= Clo*n ;
            mu = BLGlp1 (x, mu, c, Chi-Clo, b, n, Com->iw1, Com->iw2, Com->w1) ;
            if ( Clo != BLGZERO )
            {
                for (j = 0; j < n; j++) x [j] += Clo ;
            }
            return (mu) ;
        }
        mu = -mu ; /* sign change in dual function from rest of program */
        for (j = 0; j < n; j++)
        {
            if ( c [j] > mu*a [j] )
            {
                fd += Clo*a [j] ;
                x [j] = Clo ;
            }
            else
            {
                fd += Chi*a [j] ;
                x [j] = Chi ;
            }
        }
    }
    else  /* bound not scalar */
    {
        mu = -mu ; /* sign change in dual function from rest of program */
        for (j = 0; j < n; j++)
        {
            loj = BLGMAX (-x_bnd, lo [j]) ;
            hij = BLGMIN ( x_bnd, hi [j]) ;
            if ( a_is_1 == BLGTRUE )
            {
                if ( c [j] > mu )
                {
                    x [j] = loj ;
                    fd += loj ;
                }
                else
                {
                    x [j] = hij ;
                    fd += hij ;
                }
            }
            else
            {
                if ( c [j] > mu*a [j] )
                {
                    x [j] = loj ;
                    fd += loj*a [j] ;
                }
                else
                {
                    x [j] = hij ;
                    fd += hij*a [j] ;
                }
            }
        }
    }

    fd = b - fd ;

    /* dual function L(mu) = min { c'x + mu(b-a'x) : lo <= x <= hi }
       fd = slope of the dual function. Note change in sign for the
       term b-a'x from the rest of program */

    Br_value = Com->w1 ;
    Br_index = Com->iw1 ;
    Br_order = Com->iw2 ;
    iw = Com->iw3 ;        /* work array */
    nbrk = 0 ;
    if ( fd > BLGZERO ) /*increase mu (find all break points > current mu) */
    {
        if ( a_is_1 == BLGTRUE )
        {
            for (j = 0; j < n; j++)
            {
                if ( c [j] > mu )
                {
                    Br_value [nbrk] = c [j] ;
                    Br_index [nbrk] = j ;
                    nbrk++ ;
                }
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                t = a [j] ;
                s = c [j] ;
                if ( ((s >  mu*t) && (t > BLGZERO)) ||
                     ((s <= mu*t) && (t < BLGZERO)) )
                {
                    Br_value [nbrk] = s/t ;
                    Br_index [nbrk] = j ;
                    nbrk++ ;
                }
            }
        }
        BLGindexMinSort (Br_order, iw, Br_value, nbrk) ;
        if ( a_is_1 == BLGTRUE )
        {
            for (k = 0; k < nbrk; k++)
            {
                j = Br_index [Br_order [k]] ;
                loj = BLGMAX (-x_bnd, lo [j]) ;
                hij = BLGMIN ( x_bnd, hi [j]) ;
                fd += loj - hij ;
                x [j] = hij ;
                if ( fd <= BLGZERO ) break ;
            }
            x [j] += fd ;
        }
        else
        {
            if ( bnd_is_scalar == BLGTRUE )
            {
                for (k = 0; k < nbrk; k++)
                {
                    j = Br_index [Br_order [k]] ;
                    aj = a [j] ;
                    if ( aj > BLGZERO )
                    {
                        fd += (Clo - Chi)*aj ;
                        x [j] = Chi ;
                    }
                    else
                    {
                        fd += (Chi - Clo)*aj ;
                        x [j] = Clo ;
                    }
                    if ( fd <= BLGZERO ) break ;
                }
            }
            else
            {
                for (k = 0; k < nbrk; k++)
                {
                    j = Br_index [Br_order [k]] ;
                    aj = a [j] ;
                    loj = BLGMAX (-x_bnd, lo [j]) ;
                    hij = BLGMIN ( x_bnd, hi [j]) ;
                    if ( aj > BLGZERO )
                    {
                        fd += (loj - hij)*aj ;
                        x [j] = hij ;
                    }
                    else
                    {
                        fd += (hij - loj)*aj ;
                        x [j] = loj ;
                    }
                    if ( fd <= BLGZERO ) break ;
                }
            }
            x [j] += fd/aj ;
        }
    }
    else if ( fd < BLGZERO ) /* decrease mu (find break points <= current mu) */
    {
        if ( a_is_1 == BLGTRUE )
        {
            for (j = 0; j < n; j++)
            {
                if ( c [j] <= mu )
                {
                    Br_value [nbrk] = c [j] ;
                    Br_index [nbrk] = j ;
                    nbrk++ ;
                }
            }
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                t = a [j] ;
                s = c [j] ;
                if ( ((s >  mu*t) && (t < BLGZERO)) ||
                     ((s <= mu*t) && (t > BLGZERO)) )
                {
                    Br_value [nbrk] = s/t ;
                    Br_index [nbrk] = j ;
                    nbrk++ ;
                }
            }
        }
        BLGindexMaxSort (Br_order, iw, Br_value, nbrk) ;
        if ( a_is_1 == BLGTRUE )
        {
            for (k = 0; k < nbrk; k++)
            {
                j = Br_index [Br_order [k]] ;
                loj = BLGMAX (-x_bnd, lo [j]) ;
                hij = BLGMIN ( x_bnd, hi [j]) ;
                fd += hij - loj ;
                x [j] = loj ;
                if ( fd >= BLGZERO ) break ;
            }
            x [j] += fd ;
        }
        else
        {
            if ( bnd_is_scalar == BLGTRUE )
            {
                for (k = 0; k < nbrk; k++)
                {
                    j = Br_index [Br_order [k]] ;
                    aj = a [j] ;
                    if ( aj > BLGZERO )
                    {
                        fd += (Chi - Clo)*aj ;
                        x [j] = Clo ;
                    }
                    else
                    {
                        fd += (Clo - Chi)*aj ;
                        x [j] = Chi ;
                    }
                    if ( fd >= BLGZERO ) break ;
                }
            }
            else
            {
                for (k = 0; k < nbrk; k++)
                {
                    j = Br_index [Br_order [k]] ;
                    aj = a [j] ;
                    loj = BLGMAX (-x_bnd, lo [j]) ;
                    hij = BLGMIN ( x_bnd, hi [j]) ;
                    if ( aj > BLGZERO )
                    {
                        fd += (hij - loj)*aj ;
                        x [j] = loj ;
                    }
                    else
                    {
                        fd += (loj - hij)*aj ;
                        x [j] = hij ;
                    }
                    if ( fd >= BLGZERO ) break ;
                }
            }
            x [j] += fd/aj ;
        }
    }
    else return (-mu) ; /* starting guess for mu was correct */

    mu = Br_value [Br_order [k]] ;
    /* return -mu since the dual function above has the opposite sign of
       what is used elsewhere */
    return (-mu) ;
}

/* ==========================================================================
   === BLGdefault ===========================================================
   ==========================================================================
       Set default parameter values.
   ========================================================================== */
void BLGdefault
(
    BLGparm *Parm /* pointer to parameter structure */
)
{
/* =================== BLG LIMITS =========================================== */
    /* maximum number iterations in BLG algorithm */
    Parm->maxits = BLGIINT ;

/* =================== PRINTING ============================================= */
    /* T => print parameters values
       F => do not display parameter values */
    Parm->PrintParms = BLGTRUE ;

    /* T => print final statistics
       F => no printout of statistics */
    Parm->PrintFinal = BLGTRUE ;

    /* Level 0 = no printing, ... , Level 2 = maximum printing */
    Parm->PrintLevel = 0 ;

/* =================== PROBLEM DESCRIPTION ================================== */
    /* TRUE means that the objective function is quadratic */
    Parm->QuadCost = BLGFALSE ;

    /* TRUE means that when the objective function is quadratic, the
       code should compute the starting gradient Ax + c and the
       starting objective function value .5x'Ax + c'x */
    Parm->QuadLazy = BLGTRUE ;

    /* TRUE means that the a vector is identically 1 */
    Parm->a_is_1 = BLGFALSE ;

    /* TRUE means that the lower and upper bounds are scalar */
    Parm->bnd_is_scalar = BLGFALSE ;

/* =================== FEASIBILITY TESTING ================================== */
    /* problem infeasible if |b-a'x| >= feas_tol*||a||_inf */
    Parm->feas_tol = 1.e-4 ;

    /* check_feas = TRUE means to check the problem for feasibility and
       project the starting guess on the feasible set */
    Parm->check_feas = BLGTRUE ;

/* =================== BB PARAMETER ========================================= */
    /* lower bound for lambda_k */
    Parm->lambda0 = 1.e-30 ;

    /* require lambda_k >= ||g||lambda0Factor */
    Parm->lambda0Factor = 1.e-6 ;

/* =================== BB CYCLE LENGTH ====================================== */
    /* nominal CBB cycle length */
    Parm->nm = 4 ;

    /* absolute maximum number of iterations in a CBB cycle */
    Parm->MaxCycle = 6 ;

/* =================== LINE SEARCH ========================================== */
    /* T => use approximate nonmonotone Armijo line search
       F => use ordinary nonmonotone Armijo line search, switch to
            approximate Armijo when |f_r-f| < AArmijoFac*|min (f_r, f_{max})|
       Set AArmijoFac = 0 and AArmijo = BLGFALSE to ensure ordinary Armijo */
    Parm->AArmijo = BLGFALSE ;

    Parm->AArmijoFac = 1.e-8 ;

    /* T => estimated error in function value = eps*|min (f_r, f_{max}) |
       F => estimated error in function value = eps */
    Parm->PertRule = BLGFALSE ;
    /* fmax = max (f_{k-i}, i = 0, 1, ..., min (k, M-1) ) */

    Parm->M = 8 ;

    /* terminate Armijo line search when
       f(st) <= fR + st * delta * f'(0) where fR = mfax or fr */
    Parm->Armijo_delta = 1.0e-4 ;

    /* backtracking decay factor in the Armijo line search */
    Parm->Armijo_decay = 5.e-1 ;

    /* maximum number of times the Armijo line search will perform
       backtracking steps */
    Parm->max_backsteps = (int) 50 ;

    /* required decay factor (safe-guard factor) for width of bracketing
       interval when computing mu */
    Parm->gamma = .66 ;
    Parm->eps = 1.e-6 ;

    /* search for non nan function value by shrinking search interval
       at most nshrink times */
    Parm->nshrink = (int) 50 ;

    /* factor by which interval shrinks when searching for non nan value */
    Parm->nan_fac = 2.e-1 ;

     /* update fr if fmin was not improved after L iterations */
    Parm->L = 3 ;

    /* update fr if initial stepsize was accepted in previous P iterations */
    Parm->P = 40 ;

    /* update reference value fr if (fr-fmin)/(fc-fmin) > gamma1 */
    Parm->gamma1 = (BLGFLOAT) Parm->M / (BLGFLOAT) Parm->L ;

    /* update fr if (fr-f)/(fmax-f) > gamma2, np > P, and fmax > f */
    Parm->gamma2 = (BLGFLOAT) Parm->P / (BLGFLOAT) Parm->M ;

    /* use quadratic interpolation to compute Armijo step if it
       lies in the interval [armijo0*alpha, armijo1*alpha] */
    Parm->armijo0 = 1.e-10 ;
    Parm->armijo1 = 9.e-1 ;

/* =================== SUBSPACE SELECTION =================================== */
    /* TRUE means that the algorithm works in a subspace */
    Parm->Subspace = BLGFALSE ;
    Parm->nsub = -1 ;  /* subspace size */
    Parm->fsub = .01 ; /* if nsub <= 0, subspace size = n*fsub */

    Parm->d_pert = 1.e-10 ; /* include d_i in subspace if >= d_pert*grad_tol */

    /* MaxKept is an upper bound on the fraction of indices from the
       previous subspace that are retained in the current subspace */
    Parm->MaxKept = .8 ;

    /* When deciding which indices to put in the subspace, we attach greater
       weight to the indices associated with the previous subspace.
       This weight factor is denoted PreAmp. */
    Parm->PreAmp = 1. ;

    /* If subspace error > MAX (err1, err2) repeat subspace.
       Default is to not repeat subspace (set serr1 and serr2 huge).
       If there is an advantage in repeating a subspace,
       possible parameter value might be serr1 = .25 and serr2 = .10 */

    Parm->serr1 = 1.e20 ;  /* err1 = serr1 times error tolerance */
    Parm->serr2 = 1.e20 ;  /* err2 = serr2 times KKT error */

    /* If lo_j + bnd_pert < x_j < hi_j + bnd_pert, then treat x_j as free
       in the buildSubspace routine */
    Parm->bnd_pert = 0.e-3 ;

/* =================== SEARCH DIRECTION ===================================== */
    /* If FW = TRUE means use Frank-Wolfe direction (has precedence
       over GP) */
    Parm->FW = BLGFALSE ;

    /* If GP = TRUE means use gradient projection direction, otherwise use
       affine scaling direction */
    Parm->GP = BLGTRUE ;

    Parm->lpfactor = 1.e6 ; /* parameter in decision rule of QP versus LP */

    /* set components of the search direction d to zero whenever the
       associated components of x are nearly active (within the cutoffs) */
    Parm->lo_cut = 0.e-12 ; /* g_j > 0 and x_j - lo_j < lo_cut => d_j = 0 */
    Parm->hi_cut = 0.e-12 ; /* g_j < 0 and hi_j - x_j < hi_cut => d_j = 0 */

/* =================== AFFINE SCALING CONTROL =============================== */
    /* absolute error in mu computation */
    Parm->epsmu0 = 1.e-11 ;

    /* relative error in mu computation */
    Parm->epsmu1 = 1.e-10 ;

    /* maximum number of iterations when computing mu */
    Parm->max_mu_its = 100 ;
}

/* ==========================================================================
   === BLGprint_parms =======================================================
   ==========================================================================
   print values in BLGparm structure
   ========================================================================== */
void BLGprint_parms
(
    BLGparm     *Parm, /* SSMparm structure to be printed */
    BLGFLOAT grad_tol, /* gradient tolerance */
    int             n
)
{
    int nsub ;
    printf ("\n=================== BLG LIMITS ===========================\n") ;
    printf ("Stopping criterion ............................. grad_tol: %e\n",
             grad_tol) ;
    printf ("Maximum number of iterations ..................... maxits: %i\n",
             Parm->maxits) ;

    printf ("\n=================== PRINTING =============================\n") ;
    printf ("Print final cost and statistics .............. PrintFinal: ") ;
    if ( Parm->PrintFinal == BLGTRUE ) printf ("TRUE\n") ;
    else                               printf ("FALSE\n") ;
    printf ("Print parameter structure .................... PrintParms: ") ;
    if ( Parm->PrintParms == BLGTRUE ) printf ("TRUE\n") ;
    else                               printf ("FALSE\n") ;
    printf ("Print level (0 = none, 2 = maximum) .......... PrintLevel: %i\n",
             Parm->PrintLevel) ;

    printf ("\n=================== PROBLEM DESCRIPTION ==================\n") ;
    printf ("The objective function is quadratic ............ QuadCost: ") ;
    if ( Parm->QuadCost == BLGTRUE ) printf ("TRUE\n") ;
    else                             printf ("FALSE\n") ;
    if ( Parm->QuadCost == BLGTRUE ) 
    {
        printf ("    User must provide routine to compute A*d\n") ;
        if ( Parm->QuadLazy == BLGTRUE )
        {
            printf ("    The g0 argument of BLG contains the "
                    "linear cost term c.\n") ;
            printf ("    The code should compute the starting gradient "
                    "Ax + c\nand the starting cost value .5x'Ax + c'x\n") ;
        }
        else
        {
            printf ("    The g0 argument of BLG contains the "
                    "start objective function gradient.\n") ;
            printf ("    The code returns the difference between the final "
                    "(optimal) objective function value\n") ;
            printf ("    and the starting value\n") ;
        }
    }
    printf ("The a vector in the linear constraint is 1 ....... a_is_1: ") ;
    if ( Parm->a_is_1 == BLGTRUE ) printf ("TRUE\n") ;
    else                           printf ("FALSE\n") ;
    printf ("The bounds are scalar ..................... bnd_is_scalar: ") ;
    if ( Parm->bnd_is_scalar == BLGTRUE ) printf ("TRUE\n") ;
    else                                  printf ("FALSE\n") ;

    printf ("\n=================== FEASIBILITY TESTING ==================\n") ;
    printf ("Check the problem for feasibility ............ check_feas: ") ;
    if ( Parm->check_feas == BLGTRUE ) printf ("TRUE\n") ;
    else                               printf ("FALSE\n") ;
    if ( Parm->check_feas == BLGTRUE )
        printf ("Feasibility tolerance .......................... feas_tol: "
             "%e\n", Parm->feas_tol) ;

    printf ("\n=================== BB PARAMETER =========================\n") ;
    printf ("Lower bound for lambda_k ........................ lambda0: %e\n",
             Parm->lambda0) ;
    printf ("Require lambda_k >= ||g||*lambda0Factor ... lambda0Factor: %e\n",
             Parm->lambda0Factor) ;

    printf ("\n=================== BB CYCLE LENGTH ======================\n") ;
    printf ("Nominal CBB cycle length ............................. nm: %i\n",
             Parm->nm) ;
    printf ("Absolute max number iteration in CBB cycle ..... MaxCycle: %i\n",
             Parm->MaxCycle) ;

    printf ("\n=================== LINE SEARCH ==========================\n") ;
    printf ("Approximate nonmonotone Armijo line search ...... AArmijo: ") ;
    if ( Parm->AArmijo == BLGTRUE ) printf ("TRUE\n") ;
    else                            printf ("FALSE\n") ;
    if ( Parm->AArmijo == BLGTRUE )
    {
        if ( Parm->PertRule == BLGTRUE )
            printf ("    Error estimate for function value is eps*|fcomp|\n") ;
        else
            printf ("    Error estimate for function value is eps\n") ;
    }
    else
    {
        if ( Parm->AArmijoFac > BLGZERO )
        {
            printf (" ... switch to approx nonmonotone Armijo line search\n") ;
            printf ("     cost change factor for transition ....... "
                    "AArmijoFac: %e\n", Parm->AArmijoFac) ;
        }
    }
    printf ("Fmax = max (f_{k-i}, i = 0, 1, ..., min (k, M-1) )..... M: %i\n",
             Parm->M) ;
    printf ("Perturbation parameter for function value............ eps: %e\n",
             Parm->eps) ;
    printf ("Max tries to find non NAN function value ........ nshrink: %i\n",
             Parm->nshrink) ;
    printf ("Interval decay factor in NAN search ..............nan_fac: %e\n",
             Parm->nan_fac) ;
    printf ("Update fr if fmin not improved after L iterations ..... L: %i\n",
             Parm->L) ;
    printf ("Update fr if P previous initial stepsizes accepted .... P: %i\n",
             Parm->P) ;
    printf ("Armijo line search parameter ............... Armijo_delta: %e\n",
             Parm->Armijo_delta) ;
    printf ("Armijo decay factor .........................armijo_decay: %e\n",
             Parm->Armijo_decay) ;
    printf ("Max number of Armijo backtracking steps ... max_backsteps: %i\n",
             Parm->max_backsteps) ;
    printf ("Criterion for updating reference value fr......... gamma1: %e\n",
             Parm->gamma1) ;
    printf ("Criterion for updating reference value fr......... gamma2: %e\n",
             Parm->gamma2) ;
    printf ("Criterion for Q interpolation,  cbb line search,  armijo0: %e\n",
             Parm->armijo0) ;
    printf ("Criterion for Q interpolation,  cbb line search,  armijo1: %e\n",
             Parm->armijo1) ;

    printf ("\n=================== SUBSPACE SELECTION ===================\n") ;
    printf ("The algorithm works in a subspace .............. Subspace: ") ;
    if ( Parm->Subspace == BLGTRUE ) printf ("TRUE\n") ;
    else                             printf ("FALSE\n") ;
    if ( Parm->Subspace == BLGTRUE )
    {
        nsub = Parm->nsub ;
        if ( nsub <= 0 )
        {
            nsub = n*Parm->fsub ;
            nsub = BLGMAX (2, nsub) ;
        }
        printf ("The dimension of the subspace is up to ..................: "
                "%i\n", nsub) ;
        printf ("Include d_i in subspace if >= d_pert*grad_tol .... d_pert: "
                "%e\n", Parm->d_pert) ;
        printf ("Amplification factor for value of prior indices .. PreAmp: "
                "%e\n", Parm->PreAmp) ;
        printf ("Maximum fraction of prior indices in new space .. MaxKept: "
                "%e\n", Parm->MaxKept) ;
        printf ("Repeat subspace if subspace error > MAX (err1, err2)\n") ;
        printf ("    where err1 = serr1*grad_tol and err2 = serr2*KKT error\n");
        printf ("................................................... serr1: "
                "%e\n", Parm->serr1) ;
        printf ("................................................... serr2: "
                "%e\n", Parm->serr2) ;
        printf ("Free variables satisfy the following relation:\n") ;
        printf ("    hi_j - %e > x_j > lo_j + %e => x_j free\n",
                Parm->bnd_pert, Parm->bnd_pert) ;
    }

    printf ("\n=================== SEARCH DIRECTION =====================\n") ;
    if ( Parm->FW == BLGTRUE )
    {
        printf ("Use the Frank-Wolfe search direction ................. FW: ") ;
        printf ("TRUE\n") ;
    }
    else if ( Parm->GP == BLGTRUE )
    {
        printf ("Use the gradient projection search direction ......... GP: ") ;
        printf ("TRUE\n") ;
    }
    else
    {
        printf ("Use the affine scaling search direction .............. AS: ") ;
        printf ("TRUE\n") ;
    }
    printf ("Parameter for replacing GP or AS with FW ....... lpfactor: %e\n",
             Parm->lpfactor) ;
    printf ("Set d_j to 0 when x_j nearly active:\n") ;
    printf ("    g_j > 0 and x_j - lo_j < %e => d_j = 0\n", Parm->lo_cut) ;
    printf ("    g_j < 0 and hi_j - x_j < %e => d_j = 0\n", Parm->hi_cut) ;

    printf ("\n=================== AFFINE SCALING CONTROL ===============\n") ;
    printf ("error tolerance for mu computation (B-A <= .. |g|) epsmu0: %e\n",
             Parm->epsmu0) ;
    printf ("bound for (B-A)/(|A|+|B|) in mu computation ...... epsmu1: %e\n",
             Parm->epsmu1) ;
    printf ("required decay factor in mu computation ........... gamma: %e\n",
             Parm->gamma) ;
    printf ("maximum number iterations when computing mu .. max_mu_its: %i\n",
             Parm->max_mu_its) ;
}
/*
Version 1.1 Change:
    1. Force a final projection on the feasible set, even when use_gp is TRUE.
       This ensures that the returned solution satisfies the constraints
       to within machine precision
    2. Correct bug in line 5266 (BLGnapdown). hi should have been lo
Version 1.2 Change:
    1. Fix some comment statement
*/
