/* ====================================================
 * CUTEr interface for BLG     Oct. 10, 2010
 *
 * H. Zhang & W. Hager
 *
 * (Based on gencma.c of D. Orban, Feb 3, 2003)
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BLGMA

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "../../include/cuter.h"
#include "../../include/BLG_user.h"

#ifdef Isg95
#define MAINENTRY MAIN_
#else
#define MAINENTRY main
#endif

/* prototypes */
void BLG_value
(
    BLGobjective *user
) ;

void BLG_grad
(
    BLGobjective *user
) ;

void BLG_valgrad
(
    BLGobjective *user
) ;

/* global variables */
    integer CUTEr_nvar;        /* number of variables */
    integer CUTEr_ncon;        /* number of constraints */
    integer CUTEr_nnzj;        /* number of nonzeros in Jacobian */
    integer CUTEr_nnzh;        /* number of nonzeros in upper triangular
                                  part of the Hessian of the Lagrangian */

/* main program */
    int MAINENTRY( void ) {

        char *fname = "OUTSDIF.d"; /* CUTEr data file */
        integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
        integer iout = 6;          /* FORTRAN unit number for error output */
        integer ierr;              /* Exit flag from OPEN and CLOSE */

        VarTypes vtypes;

        integer    *indvar = NULL, *indfun = NULL, ncon_dummy ;
        doublereal *x, *lo, *hi, *dummy1, *dummy2 ;
        doublereal *v = NULL, *cl = NULL, *cu = NULL, *c = NULL, *cjac = NULL;
        logical    *equatn = NULL, *linear = NULL;
        char       *pname, *vnames, *gnames;
        logical     efirst = TRUE_, lfirst = TRUE_, nvfrst = FALSE_, grad;
        logical     constrained = FALSE_;

        real        calls[7], cpu[2];
        integer     nlin = 0, nbnds = 0, neq = 0;
        doublereal  dummy;
        integer     ExitCode;
        int         i;

        BLGparm Parm ;
        int n, status ;
        doublereal *a, b, *xtemp ;
        long int index ;


        /* Open problem description file OUTSDIF.d */
        ierr = 0;
        FORTRAN_OPEN( &funit, fname, &ierr );
        if( ierr ) {
            printf("Error opening file OUTSDIF.d.\nAborting.\n");
            exit(1);
        }

        /* Determine problem size */
        CDIMEN( &funit, &CUTEr_nvar, &CUTEr_ncon );

        /* Determine whether to call constrained or unconstrained tools */
        if( CUTEr_ncon ) constrained = TRUE_;

        /* Seems to be needed for some Solaris C compilers */
        ncon_dummy = CUTEr_ncon + 1;

        /* Reserve memory for variables, bounds, and multipliers */
        /* and call appropriate initialization routine for CUTEr */
        MALLOC( x,      CUTEr_nvar, doublereal );
        MALLOC( lo,     CUTEr_nvar, doublereal );
        MALLOC( hi,     CUTEr_nvar, doublereal );
        if( constrained ) {
            MALLOC( equatn, CUTEr_ncon+1, logical    );
            MALLOC( linear, CUTEr_ncon+1, logical    );
            MALLOC( v,      CUTEr_ncon+CUTEr_nvar+1, doublereal );
            MALLOC( cl,     CUTEr_ncon+1, doublereal );
            MALLOC( cu,     CUTEr_ncon+1, doublereal );
            CSETUP( &funit, &iout, &CUTEr_nvar, &CUTEr_ncon, x, lo, hi,
                    &CUTEr_nvar, equatn, linear, v, cl, cu, &ncon_dummy,
                    &efirst, &lfirst, &nvfrst );
            if ( CUTEr_ncon > 1)
            {
                printf(" There is more than one constraint! \n") ;
                exit(1) ;
            }
            else if ( ! linear[0])
            {
                printf(" The constraint is not linear! \n") ;
                exit(1) ;
            }
            else if ( ! equatn[0])
            {
                printf(" The constraint is not equality constraint! \n") ;
                exit(1) ;
            }
        } else {
            MALLOC( equatn, 1, logical    );
            MALLOC( linear, 1, logical    );
            MALLOC( cl, 1, doublereal );
            MALLOC( cu, 1, doublereal );
            MALLOC( v, CUTEr_nvar+1, doublereal );
            USETUP( &funit, &iout, &CUTEr_nvar, x, lo, hi, &CUTEr_nvar );
        }

        /* Get problem, variables and constraints names */
        MALLOC( pname, FSTRING_LEN+1, char );
        MALLOC( vnames, CUTEr_nvar*FSTRING_LEN, char );
        if( constrained ) {
           MALLOC( gnames, CUTEr_ncon*FSTRING_LEN, char );
           CNAMES( &CUTEr_nvar, &CUTEr_ncon, pname, vnames, gnames );
           FREE( gnames );
        } else {
           UNAMES( &CUTEr_nvar, pname, vnames );
        }
        FREE( vnames );

        /* Make sure to null-terminate problem name */
        pname[FSTRING_LEN] = '\0';
        i = FSTRING_LEN - 1;
        while( i-- > 0 && pname[i] == ' ') {
            pname[i] = '\0';
        }

/*  Set up the problem:
                     min f(x)

     subject to

             a' x = b,   a is n by 1 vector
             lo <= x <= hi
*/
        n = CUTEr_nvar ;
        if ( constrained ) /* there is one linear constraint*/
        {
           MALLOC( a,  CUTEr_nvar, doublereal );
           MALLOC( xtemp,  CUTEr_nvar, doublereal );
           for (i=0; i < CUTEr_nvar; i++) xtemp [i] = 0. ;
           grad = TRUE_ ;
           index = 1 ;
           CCIFG(&CUTEr_nvar, &index, xtemp, &b, a, &grad) ; /* get a */
           b = cl[0] - b ; /* get right hand side b; note cl[0]=cu[0] since equality constraint */
        }
        else
        {
           a = NULL ;
           b = 0. ;
        }
        /* Set parameters */
        BLGdefault(&Parm) ;
        Parm.PrintFinal = BLGTRUE ;
        /*  Parm.Subspace = BLGTRUE ;
            Parm.nsub = 10 ;*/
        /*  use affine scaling search direction, set the following Parm.GP = 1 */
        /*  Parm.GP = 1 ; */

        /* Call the optimizer, with some loss of efficiency, you could omit the valgrad routine */
        BLG (x, NULL, lo, hi, a, b, NULL, n, 1.e-6, NULL, NULL, NULL, &Parm,
             NULL, BLG_value, BLG_grad, BLG_valgrad) ;

        ExitCode = 0;

        /* Get CUTEr statistics */
        CREPRT( calls, cpu );

        printf ("status: %5i CPU: %8.2f\n", status, cpu [1]) ;
/*
        printf("\n\n ************************ CUTEr statistics ************************\n\n");
        printf(" Code used               : BLG\n");
        printf(" Problem                 : %-s\n", pname);
        printf(" # variables             = %-10d\n", CUTEr_nvar);
        printf(" # bound constraints     = %-10d\n", vtypes.nbnds);
        printf(" # objective functions   = %-15.7g\n", calls[0]);
        printf(" # objective gradients   = %-15.7g\n", calls[1]);
        printf(" # objective Hessians    = %-15.7g\n", calls[2]);
        printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]);
        printf(" Exit code               = %-10d\n", ExitCode);
        printf(" Final f                 = %-15.7g\n",dummy);
        printf(" Set up time             = %-10.2f seconds\n", cpu[0]);
        printf(" Solve time              = %-10.2f seconds\n", cpu[1]);
        printf(" ******************************************************************\n\n");
*/

        ierr = 0;
        FORTRAN_CLOSE( &funit, &ierr );
        if( ierr ) {
            printf( "Error closing %s on unit %d.\n", fname, funit );
            printf( "Trying not to abort.\n" );
        }

        /* Free workspace */
        FREE( pname );
        FREE( x ); FREE( lo ); FREE( hi );
        FREE( v ); FREE( cl ); FREE( cu );
        FREE(equatn) ; FREE(linear) ; FREE(a) ;

        return 0;
    }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
void BLG_value
(
    BLGobjective *user
)
{
    long int N ;
    double f ;
    N = (long int) user->n ;
    UFN (&N, user->x, &f) ;
    user->f = f ;
    return ;
}

void BLG_grad
(
    BLGobjective *user
)
{
    long int N ;
    N = (long int) user->n ;
    UGR (&N, user->x, user->g) ;
    return ;
}

void BLG_valgrad
(
    BLGobjective *user
)
{
    long int N ;
    long int grad ;
    double f ;
    N = (long int) user->n ;
    grad = (long int) 1 ;
    UOFG (&N, user->x, &f, user->g, &grad) ;
    user->f = f ;
    return ;
}
