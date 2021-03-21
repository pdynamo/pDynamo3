/*==================================================================================================================================
! . Dense matrix linear equation solvers.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "DenseLinearEquationSolvers.h"
# include "Memory.h"

# include "cblas.h"
# include "f2clapack.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _RConditionTolerance 1.0e-15

/*----------------------------------------------------------------------------------------------------------------------------------
! . Square matrix solver.
! . rhs and solution can be the same.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SquareMatrix_LinearEquationsSolve ( RealArray2D *self          ,
                                         RealArray1D *rhs           ,
                                         Boolean      preserveInput ,
                                         RealArray1D *solution      ,
                                         Status      *status        )
{
    if ( ( self != NULL ) && ( rhs != NULL ) && ( solution != NULL ) )
    {
        auto integer iFail = -2, n ;
        n = View2D_Rows ( self ) ;
        if ( ( n == View2D_Columns ( self ) ) && ( n == View1D_Extent ( rhs ) ) && ( n == View1D_Extent ( solution ) ) )
        {
            auto char         equed ;
            auto integer     *iPiv, *iWork, nRHS ;
            auto Real        *af, *bErr, *c, *fErr, *r, rCond, *t, *work, *x ;
            auto RealArray2D *a ;
            /* . Allocate space. */
            iFail = -1 ;
            nRHS  =  1 ;
            iPiv  = ( integer * ) calloc ( n, sizeof ( integer ) ) ;
            iWork = ( integer * ) calloc ( n, sizeof ( integer ) ) ;
            af    = Memory_AllocateArrayOfTypes ( ( Integer ) ( n * n ), Real ) ;
            bErr  = Memory_AllocateArrayOfTypes ( ( Integer )       n  , Real ) ;
            c     = Memory_AllocateArrayOfTypes ( ( Integer )       n  , Real ) ;
            fErr  = Memory_AllocateArrayOfTypes ( ( Integer )       n  , Real ) ;
            r     = Memory_AllocateArrayOfTypes ( ( Integer )       n  , Real ) ;
            t     = Memory_AllocateArrayOfTypes ( ( Integer )       n  , Real ) ;
            work  = Memory_AllocateArrayOfTypes ( ( Integer )   4 * n  , Real ) ;
            x     = Memory_AllocateArrayOfTypes ( ( Integer )       n  , Real ) ;
            if ( preserveInput ) a = RealArray2D_CloneDeep ( self, NULL ) ;
            else                 a = self ;
            /* . Solve. */
            if ( ( iPiv != NULL ) && ( iWork != NULL ) && ( af != NULL ) && ( bErr != NULL ) && ( c != NULL ) &&
                 ( fErr != NULL ) && ( r     != NULL ) && ( t  != NULL ) && ( work != NULL ) && ( x != NULL ) && ( a != NULL ) )
            {
	        /* . Transfer rhs to t (needed as no vector increment allowed). */
	        cblas_dcopy ( n, rhs->data, rhs->stride, t, 1 ) ;
                /* . Solve for the transpose as LAPACK expects column-major order. */
                dgesvx_ ( "E", "T", &n, &nRHS, a->data, &n, af, &n, iPiv, &equed, r, c, t, &n, x, &n, &rCond, fErr, bErr, work, iWork, &iFail ) ;
                /* . Copy solution to rhs. */
	        cblas_dcopy ( n, x, 1, solution->data, solution->stride ) ;
                /* . Check the condition number. */
                if ( fabs ( rCond ) < _RConditionTolerance ) iFail = n + 1 ;
            }
            /* . Finish up. */
            Memory_Deallocate ( af    ) ;
            Memory_Deallocate ( bErr  ) ;
            Memory_Deallocate ( c     ) ;
            Memory_Deallocate ( fErr  ) ;
            Memory_Deallocate ( iPiv  ) ;
            Memory_Deallocate ( iWork ) ;
            Memory_Deallocate ( r     ) ;
            Memory_Deallocate ( t     ) ;
            Memory_Deallocate ( work  ) ;
            Memory_Deallocate ( x     ) ;
            if ( preserveInput ) RealArray2D_Deallocate ( &a ) ;
        }
        /* . Status. */
             if ( iFail == -2 ) { Status_Set ( status, Status_InvalidArgument ) ; }
        else if ( iFail == -1 ) { Status_Set ( status, Status_OutOfMemory     ) ; }
        else if ( iFail ==  0 ) { Status_Set ( status, Status_OK              ) ; }
        else                    { Status_Set ( status, Status_AlgorithmError  ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Symmetric matrix solver.
! . rhs and solution can be the same.
! . self is preserved on output.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_LinearEquationsSolve ( SymmetricMatrix *self     ,
                                            RealArray1D     *rhs      ,
                                            RealArray1D     *solution ,
                                            Status          *status   )
{
    if ( ( self != NULL ) && ( rhs != NULL ) && ( solution != NULL ) )
    {
        auto integer iFail = -2, n   ;
        n = self->extent ;
        /* . Check dimensions. */
        if ( ( n == View1D_Extent ( rhs ) ) && ( n == View1D_Extent ( solution ) ) )
        {
            auto integer *iPiv, *iWork, nRHS ;
            auto Real    *afp, *bErr, *fErr, rCond, *t, *work, *x ;
            /* . Allocate space. */
            iFail = -1 ;
            nRHS  =  1 ;
            iPiv  = ( integer * ) calloc ( n, sizeof ( integer ) ) ;
            iWork = ( integer * ) calloc ( n, sizeof ( integer ) ) ;
            afp   = Memory_AllocateArrayOfTypes ( ( Integer ) ( ( n * ( n + 1 ) ) / 2 ), Real ) ;
            bErr  = Memory_AllocateArrayOfTypes ( ( Integer )       n  , Real ) ;
            fErr  = Memory_AllocateArrayOfTypes ( ( Integer )       n  , Real ) ;
            work  = Memory_AllocateArrayOfTypes ( ( Integer ) ( 3 * n ), Real ) ;
            x     = Memory_AllocateArrayOfTypes ( ( Integer )       n  , Real ) ;
	    /* . Transfer rhs to t (needed as no vector increment allowed). */
	    if ( rhs->stride == 1 )
	    {
	        t = rhs->data ;
            }
	    else
	    {
                t = Memory_AllocateArrayOfTypes ( ( Integer ) n, Real ) ;
   	        cblas_dcopy ( n, rhs->data, rhs->stride, t, 1 ) ;
	    }
            /* . Solve. */
            if ( ( iPiv != NULL ) && ( iWork != NULL ) && ( afp != NULL ) && ( bErr != NULL ) && ( fErr != NULL ) && ( work != NULL ) && ( x != NULL ) && ( t != NULL ) )
            {
                dspsvx_ ( "N", "U", &n, &nRHS, self->data, afp, iPiv, t, &n, x, &n, &rCond, fErr, bErr, work, iWork, &iFail ) ;
                /* . Copy solution to solution. */
	        cblas_dcopy ( n, x, 1, solution->data, solution->stride ) ;
                /* . Check the condition number. */
                if ( fabs ( rCond ) < _RConditionTolerance ) iFail = n + 1 ;
            }
            /* . Finish up. */
            Memory_Deallocate ( afp   ) ;
            Memory_Deallocate ( bErr  ) ;
            Memory_Deallocate ( fErr  ) ;
            Memory_Deallocate ( iPiv  ) ;
            Memory_Deallocate ( iWork ) ;
            Memory_Deallocate ( work  ) ;
            Memory_Deallocate ( x     ) ;
	    if ( rhs->stride != 1 ) Memory_Deallocate ( t ) ;
        }
        /* . Status. */
             if ( iFail == -2 ) { Status_Set ( status, Status_InvalidArgument ) ; }
        else if ( iFail == -1 ) { Status_Set ( status, Status_OutOfMemory     ) ; }
        else if ( iFail ==  0 ) { Status_Set ( status, Status_OK              ) ; }
        else                    { Status_Set ( status, Status_AlgorithmError  ) ; }
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/
# undef _RConditionTolerance
