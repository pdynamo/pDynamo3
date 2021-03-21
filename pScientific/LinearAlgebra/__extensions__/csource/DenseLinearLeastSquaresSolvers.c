/*==================================================================================================================================
! . Dense matrix linear least squares solvers.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "DenseLinearLeastSquaresSolvers.h"
# include "Memory.h"
# include "NumericalMacros.h"

# include "f2clapack.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _RelativeTolerance 1.0e-12

/*----------------------------------------------------------------------------------------------------------------------------------
! . General matrix solver.
! . rhs and solution can be the same if the matrix is square.
!---------------------------------------------------------------------------------------------------------------------------------*/
void LinearLeastSquaresSVDSolve ( RealArray2D *self              ,
                                  RealArray1D *rhs               ,
                                  Boolean      preserveInput     ,
                                  Real        *relativeTolerance ,
                                  RealArray1D *solution          ,
                                  Integer     *rank              ,
                                  Real        *condition         ,
                                  Status      *status            )
{
    if ( ( self != NULL ) && ( rhs != NULL ) && ( solution != NULL ) )
    {
        auto integer iFail = -2, m, n ;
        m = View2D_Rows    ( self ) ;
        n = View2D_Columns ( self ) ;
        /* . Check dimensions. */
        if ( ( m == View1D_Extent ( rhs ) ) && ( n == View1D_Extent ( solution ) ) )
        {
            auto integer      lWork ;
            auto Integer      eRank, i, mN, nS ;
            auto Real         cutOff, maxS, minS, tolerance, value, *work ;
            auto RealArray1D *s, slice, *temp ;
            auto RealArray2D *a, *ut, *v ;
            /* . Options. */
            if ( relativeTolerance == NULL ) tolerance = _RelativeTolerance   ;
            else                             tolerance = (*relativeTolerance) ;
            /* . Allocate space. */
            iFail = -1 ;
            mN    = Maximum ( m, n ) ;
            nS    = Maximum ( 1, Minimum ( m, n ) ) ;
            lWork = Maximum ( 3 * nS + mN, 5 * nS ) ;
            s     = RealArray1D_AllocateWithExtent  (             nS              , status ) ;
            temp  = RealArray1D_AllocateWithExtent  (             mN              , status ) ;
            ut    = RealArray2D_AllocateWithExtents ( ( Integer ) m, ( Integer ) m, status ) ;
            v     = RealArray2D_AllocateWithExtents ( ( Integer ) n, ( Integer ) n, status ) ;
            work  = Memory_AllocateArrayOfTypes ( ( Integer ) lWork, Real ) ;
            if ( preserveInput ) { a = RealArray2D_TransposeClone ( self, status ) ; } /* . Transpose as LAPACK expects column-major order. */
            else
            {
                if ( View2D_IsSquare ( a ) ) RealArray2D_TransposeSquare  ( a, status ) ;
                else                         RealArray2D_TransposeGeneral ( a, status ) ; /* . Will give error if not uniform. */
            }
            /* . Solve. */
            if ( Status_IsOK ( status ) )
            {
                dgesvd_ ( "A", "A", &m, &n, a->data, &m, s->data, ut->data, &m, v->data, &n, work, &lWork, &iFail ) ;
                /* . Treat the singular values. */
                maxS   = RealArray1D_AbsoluteMaximum ( s ) ;
                minS   = maxS ;
                cutOff = tolerance * maxS ;
                for ( eRank = i = 0 ; i < nS ; i++ )
                {
                    value = Array1D_Item ( s, i ) ;
                    if ( fabs ( value ) > cutOff )
                    {
                        Array1D_Item ( s, i ) = 1.0 / value ;
                        eRank ++ ;
                        minS = Minimum ( minS, fabs ( value ) ) ;
                    }
                    else
                    {
                        Array1D_Item ( s, i ) = 0.0e+00 ;
                        minS = cutOff ;
                    }
                }
                if ( condition != NULL ) (*condition) = maxS / minS ;
                if ( rank      != NULL ) (*rank)      = eRank       ;
                /* . Form V * S * U^T * RHS in temp. */
                RealArray1D_View           ( temp, 0, m , 1, False, &slice, NULL ) ;
                RealArray2D_VectorMultiply ( False, 1.0, ut, rhs, 0.0, &slice, NULL ) ;
                RealArray1D_View           ( temp, 0, nS, 1, False, &slice, NULL ) ;
                RealArray1D_Multiply       ( &slice, s, NULL ) ;
                RealArray1D_View           ( temp, 0, n , 1, False, &slice, NULL ) ;
                RealArray2D_VectorMultiply ( False, 1.0, v, &slice, 0.0, solution, NULL ) ;
            }
            /* . Finish up. */
            RealArray1D_Deallocate ( &s    ) ;
            RealArray1D_Deallocate ( &temp ) ;
            RealArray2D_Deallocate ( &ut   ) ;
            RealArray2D_Deallocate ( &v    ) ;
            Memory_Deallocate ( work ) ;
            if ( preserveInput ) RealArray2D_Deallocate ( &a ) ;
        }
        /* . Status. */
             if ( iFail == -2 ) { Status_Set ( status, Status_InvalidArgument ) ; }
        else if ( iFail == -1 ) { Status_Set ( status, Status_OutOfMemory     ) ; }
        else if ( iFail ==  0 ) { Status_Set ( status, Status_OK              ) ; }
        else                    { Status_Set ( status, Status_AlgorithmError  ) ; }
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/
# undef _RelativeTolerance
