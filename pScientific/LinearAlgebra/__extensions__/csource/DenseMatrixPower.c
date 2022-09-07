/*==================================================================================================================================
! . Dense matrix power.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "Array_Macros.h"
# include "DenseEigenvalueSolvers.h"
# include "DenseMatrixPower.h"
# include "Integer.h"
# include "Iterator.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "RealIterator.h"

# include "cblas.h"
# include "f2clapack.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _RelativeTolerance 1.0e-12

/*----------------------------------------------------------------------------------------------------------------------------------
! . General matrix pseudo-inverse.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Matrix_PseudoInverse (       RealArray2D *self              ,
                            const Boolean      preserveInput     ,
                                  Real        *relativeTolerance ,
                                  RealArray2D *inverse           ,
                                  Integer     *rank              ,
                                  Real        *condition         ,
                                  Status      *status            )
{
    if ( ( self != NULL ) && ( inverse != NULL ) )
    {
        auto integer iFail = -2, m, n ;
        m = View2D_Rows    ( self ) ;
        n = View2D_Columns ( self ) ;
        /* . Inverse has the same dimensions as the transpose of self. */
        if ( ( m == View2D_Columns ( inverse ) ) && ( n == View2D_Rows ( inverse ) ) )
        {
            auto integer      lWork ;
            auto Integer      eRank, i, mN, nS ;
            auto Real         cutOff, maxS, minS, tolerance, value, *work ;
            auto RealArray1D *s, vView ;
            auto RealArray2D *a, mView, *ut, *v ;
            /* . Options. */
            if ( relativeTolerance == NULL ) tolerance = _RelativeTolerance   ;
            else                             tolerance = (*relativeTolerance) ;
            /* . Allocate space. */
            iFail = -1 ;
            mN    = Maximum ( m, n ) ;
            nS    = Maximum ( 1, Minimum ( m, n ) ) ;
            lWork = Maximum ( 3 * nS + mN, 5 * nS ) ;
            s     = RealArray1D_AllocateWithExtent  (             nS              , status ) ;
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
                /* . Get the SVD of the matrix = U * S * V^T. */
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
                /* . Form V * S * U^T in inverse. */
                /* . ( n * n ) * ( n * m ) * ( m * m ). */
                if ( n >= m )
                {
                    /* . S * U^T first then premultiply by V. */
                    for ( i = 0 ; i < m ; i++ )
                    {
                        RealArray2D_RowView ( ut, i, False, &vView, NULL ) ;
                        RealArray1D_Scale   ( &vView, Array1D_Item ( s, i ) ) ;
                    }
                    if ( n == m )
                    {
                        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, v, ut, 0.0, inverse, NULL ) ;
                    }
                    else
                    {
                        RealArray2D_View ( v, 0, 0, n, m, 1, 1, False, &mView, NULL ) ;
                        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, &mView, ut, 0.0, inverse, NULL ) ;
                    }
                }
                else
                {
                    /* . V * S first and then postmultiply by U^T. */
                    for ( i = 0 ; i < n ; i++ )
                    {
                        RealArray2D_ColumnView ( v, i, False, &vView, NULL ) ;
                        RealArray1D_Scale      ( &vView, Array1D_Item ( s, i ) ) ;
                    }
                    RealArray2D_View ( ut, 0, 0, n, m, 1, 1, False, &mView, NULL ) ;
                    RealArray2D_MatrixMultiply ( False, False, 1.0e+00, v, &mView, 0.0, inverse, NULL ) ;
                }
            }
            /* . Finish up. */
            RealArray1D_Deallocate ( &s    ) ;
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

/*----------------------------------------------------------------------------------------------------------------------------------
! . Symmetric matrix inverse power.
! . Power is positive!
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_InversePower (       SymmetricMatrix *self          , 
                                    const Boolean          preserveInput ,
                                    const Real             power         ,
                                    const Real             tolerance     ,
                                          SymmetricMatrix *result        ,
                                          Status          *status        )
{
    if ( ( self != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        if ( power > 0.0e+00 )
        {
            auto Integer n = SymmetricMatrix_Extent ( self ) ;
            if ( n == SymmetricMatrix_Extent ( result ) )
            {
                Iterator        *iterator     ;
                RealArray1D     *eigenValues  ;
                RealArray2D     *eigenVectors ;
                SymmetricMatrix *target       ;
                if ( preserveInput ) target = SymmetricMatrix_CloneDeep ( self, status ) ;
                else                 target = self ;
                eigenValues  = RealArray1D_AllocateWithExtent  ( n,    status ) ;
                eigenVectors = RealArray2D_AllocateWithExtents ( n, n, status ) ;
                SymmetricMatrix_EigenvaluesSolve ( target, preserveInput, 0, n, eigenValues, eigenVectors, False, status ) ;
                iterator = View1D_MakeIterator ( ( View1D * ) eigenValues, status ) ;
                if ( power == 1.0e+00 ) RealIterator_Reciprocate      ( iterator, Array_BlockDataPointer ( eigenValues, eigenValues->offset ),        tolerance, 0.0e+00, status ) ;
                else                    RealIterator_ReciprocatePower ( iterator, Array_BlockDataPointer ( eigenValues, eigenValues->offset ), power, tolerance, 0.0e+00, status ) ;
                SymmetricMatrix_MakeFromEigensystem ( result, True, n, eigenValues, eigenVectors, status ) ;
                if ( preserveInput ) SymmetricMatrix_Deallocate ( &target ) ;
                Iterator_Deallocate    ( &iterator     ) ;
                RealArray1D_Deallocate ( &eigenValues  ) ;
                RealArray2D_Deallocate ( &eigenVectors ) ;
            }
            else Status_Set ( status, Status_NonConformableArrays ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Symmetric matrix power.
! . Power can take any value.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Power (       SymmetricMatrix *self          , 
                             const Boolean          preserveInput ,
                             const Real             power         ,
                                   SymmetricMatrix *result        ,
                                   Status          *status        )
{
    if ( ( self != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = SymmetricMatrix_Extent ( self ) ;
        if ( n == SymmetricMatrix_Extent ( result ) )
        {
            Iterator        *iterator     ;
            RealArray1D     *eigenValues  ;
            RealArray2D     *eigenVectors ;
            SymmetricMatrix *target       ;
            if ( preserveInput ) target = SymmetricMatrix_CloneDeep ( self, status ) ;
            else                 target = self ;
            eigenValues  = RealArray1D_AllocateWithExtent  ( n,    status ) ;
            eigenVectors = RealArray2D_AllocateWithExtents ( n, n, status ) ;
            SymmetricMatrix_EigenvaluesSolve ( target, preserveInput, 0, n, eigenValues, eigenVectors, False, status ) ;
            iterator = View1D_MakeIterator ( ( View1D * ) eigenValues, status ) ;
            RealIterator_Power ( iterator, Array_BlockDataPointer ( eigenValues, eigenValues->offset ), power, status ) ;
            SymmetricMatrix_MakeFromEigensystem ( result, True, n, eigenValues, eigenVectors, status ) ;
            if ( preserveInput ) SymmetricMatrix_Deallocate ( &target ) ;
            Iterator_Deallocate    ( &iterator     ) ;
            RealArray1D_Deallocate ( &eigenValues  ) ;
            RealArray2D_Deallocate ( &eigenVectors ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/
# undef _RelativeTolerance
