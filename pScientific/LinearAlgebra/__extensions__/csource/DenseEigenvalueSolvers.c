/*==================================================================================================================================
! . Dense matrix eigenvalue solvers.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "DenseEigenvalueSolvers.h"
# include "Memory.h"
# include "Slice.h"

# include "f2clapack.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Symmetric matrix solver.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_EigenvaluesSolve (       SymmetricMatrix *self          ,
                                        const Boolean          preserveInput ,
                                        const Integer          lower         ,
                                        const Integer          upper         ,
                                              RealArray1D     *eigenValues   ,
                                              RealArray2D     *eigenVectors  ,
                                        const Boolean          isColumnMajor ,
                                              Status          *status        )
{
    if ( ( self != NULL ) && ( eigenValues != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean doEigenVectors ;
        auto Integer d, numberOfEigenValues, stride = 1 ;
        doEigenVectors = ( eigenVectors != NULL ) ;
        d              = self->extent ;
        /* . Check that lower and upper are within range. */
        numberOfEigenValues = SliceIndices_CheckSlice ( &lower, &upper, &stride, d, NULL, NULL, NULL, status ) ;
        if ( numberOfEigenValues > 0 )
        {
            /* . Various checks - self and eigenValues need to be compact as does extent1 of eigenVectors. */
                 if ( (                      ! View1D_IsCompact  ( eigenValues  ) ) ||
                      ( doEigenVectors  && ( ! View2D_IsCompact1 ( eigenVectors ) ) ) ) Status_Set ( status, Status_InvalidArrayOperation ) ;
            else if   ( doEigenVectors  &&
                        ( (     isColumnMajor   && ( ( View2D_Rows ( eigenVectors ) < numberOfEigenValues ) || ( View2D_Columns ( eigenVectors ) < d                   ) ) ) ||
                          ( ( ! isColumnMajor ) && ( ( View2D_Rows ( eigenVectors ) < d                   ) || ( View2D_Columns ( eigenVectors ) < numberOfEigenValues ) ) ) ) )
                                                                                                     Status_Set ( status, Status_NonConformableArrays ) ;
            else
            {
                auto char            *option ;
                auto integer          dummy, *iFail = NULL, il, info, iu, *iWork = NULL, m, n ;
                auto Real             absTol ;
                auto RealArray1D     *rWork = NULL ;
                auto SymmetricMatrix *work  = NULL ;
                /* . Initialization. */
                m = ( integer ) numberOfEigenValues ;
                n = ( integer ) d ;
                /* . Allocate space. */
                if ( doEigenVectors ) iFail = ( integer * ) calloc ( n, sizeof ( integer ) ) ;
                iWork = ( integer * ) calloc ( 5 * n, sizeof ( integer ) ) ;
                rWork = RealArray1D_AllocateWithExtent ( 8 * n, status ) ;
                if ( preserveInput ) work = SymmetricMatrix_CloneDeep ( self, status ) ;
                else                 work = self ;
                if ( ( ( doEigenVectors && ( iFail != NULL ) ) || ( ! doEigenVectors ) ) && ( iWork != NULL ) && ( rWork != NULL ) && ( work != NULL ) )
                {
                    absTol = 2.0e+00 * dlamch_ ( "S" ) ; /* . The recommended value. */
                    info   = 0 ;
                    /* . Set option flag. */
                    if ( m == d ) option = "A" ;
                    else          option = "I" ;
                    /* . Reset indices (Fortran convention). */
                    il = ( integer ) ( lower + 1 ) ;
                    iu = ( integer )   upper       ;
                    if ( doEigenVectors )
                    {
                        /* . LAPACK expects column-major order. */
                        /* . For non-square matrices massage the view parameters ready for the transpose later. */
                        if ( ( ! isColumnMajor ) && ( ! View2D_IsSquare ( eigenVectors ) ) )
                        {
                            auto Integer t ;
                            t                     = eigenVectors->extent0 ;
                            eigenVectors->extent0 = eigenVectors->extent1 ;
                            eigenVectors->extent1 = t                 ;
                            eigenVectors->stride0 = t * eigenVectors->stride1 ; /* . stride1 stays the same. */
                        }
                        dummy = ( integer ) eigenVectors->stride0 ;
                        dspevx_ ( "V"                               ,
                                  option                            ,
                                  "U"                               ,
                                  &n                                ,
                                  SymmetricMatrix_Data ( work )     ,
                                  NULL                              ,
                                  NULL                              ,
                                  &il                               ,
                                  &iu                               ,
                                  &absTol                           ,
                                  &m                                ,
                                  Array1D_Data ( eigenValues  )     ,
                                  Array2D_Data ( eigenVectors )     ,
                                  &dummy                            ,
                                  Array1D_Data ( rWork        )     ,
                                  iWork                             ,
                                  iFail                             ,
                                  &info                             ) ;
                        if ( ! isColumnMajor )
                        {
                            if ( View2D_IsSquare ( eigenVectors ) ) RealArray2D_TransposeSquare  ( eigenVectors, status ) ;
                            else                                    RealArray2D_TransposeGeneral ( eigenVectors, status ) ; /* . Will give error if not uniform. */
                        }
                    }
                    else
                    {
                        dummy = 1 ;
                        dspevx_ ( "N"                               ,
                                  option                            ,
                                  "U"                               ,
                                  &n                                ,
                                  SymmetricMatrix_Data ( work )     ,
                                  NULL                              ,
                                  NULL                              ,
                                  &il                               ,
                                  &iu                               ,
                                  &absTol                           ,
                                  &m                                ,
                                  Array1D_Data ( eigenValues )      ,
                                  NULL                              ,
                                  &dummy                            ,
                                  Array1D_Data ( rWork )            ,
                                  iWork                             ,
                                  NULL                              ,
                                  &info                             ) ;
                    }
                    if ( info != 0 ) { Status_Set ( status, Status_AlgorithmError ) ; printf ( "\nSymmetric Matrix Diagonalization Error = %d\n", info ) ; }
                }
                else Status_Set ( status, Status_OutOfMemory ) ;
                Memory_Deallocate      (  iFail ) ;
                Memory_Deallocate      (  iWork ) ;
                RealArray1D_Deallocate ( &rWork ) ;
                if ( preserveInput ) SymmetricMatrix_Deallocate ( &work ) ;
            }
        }
    }
}
