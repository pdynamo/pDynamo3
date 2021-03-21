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
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "RealIterator.h"

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
