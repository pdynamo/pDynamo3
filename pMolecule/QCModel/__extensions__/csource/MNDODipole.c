/*==================================================================================================================================
! . MNDO dipole integrals.
!=================================================================================================================================*/

# include <math.h>

# include "Integer.h"
# include "MNDODefinitions.h"
# include "MNDODipole.h"
# include "MNDOParameters.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole moment integrals - no dimension checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Inefficient but in same form as more general integral types. */
void MNDO_DipoleIntegrals ( const MNDOParametersContainer *parameters   ,
                            const IntegerArray1D          *basisIndices ,
                            const Coordinates3            *coordinates3 ,
                            const Vector3                 *center       ,
                                  SymmetricMatrix         *dX           ,
                                  SymmetricMatrix         *dY           ,
                                  SymmetricMatrix         *dZ           ,
                                  Status                  *status       )
{
    if ( ( parameters   != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( dX           != NULL ) &&
         ( dY           != NULL ) &&
         ( dZ           != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer         c, i, i0, nI ;
        auto MNDOParameters *iData ;
        auto Real            factor, h, x, xC = 0.0e+00, y, yC = 0.0e+00, z, zC = 0.0e+00 ;
        factor = 1.0e+00 / sqrt ( 3.0e+00 ) ;
        if ( center != NULL )
        {
            xC = Vector3_Item ( center, 0 ) ;
            yC = Vector3_Item ( center, 1 ) ;
            zC = Vector3_Item ( center, 2 ) ;
        }
        SymmetricMatrix_Set ( dX, 0.0e+00 ) ;
        SymmetricMatrix_Set ( dY, 0.0e+00 ) ;
        SymmetricMatrix_Set ( dZ, 0.0e+00 ) ;
        for ( i = 0 ; i < Coordinates3_Rows  ( coordinates3 ) ; i++ )
        {
            iData = parameters->entries[i] ;
            i0    = Array1D_Item ( basisIndices, i ) ;
            nI    = iData->norbitals ;
            Coordinates3_GetRow ( coordinates3, i, x, y, z ) ;
            for ( c = i0 ; c < ( i0 + nI ) ; c++ )
            {
                SymmetricMatrix_Item ( dX, c, c ) = ( x - xC ) ;
                SymmetricMatrix_Item ( dY, c, c ) = ( y - yC ) ;
                SymmetricMatrix_Item ( dZ, c, c ) = ( z - zC ) ;
            }
            if ( iData->norbitals >= 4 )
            {
                h = iData->ddp[1] ; /* . originally ddp[2] or dd */
                SymmetricMatrix_Item ( dX, i0 + PX, i0 ) = h ;
                SymmetricMatrix_Item ( dY, i0 + PY, i0 ) = h ;
                SymmetricMatrix_Item ( dZ, i0 + PZ, i0 ) = h ;
            }
            if ( iData->norbitals >= 9 )
            {
                h = iData->ddp[4] ; /* . originally ddp[5] */
                SymmetricMatrix_Item ( dX, i0 + DXZ,   i0 + PZ ) =   h ;
                SymmetricMatrix_Item ( dX, i0 + DX2Y2, i0 + PX ) =   h ;
                SymmetricMatrix_Item ( dX, i0 + DXY,   i0 + PY ) =   h ;
                SymmetricMatrix_Item ( dX, i0 + DZ2,   i0 + PX ) = - h * factor ;
                SymmetricMatrix_Item ( dY, i0 + DYZ,   i0 + PZ ) =   h ;
                SymmetricMatrix_Item ( dY, i0 + DX2Y2, i0 + PY ) = - h ;
                SymmetricMatrix_Item ( dY, i0 + DXY,   i0 + PX ) =   h ;
                SymmetricMatrix_Item ( dY, i0 + DZ2,   i0 + PY ) = - h * factor ;
                SymmetricMatrix_Item ( dZ, i0 + DXZ,   i0 + PX ) =   h ;
                SymmetricMatrix_Item ( dZ, i0 + DYZ,   i0 + PY ) =   h ;
                SymmetricMatrix_Item ( dZ, i0 + DZ2,   i0 + PZ ) =   h * factor * 2.0e+00 ;
            }
        }
    }
}
