/*==================================================================================================================================
! . Container integrals - 1 basis, 1 electron.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisContainerIntegrals_b1e1n0.h"
# include "GaussianBasisIntegrals_b1e1n0.h"
# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the self-overlap integrals.
! . SelfOverlap is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_SelfOverlap ( const GaussianBasisContainer *self         ,
                                                   const IntegerArray1D         *basisIndices ,
                                                         RealArray1D            *selfOverlap  ,
                                                         Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( basisIndices != NULL ) &&
         ( selfOverlap  != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( View1D_Extent ( selfOverlap ) == Array1D_Item ( basisIndices, self->capacity ) )
        {
            auto Integer     i, start, stop ;
            auto RealArray1D view ;
            RealArray1D_Set ( selfOverlap, 0.0e+00 ) ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                start = Array1D_Item ( basisIndices, i   ) ;
                stop  = Array1D_Item ( basisIndices, i+1 ) ;
                RealArray1D_View ( selfOverlap, start, ( stop - start ), 1, False, &view, NULL ) ;
                GaussianBasisIntegrals_SelfOverlap ( self->entries[i], &view ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}
