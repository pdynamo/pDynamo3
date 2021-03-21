/*==================================================================================================================================
! . Container integrals - 1 basis, 0 electrons, 1 nucleus/point.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "GaussianBasisContainerIntegrals_b1e0n1.h"
# include "GaussianBasisIntegrals_b1e0n1.h"
# include "Integer.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at grid points.
! . Values is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_Grid ( const GaussianBasisContainer *self         ,
                                            const IntegerArray1D         *basisIndices ,
                                            const Coordinates3           *coordinates3 ,
                                            const Coordinates3           *rG           ,
                                                  RealArray2D            *values       ,
                                                  Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( rG           != NULL ) &&
         ( values       != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer  b, g ;
        b = View2D_Rows    ( values ) ; /* . N * G. */
        g = View2D_Columns ( values ) ;
        if ( ( Coordinates3_Rows ( rG                           ) == g              ) &&
             ( Coordinates3_Rows ( coordinates3                 ) == self->capacity ) &&
             ( Array1D_Item      ( basisIndices, self->capacity ) == b              ) )
        {
            auto Integer      i, start, stop ;
            auto Real        *rI ;
            auto RealArray2D  view ;
            RealArray2D_Set ( values, 0.0e+00 ) ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                rI    = Coordinates3_RowPointer ( coordinates3, i   ) ;
                start = Array1D_Item    ( basisIndices, i   ) ;
                stop  = Array1D_Item    ( basisIndices, i+1 ) ;
                RealArray2D_View ( values, start, 0, ( stop - start ), g, 1, 1, False, &view, NULL ) ;
                GaussianBasisIntegrals_Grid ( self->entries[i] ,
                                              rI               ,
                                              rG               ,
                                              &view            ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions and, optionally, their derivatives at grid points.
! . The results are put in a grid function data block.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_GridFunctionDataBlock ( const GaussianBasisContainer *self         ,
                                                             const IntegerArray1D         *basisIndices ,
                                                             const Coordinates3           *coordinates3 ,
                                                             const Coordinates3           *rG           ,
                                                             const Boolean                 resize       ,
                                                             const Real                   *tolerance    ,
                                                                   GridFunctionDataBlock  *data         ,
                                                                   Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( rG           != NULL ) &&
         ( data         != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer  b, g ;
        b = data->numberOfFunctions ; /* . N * G. */
        g = data->numberOfPoints    ;
        if ( ( Coordinates3_Rows ( rG                           ) <= g              ) &&
             ( Coordinates3_Rows ( coordinates3                 ) == self->capacity ) &&
             ( Array1D_Item      ( basisIndices, self->capacity ) <= b              ) )
        {
            auto GaussianBasis *iBasis ;
            auto Integer        f0, i, nF ;
            auto Real          *rI ;
            auto RealArray2D    f, fX, fY, fZ, fXX, fXY, fXZ, fYY, fYZ, fZZ,
                                fXXX, fXXY, fXXZ, fXYY, fXYZ, fXZZ, fYYY, fYYZ, fYZZ, fZZZ ;
            GridFunctionDataBlock_Initialize ( data ) ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                rI     = Coordinates3_RowPointer ( coordinates3, i   ) ;
                f0     = Array1D_Item            ( basisIndices, i   ) ;
                nF     = Array1D_Item            ( basisIndices, i+1 ) - f0 ;
                RealArray2D_View ( data->f, f0, 0, nF, g, 1, 1, False, &f, NULL ) ;
                if ( data->order == 0 ) GaussianBasisIntegrals_Grid ( iBasis, rI, rG, &f ) ;
                else
                {
                    RealArray2D_View ( data->fX, f0, 0, nF, g, 1, 1, False, &fX, NULL ) ;
                    RealArray2D_View ( data->fY, f0, 0, nF, g, 1, 1, False, &fY, NULL ) ;
                    RealArray2D_View ( data->fZ, f0, 0, nF, g, 1, 1, False, &fZ, NULL ) ;
                    if ( data->order == 1 ) GaussianBasisIntegrals_GridD ( iBasis, rI, rG, &f, &fX, &fY, &fZ ) ;
                    else
                    {
                        RealArray2D_View ( data->fXX, f0, 0, nF, g, 1, 1, False, &fXX, NULL ) ;
                        RealArray2D_View ( data->fXY, f0, 0, nF, g, 1, 1, False, &fXY, NULL ) ;
                        RealArray2D_View ( data->fXZ, f0, 0, nF, g, 1, 1, False, &fXZ, NULL ) ;
                        RealArray2D_View ( data->fYY, f0, 0, nF, g, 1, 1, False, &fYY, NULL ) ;
                        RealArray2D_View ( data->fYZ, f0, 0, nF, g, 1, 1, False, &fYZ, NULL ) ;
                        RealArray2D_View ( data->fZZ, f0, 0, nF, g, 1, 1, False, &fZZ, NULL ) ;
                        if ( data->order == 2 ) GaussianBasisIntegrals_GridD2 ( iBasis, rI, rG, &f, &fX, &fY, &fZ, &fXX, &fXY, &fXZ, &fYY, &fYZ, &fZZ ) ;
                        else
                        {
                            RealArray2D_View ( data->fXXX, f0, 0, nF, g, 1, 1, False, &fXXX, NULL ) ;
                            RealArray2D_View ( data->fXXY, f0, 0, nF, g, 1, 1, False, &fXXY, NULL ) ;
                            RealArray2D_View ( data->fXXZ, f0, 0, nF, g, 1, 1, False, &fXXZ, NULL ) ;
                            RealArray2D_View ( data->fXYY, f0, 0, nF, g, 1, 1, False, &fXYY, NULL ) ;
                            RealArray2D_View ( data->fXYZ, f0, 0, nF, g, 1, 1, False, &fXYZ, NULL ) ;
                            RealArray2D_View ( data->fXZZ, f0, 0, nF, g, 1, 1, False, &fXZZ, NULL ) ;
                            RealArray2D_View ( data->fYYY, f0, 0, nF, g, 1, 1, False, &fYYY, NULL ) ;
                            RealArray2D_View ( data->fYYZ, f0, 0, nF, g, 1, 1, False, &fYYZ, NULL ) ;
                            RealArray2D_View ( data->fYZZ, f0, 0, nF, g, 1, 1, False, &fYZZ, NULL ) ;
                            RealArray2D_View ( data->fZZZ, f0, 0, nF, g, 1, 1, False, &fZZZ, NULL ) ;
                            GaussianBasisIntegrals_GridD3 ( iBasis, rI, rG, &f, &fX, &fY, &fZ, &fXX, &fXY, &fXZ, &fYY, &fYZ, &fZZ ,
                                                             &fXXX, &fXXY, &fXXZ, &fXYY, &fXYZ, &fXZZ, &fYYY, &fYYZ, &fYZZ, &fZZZ ) ;
                        }
                    }
                }
                /* . Increment the function number. */
                data->numberOfFunctions += iBasis->nbasisw ;
            }
            if ( ( resize ) && ( tolerance != NULL ) && ( (*tolerance) > 0.0e+00 ) )
            {
                GridFunctionDataBlock_FilterValues ( data, 0, tolerance ) ;
                GridFunctionDataBlock_Resize ( data, data->numberOfFunctions, status ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}
