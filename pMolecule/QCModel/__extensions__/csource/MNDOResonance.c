/*==================================================================================================================================
! . MNDO resonance interactions.
!=================================================================================================================================*/

# include <math.h>

# include "GaussianBasis.h"
# include "GaussianBasisIntegrals_b2e1n0.h"
# include "GaussianBasisTransform.h"
# include "Integer.h"
# include "MNDOResonance.h"
# include "MNDOParameters.h"
# include "Real.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . The resonance gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDO_ResonanceGradients ( const MNDOParametersContainer *parameters   ,
                               const GaussianBasisContainer  *bases        ,
                               const IntegerArray1D          *basisIndices ,
                               const Coordinates3            *coordinates3 ,
                               const SymmetricMatrix         *dTotal       ,
                                     Coordinates3            *gradients3   )
{
    if ( ( parameters   != NULL ) &&
         ( bases        != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( dTotal       != NULL ) &&
         ( gradients3   != NULL ) )
    {
        auto Integer         i, i0, j, j0, nI, nJ, u, uv, v ;
        auto GaussianBasis  *iBasis, *jBasis ;
        auto MNDOParameters *iData , *jData  ;
        auto Real            b, gX, gY, gZ, *xI, *xJ ;
        auto RealArray2D    *sX, *sY, *sZ ;
        for ( i = 0 ; i < Coordinates3_Rows ( coordinates3 ) ; i++ )
        {
            iBasis = bases->entries[i] ;
            iData  = parameters->entries[i] ;
            i0     = Array1D_Item ( basisIndices, i ) ;
            nI     = iData->norbitals ;
            xI     = Coordinates3_RowPointer ( coordinates3, i ) ;
            for ( j = 0 ; j < i ; j++ )
            {
                jBasis = bases->entries[j] ;
                jData  = parameters->entries[j] ;
                j0     = Array1D_Item ( basisIndices, j ) ;
                nJ     = jData->norbitals ;
                xJ     = Coordinates3_RowPointer ( coordinates3, j ) ;
                sX     = RealArray2D_AllocateWithExtents ( iBasis->nbasisw, jBasis->nbasisw, NULL ) ;
                sY     = RealArray2D_AllocateWithExtents ( iBasis->nbasisw, jBasis->nbasisw, NULL ) ;
                sZ     = RealArray2D_AllocateWithExtents ( iBasis->nbasisw, jBasis->nbasisw, NULL ) ;
                GaussianBasisIntegrals_2OverlapD ( iBasis, xI, jBasis, xJ, sX, sY, sZ ) ;
                GaussianBasis_TransformIntegrals2 ( &sX, iBasis->c2o, jBasis->c2o ) ;
                GaussianBasis_TransformIntegrals2 ( &sY, iBasis->c2o, jBasis->c2o ) ;
                GaussianBasis_TransformIntegrals2 ( &sZ, iBasis->c2o, jBasis->c2o ) ;
                if ( ( sX != NULL ) && ( sY != NULL ) && ( sZ != NULL ) )
                {
                    gX = gY = gZ = 0.0e+00 ;
                    for ( u = uv = 0 ; u < nI ; u++ )
                    {
                        uv = ( ( u + i0 ) * ( ( u + i0 ) + 1 ) ) / 2 + j0 ;
                        for ( v = 0 ; v < nJ ; uv++, v++ )
                        {
                            b   = ( iData->beta[u] + jData->beta[v] ) * dTotal->data[uv] ; /* . Note the implicit factor of 2 here. */
                            gX += b * Array2D_Item ( sX, u, v ) ;
                            gY += b * Array2D_Item ( sY, u, v ) ;
                            gZ += b * Array2D_Item ( sZ, u, v ) ;
                        }
                    }
                    Coordinates3_IncrementRow ( gradients3, i, gX, gY, gZ ) ;
                    Coordinates3_DecrementRow ( gradients3, j, gX, gY, gZ ) ;
                }
                RealArray2D_Deallocate ( &sX ) ;
                RealArray2D_Deallocate ( &sY ) ;
                RealArray2D_Deallocate ( &sZ ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The resonance integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDO_ResonanceIntegrals ( const MNDOParametersContainer *parameters        ,
                               const GaussianBasisContainer  *bases             ,
                               const IntegerArray1D          *basisIndices      ,
                               const Coordinates3            *coordinates3      ,
                                     SymmetricMatrix         *oneElectronMatrix )
{
    if ( ( parameters        != NULL ) &&
         ( bases             != NULL ) &&
         ( basisIndices      != NULL ) &&
         ( coordinates3      != NULL ) &&
         ( oneElectronMatrix != NULL ) )
    {
        auto Integer         i, i0, j, j0, nI, nJ, u, v ;
        auto GaussianBasis  *iBasis, *jBasis ;
        auto MNDOParameters *iData , *jData  ;
        auto Real           *xI, *xJ ;
        auto RealArray2D    *s ;
        for ( i = 0 ; i < Coordinates3_Rows ( coordinates3 ) ; i++ )
        {
            iBasis = bases->entries[i] ;
            iData  = parameters->entries[i] ;
            i0     = Array1D_Item ( basisIndices, i ) ;
            nI     = iData->norbitals ;
            xI     = Coordinates3_RowPointer ( coordinates3, i ) ;
            for ( u = 0 ; u < nI ; u++ ) SymmetricMatrix_Item ( oneElectronMatrix, u+i0, u+i0 ) += iData->uspd[u] ;
            for ( j = 0 ; j < i  ; j++ )
            {
                jBasis = bases->entries[j] ;
                jData  = parameters->entries[j] ;
                j0     = Array1D_Item ( basisIndices, j ) ;
                nJ     = jData->norbitals ;
                xJ     = Coordinates3_RowPointer ( coordinates3, j ) ;
                s      = RealArray2D_AllocateWithExtents ( iBasis->nbasisw, jBasis->nbasisw, NULL ) ;
                GaussianBasisIntegrals_2Overlap ( iBasis, xI, jBasis, xJ, s ) ;
                GaussianBasis_TransformIntegrals2 ( &s, iBasis->c2o, jBasis->c2o ) ;
                if ( s != NULL )
                {
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        for ( v = 0 ; v < nJ ; v++ )
                        {
                            SymmetricMatrix_Item ( oneElectronMatrix, u+i0, v+j0 ) += ( 0.5e+00 * ( iData->beta[u] + jData->beta[v] ) * Array2D_Item ( s, u, v ) ) ;
                        }
                    }
                    RealArray2D_Deallocate ( &s ) ;
                }
            }
        }
    }
}
