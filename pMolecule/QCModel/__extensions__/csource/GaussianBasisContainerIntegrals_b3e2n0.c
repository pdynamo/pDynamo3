/*==================================================================================================================================
! . Container integrals - 3 basis, 2 electrons, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "BlockStorage.h"
# include "Integer.h"
# include "GaussianBasisContainerIntegrals_b3e2n0.h"
# include "GaussianBasisIntegrals_b3e2n0.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void ProcessFitIntegrals  ( const Integer          i0            ,
                                   const Integer          j0            ,
                                   const Integer          f0            ,
                                         Block           *block         ,
                                         BlockStorage    *fitIntegrals  ,
                                         Status          *status        ) ;
static void ProcessFitIntegralsD ( const Integer          i             ,
                                   const Integer          j             ,
                                   const Integer          f             ,
                                   const Integer          i0            ,
                                   const Integer          j0            ,
                                   const Integer          f0            ,
                                   const SymmetricMatrix *sDensity      ,
                                   const RealArray1D     *oDensity      ,
                                         Block           *block         ,
                                         Coordinates3    *gradients3    ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the two-electron integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _FitIntegrals_BlockSize 1024
# define _FitIntegrals_UnderFlow 1.0e-12
void GaussianBasisContainerIntegrals_ElectronFit ( const GaussianBasisContainer *self         ,
                                                   const IntegerArray1D         *selfIndices  ,
                                                   const GaussianBasisContainer *other        ,
                                                   const IntegerArray1D         *otherIndices ,
                                                   const Coordinates3           *coordinates3 ,
                                                         BlockStorage           *fitIntegrals ,
                                                         Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( selfIndices  != NULL ) &&
         ( other        != NULL ) &&
         ( otherIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( fitIntegrals != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Block         *block ;
        auto Integer        c, f, f0, i, i0, j, j0, m, n ;
        auto GaussianBasis *fBasis, *iBasis, *jBasis ;
        auto Real           d, *rF, *rI, rIJ[3], rIJ2, *rJ ;
        /* . Initialization. */
        BlockStorage_Empty ( fitIntegrals ) ;
        fitIntegrals->blockSize      = _FitIntegrals_BlockSize ;
        fitIntegrals->checkUnderFlow = True ;
        fitIntegrals->nIndices16     = 1 ;
        fitIntegrals->nIndices32     = 1 ;
        fitIntegrals->nReal          = 1 ;
        fitIntegrals->underFlow      = _FitIntegrals_UnderFlow ;
        m     = GaussianBasisContainer_LargestBasis ( self , True ) ;
        n     = GaussianBasisContainer_LargestBasis ( other, True ) ;
        block = Block_Allocate ( m*m*n, 3, 1, 1, status ) ;
        if ( ! Status_IsOK ( status ) ) goto FinishUp ;
        /* . Triple loop over centers. */
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            iBasis = self->entries[i] ;
            i0     = Array1D_Item ( selfIndices, i ) ;
            rI     = Coordinates3_RowPointer ( coordinates3, i ) ;
            for ( j = 0 ; j <= i ; j++ )
            {
                jBasis = self->entries[j] ;
                j0     = Array1D_Item ( selfIndices, j ) ;
                rJ     = Coordinates3_RowPointer ( coordinates3, j ) ;
                for ( c = 0, rIJ2 = 0.0e+00 ; c < 3 ; c++ ) { d = rI[c] - rJ[c] ; rIJ[c] = d ; rIJ2 += d * d ; }
                for ( f = 0 ; f < other->capacity ; f++ )
                {
                    fBasis = other->entries[f] ;
                    f0     = Array1D_Item ( otherIndices, f ) ;
                    rF     = Coordinates3_RowPointer ( coordinates3, f ) ;
                    GaussianBasisIntegrals_ElectronFit ( iBasis, rI, jBasis, rJ, rIJ, rIJ2, fBasis, rF, block ) ;
                    ProcessFitIntegrals ( i0, j0, f0, block, fitIntegrals, status ) ;
                    if ( ! Status_IsOK ( status ) ) goto FinishUp ;
                }
            }
        }
FinishUp:
        if ( ! Status_IsOK ( status ) ) BlockStorage_Empty ( fitIntegrals ) ;
        Block_Deallocate ( &block ) ;
    }
}
# undef _FitIntegrals_BlockSize
# undef _FitIntegrals_UnderFlow

/*----------------------------------------------------------------------------------------------------------------------------------
! . The two-electron integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . oDensity = ( fPotential + wVector ). */
void GaussianBasisContainerIntegrals_ElectronFitD ( const GaussianBasisContainer *self         ,
                                                    const IntegerArray1D         *selfIndices  ,
                                                    const GaussianBasisContainer *other        ,
                                                    const IntegerArray1D         *otherIndices ,
                                                    const Coordinates3           *coordinates3 ,
                                                    const SymmetricMatrix        *sDensity     ,
                                                    const RealArray1D            *oDensity     ,
                                                          Coordinates3           *gradients3   ,
                                                          Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( selfIndices  != NULL ) &&
         ( other        != NULL ) &&
         ( otherIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( sDensity     != NULL ) &&
         ( oDensity     != NULL ) &&
         ( gradients3   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Block         *block ;
        auto Integer        c, f, f0, i, i0, j, j0, m, n ;
        auto GaussianBasis *fBasis, *iBasis, *jBasis ;
        auto Real           d, *rF, *rI, rIJ[3], rIJ2, *rJ ;
        m     = GaussianBasisContainer_LargestBasis ( self , True ) ;
        n     = GaussianBasisContainer_LargestBasis ( other, True ) ;
        block = Block_Allocate ( m*m*n, 3, 0, 6, status ) ;
        if ( ! Status_IsOK ( status ) ) goto FinishUp ;
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            iBasis = self->entries[i] ;
            i0     = Array1D_Item ( selfIndices, i ) ;
            rI     = Coordinates3_RowPointer ( coordinates3, i ) ;
            for ( j = 0 ; j <= i ; j++ )
            {
                jBasis = self->entries[j] ;
                j0     = Array1D_Item ( selfIndices, j ) ;
                rJ     = Coordinates3_RowPointer ( coordinates3, j ) ;
                for ( c = 0, rIJ2 = 0.0e+00 ; c < 3 ; c++ ) { d = rI[c] - rJ[c] ; rIJ[c] = d ; rIJ2 += d * d ; }
                for ( f = 0 ; f < other->capacity ; f++ )
                {
                    fBasis = other->entries[f] ;
                    f0     = Array1D_Item ( otherIndices, f ) ;
                    rF     = Coordinates3_RowPointer ( coordinates3, f ) ;
                    if ( ( rI == rJ ) && ( rI == rF ) ) continue ; /* . Triple diagonal terms are zero. */
                    GaussianBasisIntegrals_ElectronFitD ( iBasis, rI, jBasis, rJ, rIJ, rIJ2, fBasis, rF, block ) ;
                    ProcessFitIntegralsD ( i, j, f, i0, j0, f0, sDensity, oDensity, block, gradients3 ) ;
                }
            }
        }
FinishUp:
        Block_Deallocate ( &block ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Process the FitIntegrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define BFINDEX(i) ( i * ( i + 1 ) ) / 2
static void ProcessFitIntegrals ( const Integer       i0           ,
                                  const Integer       j0           ,
                                  const Integer       f0           ,
                                        Block        *block        ,
                                        BlockStorage *fitIntegrals ,
                                        Status       *status       )
{
    if ( block->count > 0 )
    {
        auto Integer     c, f, i, ij, j, m3 ;
        auto Cardinal16 *indices16 = block->indices16 ;
        auto Cardinal32 *indices32 = block->indices32 ;
        auto Real       *integrals = block->data      ;
        for ( c = 0 ; c < block->count ; c++ )
        {
            m3 = 3 * c ;
            i  = indices16[m3  ] + i0 ;
            j  = indices16[m3+1] + j0 ;
            f  = indices16[m3+2] + f0 ;
            if ( i >= j ) ij = BFINDEX ( i ) + j ;
            else          ij = BFINDEX ( j ) + i ;
            indices16[c] = f  ;
            indices32[c] = ij ;
        }
        BlockStorage_AddData ( fitIntegrals, block->count, integrals, indices16, indices32, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Process the FitIntegral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void ProcessFitIntegralsD ( const Integer          i          ,
                                   const Integer          j          ,
                                   const Integer          f          ,
                                   const Integer          i0         ,
                                   const Integer          j0         ,
                                   const Integer          f0         ,
                                   const SymmetricMatrix *sDensity   ,
                                   const RealArray1D     *oDensity   ,
                                         Block           *block      ,
                                         Coordinates3    *gradients3 )
{
    if ( block->count > 0 )
    {
        auto Integer     c, ff, i1, i2, m3, m6, t ;
        auto Real        d, dIx, dIy, dIz, dJx, dJy, dJz ;
        auto Cardinal16 *indices16 = block->indices16 ;
        auto Real       *integrals = block->data      ;
        dIx = dIy = dIz = dJx = dJy = dJz = 0.0e+00 ;
        for ( c = 0 ; c < block->count ; c++ )
        {
            m3 = 3 * c ;
            m6 = 6 * c ;
            i1 = indices16[m3  ] + i0 ;
            i2 = indices16[m3+1] + j0 ;
            ff = indices16[m3+2] + f0 ;
            if ( i1 < i2 ) { t = i1 ; i1 = i2 ; i2 = t ; }
            d    = SymmetricMatrix_Item ( sDensity, i1, i2 ) * Array1D_Item ( oDensity, ff ) ;
            dIx += d * integrals[m6  ] ;
            dIy += d * integrals[m6+1] ;
            dIz += d * integrals[m6+2] ;
            dJx += d * integrals[m6+3] ;
            dJy += d * integrals[m6+4] ;
            dJz += d * integrals[m6+5] ;
        }
        Coordinates3_IncrementRow ( gradients3, i, dIx, dIy, dIz ) ;
        Coordinates3_IncrementRow ( gradients3, j, dJx, dJy, dJz ) ;
        Coordinates3_DecrementRow ( gradients3, f, dIx, dIy, dIz ) ;
        Coordinates3_DecrementRow ( gradients3, f, dJx, dJy, dJz ) ;
    }
}
# undef BFINDEX
