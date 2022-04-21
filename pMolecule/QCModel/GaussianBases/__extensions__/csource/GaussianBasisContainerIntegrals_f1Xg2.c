/*==================================================================================================================================
! . Container integrals - 3 basis, 2 electrons, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "BlockStorage.h"
# include "Integer.h"
# include "IntegerUtilities.h"
# include "GaussianBasisContainerIntegrals_f1Xg2.h"
# include "GaussianBasisIntegrals_f1Xg2.h"
# include "NumericalMacros.h"
# include "RealUtilities.h"

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
                                   const SymmetricMatrix *density       ,
                                   const RealArray1D     *xVector       ,
                                         Block           *block         ,
                                         Coordinates3    *gradients3    ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the two-electron integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _FitIntegrals_BlockSize 1024
# define _FitIntegrals_UnderFlow 1.0e-12
void GaussianBasisContainerIntegrals_f1Xg2i ( const GaussianBasisContainer *self         ,
                                              const GaussianBasisContainer *other        ,
                                              const Coordinates3           *coordinates3 ,
                                              const GaussianBasisOperator   operator     ,
                                                    BlockStorage           *fitIntegrals ,
                                                    Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( other        != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( fitIntegrals != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Block         *block ;
        auto Integer        c, f, f0, i, i0, *iWork, j, j0, m, n, s3 ;
        auto GaussianBasis *fBasis, *iBasis, *jBasis ;
        auto Real           d, *rF, *rI, rIJ[3], rIJ2 = 0.0e+00, *rJ, *rWork ;
        /* . Initialization. */
        BlockStorage_Empty ( fitIntegrals ) ;
        fitIntegrals->blockSize      = _FitIntegrals_BlockSize ;
        fitIntegrals->checkUnderFlow = True ;
        fitIntegrals->nIndices16     = 1 ;
        fitIntegrals->nIndices32     = 1 ;
        fitIntegrals->nReal          = 1 ;
        fitIntegrals->underFlow      = _FitIntegrals_UnderFlow ;
        m     = GaussianBasisContainer_LargestBasis ( self , False ) ;
        n     = GaussianBasisContainer_LargestBasis ( other, False ) ;
        block = Block_Allocate ( m*m*n, 3, 1, 1, status ) ;
        m     = GaussianBasisContainer_LargestShell ( self , True ) ;
        n     = GaussianBasisContainer_LargestShell ( other, True ) ;
        s3    = m*m*n ;
        if ( operator == GaussianBasisOperator_Overlap )
        {
            iWork = NULL ;
            rWork = Real_Allocate    ( 2*s3, status ) ;
        }
        else
        {
            iWork = Integer_Allocate ( 3*s3, status ) ;
            rWork = Real_Allocate    ( 3*s3, status ) ;
        }
        if ( ! Status_IsOK ( status ) ) goto FinishUp ;
        /* . Triple loop over centers. */
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            iBasis = self->entries[i] ;
            i0     = Array1D_Item ( self->centerFunctionPointers, i ) ;
            rI     = Coordinates3_RowPointer ( coordinates3, i ) ;
            for ( j = 0 ; j <= i ; j++ )
            {
                jBasis = self->entries[j] ;
                j0     = Array1D_Item ( self->centerFunctionPointers, j ) ;
                rJ     = Coordinates3_RowPointer ( coordinates3, j ) ;
                if ( operator != GaussianBasisOperator_Overlap ) 
                {
                    for ( c = 0, rIJ2 = 0.0e+00 ; c < 3 ; c++ ) { d = rI[c] - rJ[c] ; rIJ[c] = d ; rIJ2 += d * d ; }
                }
                for ( f = 0 ; f < other->capacity ; f++ )
                {
                    fBasis = other->entries[f] ;
                    f0     = Array1D_Item ( other->centerFunctionPointers, f ) ;
                    rF     = Coordinates3_RowPointer ( coordinates3, f ) ;
                         if ( operator == GaussianBasisOperator_AntiCoulomb ) GaussianBasisIntegrals_f1Ag2i ( iBasis, rI, jBasis, rJ, rIJ, rIJ2, fBasis, rF, s3, iWork, rWork, block ) ;
                    else if ( operator == GaussianBasisOperator_Coulomb     ) GaussianBasisIntegrals_f1Cg2i ( iBasis, rI, jBasis, rJ, rIJ, rIJ2, fBasis, rF, s3, iWork, rWork, block ) ;
                    else if ( operator == GaussianBasisOperator_Overlap     ) GaussianBasisIntegrals_f1Og2i ( iBasis, rI, jBasis, rJ,            fBasis, rF, s3,        rWork, block ) ;
/*printf ( "\nProcessing Fit Integrals: %u %d %d %d %d | %u %d %d %d %d | %u %d %d %d %d\n",
                       iBasis->isSpherical, i, i0, iBasis->nBasis, iBasis->nCBF ,
                       jBasis->isSpherical, j, j0, jBasis->nBasis, jBasis->nCBF ,
                       fBasis->isSpherical, f, f0, fBasis->nBasis, fBasis->nCBF ) ; */
                    ProcessFitIntegrals ( i0, j0, f0, block, fitIntegrals, status ) ;
                    if ( ! Status_IsOK ( status ) ) goto FinishUp ;
                }
            }
        }
FinishUp:
        if ( ! Status_IsOK ( status ) ) BlockStorage_Empty ( fitIntegrals ) ;
        Block_Deallocate   ( &block ) ;
        Integer_Deallocate ( &iWork ) ;
        Real_Deallocate    ( &rWork ) ;
    }
}
# undef _FitIntegrals_BlockSize
# undef _FitIntegrals_UnderFlow

/*----------------------------------------------------------------------------------------------------------------------------------
! . The two-electron integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . density is orbital density matrix, xVector is one of A, D, W, A+D or A+W. */
void GaussianBasisContainerIntegrals_f1Xg2R1 ( const GaussianBasisContainer *self         ,
                                               const GaussianBasisContainer *other        ,
                                               const Coordinates3           *coordinates3 ,
                                               const SymmetricMatrix        *density      ,
                                               const RealArray1D            *xVector      ,
                                               const GaussianBasisOperator   operator     ,
                                                     Coordinates3           *gradients3   ,
                                                     Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( other        != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( density      != NULL ) &&
         ( xVector      != NULL ) &&
         ( gradients3   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Block         *block ;
        auto Integer        c, f, f0, i, i0, *iWork, j, j0, m, n, s3 ;
        auto GaussianBasis *fBasis, *iBasis, *jBasis ;
        auto Real           d, *rF, *rI, rIJ[3], rIJ2 = 0.0e+00, *rJ, *rWork ;
        m     = GaussianBasisContainer_LargestBasis ( self , False ) ;
        n     = GaussianBasisContainer_LargestBasis ( other, False ) ;
        block = Block_Allocate ( m*m*n, 3, 0, 6, status ) ;
        m     = GaussianBasisContainer_LargestShell ( self , True ) ;
        n     = GaussianBasisContainer_LargestShell ( other, True ) ;
        s3    = m*m*n ;
        if ( operator == GaussianBasisOperator_Overlap )
        {
            iWork = NULL ;
            rWork = Real_Allocate    ( 7*s3, status ) ;
        }
        else
        {
            iWork = Integer_Allocate ( 6*s3, status ) ;
            rWork = Real_Allocate    ( 8*s3, status ) ;
        }
        if ( ! Status_IsOK ( status ) ) goto FinishUp ;
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            iBasis = self->entries[i] ;
            i0     = Array1D_Item ( self->centerFunctionPointers, i ) ;
            rI     = Coordinates3_RowPointer ( coordinates3, i ) ;
            for ( j = 0 ; j <= i ; j++ )
            {
                jBasis = self->entries[j] ;
                j0     = Array1D_Item ( self->centerFunctionPointers, j ) ;
                rJ     = Coordinates3_RowPointer ( coordinates3, j ) ;
                if ( operator != GaussianBasisOperator_Overlap ) 
                {
                    for ( c = 0, rIJ2 = 0.0e+00 ; c < 3 ; c++ ) { d = rI[c] - rJ[c] ; rIJ[c] = d ; rIJ2 += d * d ; }
                }
                for ( f = 0 ; f < other->capacity ; f++ )
                {
                    fBasis = other->entries[f] ;
                    f0     = Array1D_Item ( other->centerFunctionPointers, f ) ;
                    rF     = Coordinates3_RowPointer ( coordinates3, f ) ;
                    if ( ( rI == rJ ) && ( rI == rF ) ) continue ; /* . Triple diagonal terms are zero. */
                         if ( operator == GaussianBasisOperator_AntiCoulomb ) GaussianBasisIntegrals_f1Ag2r1 ( iBasis, rI, jBasis, rJ, rIJ, rIJ2, fBasis, rF, s3, iWork, rWork, block ) ;
                    else if ( operator == GaussianBasisOperator_Coulomb     ) GaussianBasisIntegrals_f1Cg2r1 ( iBasis, rI, jBasis, rJ, rIJ, rIJ2, fBasis, rF, s3, iWork, rWork, block ) ;
                    else if ( operator == GaussianBasisOperator_Overlap     ) GaussianBasisIntegrals_f1Og2r1 ( iBasis, rI, jBasis, rJ,            fBasis, rF, s3,        rWork, block ) ;
                    ProcessFitIntegralsD ( i, j, f, i0, j0, f0, density, xVector, block, gradients3 ) ;
                }
            }
        }
FinishUp:
        Block_Deallocate   ( &block ) ;
        Integer_Deallocate ( &iWork ) ;
        Real_Deallocate    ( &rWork ) ;
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
                                   const SymmetricMatrix *density    ,
                                   const RealArray1D     *xVector    ,
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
            d    = SymmetricMatrix_Item ( density, i1, i2 ) * Array1D_Item ( xVector, ff ) ;
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
