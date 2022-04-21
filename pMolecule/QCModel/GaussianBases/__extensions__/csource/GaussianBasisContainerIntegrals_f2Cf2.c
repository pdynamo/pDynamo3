/*==================================================================================================================================
! . Container integrals - 4 basis, 2 electrons, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "GaussianBasisContainerIntegrals_f2Cf2.h"
# include "GaussianBasisIntegrals_f2Cf2.h"
# include "Integer.h"
# include "IntegerUtilities.h"
# include "RealUtilities.h"

/*

  Two-electron integrals (ij|kl) for a unique basis of size n:

  - number of different (ij) or (kl) pairs is p = 1/2*n*(n+1)
  - number of different integrals                 1/2*p*(p+1)

  Limits on basis loops (care required with corner cases when generalizing to atoms and shells):

  - i = 1 to n ; j = 1 to i ; k = 1 to i ; l = 1 to j if k = i or 1 to k otherwise.

*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void ProcessTEIs  ( const Integer          i0              ,
                           const Integer          j0              ,
                           const Integer          k0              ,
                           const Integer          l0              ,
                                 Block           *block           ,
                                 BlockStorage    *teis            ,
                                 Status          *status          ) ;
static void ProcessTEIsD ( const Boolean          doCoulomb       ,
                           const Boolean          doExchange      ,
                           const Integer          i               ,
                           const Integer          j               ,
                           const Integer          k               ,
                           const Integer          l               ,
                           const Integer          i0              ,
                           const Integer          j0              ,
                           const Integer          k0              ,
                           const Integer          l0              ,
                           const Real             exchangeScaling ,
                           const SymmetricMatrix *dTotal          ,
                           const SymmetricMatrix *dSpin           ,
                                 Block           *block           ,
                                 Coordinates3    *gradients3      ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the two-electron integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _TEIs_BlockSize 1024
# define _TEIs_UnderFlow 1.0e-12
void GaussianBasisContainerIntegrals_f2Cf2i ( const GaussianBasisContainer *self         ,
                                              const Coordinates3           *coordinates3 ,
                                                    BlockStorage           *teis         ,
                                                    Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( teis         != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Block         *block ;
        auto Integer        c, i, i0, *iWork, j, j0, k, k0, l, l0, n, s4 ;
        auto GaussianBasis *iBasis, *jBasis, *kBasis, *lBasis ;
        auto Real           d, *rI, rIJ[3], rIJ2, *rJ, *rK, rKL[3], rKL2, *rL, *rWork ;
        /* . Initialization. */
        BlockStorage_Empty ( teis ) ;
        teis->blockSize      = _TEIs_BlockSize ;
        teis->checkUnderFlow = True ;
        teis->nIndices16     = 4 ;
        teis->nReal          = 1 ;
        teis->underFlow      = _TEIs_UnderFlow ;
        n     = GaussianBasisContainer_LargestBasis ( self, False ) ;
        block = Block_Allocate ( n*n*n*n, 4, 0, 1, status ) ;
        n     = GaussianBasisContainer_LargestShell ( self, True ) ;
        s4    = n*n*n*n ;
        iWork = Integer_Allocate ( 3*s4, status ) ;
        rWork = Real_Allocate    ( 3*s4, status ) ;
        if ( ! Status_IsOK ( status ) ) goto FinishUp ;
        /* . Quadruple loop over centers. */
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
                for ( c = 0, rIJ2 = 0.0e+00 ; c < 3 ; c++ ) { d = rI[c] - rJ[c] ; rIJ[c] = d ; rIJ2 += d * d ; }
                for ( k = 0 ; k <= i ; k++ )
                {
                    kBasis = self->entries[k] ;
                    k0     = Array1D_Item ( self->centerFunctionPointers, k ) ;
                    rK     = Coordinates3_RowPointer ( coordinates3, k ) ;
                    for ( l = 0 ; l <= k ; l++ )
                    {
                        lBasis = self->entries[l] ;
                        l0     = Array1D_Item ( self->centerFunctionPointers, l ) ;
                        rL     = Coordinates3_RowPointer ( coordinates3, l ) ;
                        for ( c = 0, rKL2 = 0.0e+00 ; c < 3 ; c++ ) { d = rK[c] - rL[c] ; rKL[c] = d ; rKL2 += d * d ; }
/* . Need flag for j < l. */
                        GaussianBasisIntegrals_f2Cf2i ( iBasis, rI, jBasis, rJ, rIJ, rIJ2, kBasis, rK, lBasis, rL, rKL, rKL2, ( j < l ), s4, iWork, rWork, block ) ;
                        ProcessTEIs ( i0, j0, k0, l0, block, teis, status ) ;
                        if ( ! Status_IsOK ( status ) ) goto FinishUp ;
                    }
                }
            }
        }
FinishUp:
        if ( ! Status_IsOK ( status ) ) BlockStorage_Deallocate ( &teis ) ;
        Block_Deallocate   ( &block ) ;
        Integer_Deallocate ( &iWork ) ;
        Real_Deallocate    ( &rWork ) ;
    }
}
# undef _TEIs_BlockSize
# undef _TEIs_UnderFlow

/*----------------------------------------------------------------------------------------------------------------------------------
! . The two-electron integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f2Cf2R1 ( const GaussianBasisContainer *self            ,
                                               const Coordinates3           *coordinates3    ,
                                               const SymmetricMatrix        *dTotal          ,
                                               const SymmetricMatrix        *dSpin           ,
                                               const Boolean                 doCoulomb       ,
                                               const Real                    exchangeScaling ,
                                                     Coordinates3           *gradients3      ,
                                                     Status                 *status          )
{
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( dTotal       != NULL ) &&
         ( gradients3   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Block         *block ;
        auto Boolean        doExchange ;
        auto Integer        c, i, i0, *iWork, j, j0, k, k0, l, l0, n, s4 ;
        auto GaussianBasis *iBasis, *jBasis, *kBasis, *lBasis ;
        auto Real           d, *rI, rIJ[3], rIJ2, *rJ, *rK, rKL[3], rKL2, *rL, *rWork ;
        n     = GaussianBasisContainer_LargestBasis ( self, False ) ;
        block = Block_Allocate ( n*n*n*n, 4, 0, 9, status ) ;
        n     = GaussianBasisContainer_LargestShell ( self, True ) ;
        s4    = n*n*n*n ;
        iWork = Integer_Allocate (  6*s4, status ) ;
        rWork = Real_Allocate    ( 11*s4, status ) ;
        if ( ! Status_IsOK ( status ) ) goto FinishUp ;
        doExchange = ( exchangeScaling != 0.0e+00 ) ;
        /* . Quadruple loop over centers. */
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
                for ( c = 0, rIJ2 = 0.0e+00 ; c < 3 ; c++ ) { d = rI[c] - rJ[c] ; rIJ[c] = d ; rIJ2 += d * d ; }
                for ( k = 0 ; k <= i ; k++ )
                {
                    kBasis = self->entries[k] ;
                    k0     = Array1D_Item ( self->centerFunctionPointers, k ) ;
                    rK     = Coordinates3_RowPointer ( coordinates3, k ) ;
                    for ( l = 0 ; l <= k ; l++ )
                    {
                        lBasis = self->entries[l] ;
                        l0     = Array1D_Item ( self->centerFunctionPointers, l ) ;
                        rL     = Coordinates3_RowPointer ( coordinates3, l ) ;
                        for ( c = 0, rKL2 = 0.0e+00 ; c < 3 ; c++ ) { d = rK[c] - rL[c] ; rKL[c] = d ; rKL2 += d * d ; }
/* . Need flag for j < l. */
                        GaussianBasisIntegrals_f2Cf2r1 ( iBasis, rI, jBasis, rJ, rIJ, rIJ2, kBasis, rK, lBasis, rL, rKL, rKL2, ( j < l ), s4, iWork, rWork, block ) ;
                        ProcessTEIsD ( doCoulomb, doExchange, i, j, k, l, i0, j0, k0, l0, exchangeScaling, dTotal, dSpin, block, gradients3 ) ;
                    }
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
! . Process the TEIs.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void ProcessTEIs  ( const Integer        i0     ,
                           const Integer        j0     ,
                           const Integer        k0     ,
                           const Integer        l0     ,
                                 Block         *block  ,
                                 BlockStorage  *teis   ,
                                 Status        *status )
{
    if ( block->count > 0 )
    {
        auto Integer     i, m4 ;
        auto Cardinal16 *indices16 = block->indices16 ;
        auto Real       *integrals = block->data      ;
        for ( i = 0 ; i < block->count ; i++ )
        {
            m4 = 4 * i ;
            indices16[m4  ] += i0 ;
            indices16[m4+1] += j0 ;
            indices16[m4+2] += k0 ;
            indices16[m4+3] += l0 ;
        }
        BlockStorage_AddData ( teis, block->count, integrals, indices16, NULL, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Process the TEI derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define BFINDEX(i) ( i * ( i + 1 ) ) / 2
static void ProcessTEIsD ( const Boolean          doCoulomb       ,
                           const Boolean          doExchange      ,
                           const Integer          i               ,
                           const Integer          j               ,
                           const Integer          k               ,
                           const Integer          l               ,
                           const Integer          i0              ,
                           const Integer          j0              ,
                           const Integer          k0              ,
                           const Integer          l0              ,
                           const Real             exchangeScaling ,
                           const SymmetricMatrix *dTotal          ,
                           const SymmetricMatrix *dSpin           ,
                                 Block           *block           ,
                                 Coordinates3    *gradients3      )
{
    if ( block->count > 0 )
    {
        auto Boolean     doSpin = ( dSpin != NULL ) ;
        auto Integer     c, i1, i2, i3, i4, m4, m9, nIJ, nIK, nIL, nJK, nJL, nKL, t ;
        auto Real        d, dIx, dIy, dIz, dJx, dJy, dJz, dKx, dKy, dKz, scaling ;
        auto Cardinal16 *indices16 = block->indices16 ;
        auto Real       *integrals = block->data      ;
        dIx = dIy = dIz = dJx = dJy = dJz = dKx = dKy = dKz = 0.0e+00 ;
        for ( c = 0 ; c < block->count ; c++ )
        {
            m4      = 4 * c ;
            m9      = 9 * c ;
            i1      = indices16[m4  ] + i0 ;
            i2      = indices16[m4+1] + j0 ;
            i3      = indices16[m4+2] + k0 ;
            i4      = indices16[m4+3] + l0 ;
            scaling = 1.0e+00 ;
	    if ( i1 < i2 ) { t = i1 ; i1 = i2 ; i2 = t ; }
            if ( i3 < i4 ) { t = i3 ; i3 = i4 ; i4 = t ; }
            if ( ( i1 < i3 ) || ( ( i1 == i3 ) && ( i2 < i4  ) ) ) { t = i1 ; i1 = i3 ; i3 = t ; t = i2 ; i2 = i4 ; i4 = t ; }
	    if ( i1 == i2 ) scaling *= 0.5e+00 ;
	    if ( i3 == i4 ) scaling *= 0.5e+00 ;
            if ( ( i1 == i3 ) && ( i2 == i4 ) ) scaling *= 0.5e+00 ;
            /*
            In Fock building all off-diagonal elements are scaled by 1/2, and the energy is 1/2 tr ( P * F ).
            The latter is 1/2 * Pii*Fii + Pij*Fij (i>j) so off-diagonal elements occur twice which means the
            scaling is unnecessary. Likewise, the 1/2 removes an extra factor of 2 from the density factor.
            */
            /* . Coulomb. */
            if ( doCoulomb )
            {
                nIJ = BFINDEX ( i1 ) + i2 ;
                nKL = BFINDEX ( i3 ) + i4 ;
                d = 4.0e+00 * scaling * dTotal->data[nIJ] * dTotal->data[nKL] ;
            }
            else d = 0.0e+00 ;
            /* . Exchange. */
            if ( doExchange )
            {
                nIK = BFINDEX ( i1 ) + i3 ;
                nIL = BFINDEX ( i1 ) + i4 ;
                if ( i2 > i3 ) nJK = BFINDEX ( i2 ) + i3 ;
                else           nJK = BFINDEX ( i3 ) + i2 ;
                if ( i2 > i4 ) nJL = BFINDEX ( i2 ) + i4 ;
                else           nJL = BFINDEX ( i4 ) + i2 ;
                scaling *= exchangeScaling ;
                d -= scaling * ( dTotal->data[nIK] * dTotal->data[nJL] +
                                 dTotal->data[nIL] * dTotal->data[nJK] ) ;
                if ( doSpin ) d -= scaling * ( dSpin->data[nIK] * dSpin->data[nJL] +
                                               dSpin->data[nIL] * dSpin->data[nJK] ) ;
            }
            dIx += d * integrals[m9  ] ;
            dIy += d * integrals[m9+1] ;
            dIz += d * integrals[m9+2] ;
            dJx += d * integrals[m9+3] ;
            dJy += d * integrals[m9+4] ;
            dJz += d * integrals[m9+5] ;
            dKx += d * integrals[m9+6] ;
            dKy += d * integrals[m9+7] ;
            dKz += d * integrals[m9+8] ;
        }
        Coordinates3_IncrementRow ( gradients3, i, dIx, dIy, dIz ) ;
        Coordinates3_IncrementRow ( gradients3, j, dJx, dJy, dJz ) ;
        Coordinates3_IncrementRow ( gradients3, k, dKx, dKy, dKz ) ;
        Coordinates3_DecrementRow ( gradients3, l, dIx, dIy, dIz ) ;
        Coordinates3_DecrementRow ( gradients3, l, dJx, dJy, dJz ) ;
        Coordinates3_DecrementRow ( gradients3, l, dKx, dKy, dKz ) ;
    }
}
# undef BFINDEX
