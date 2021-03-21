/*==================================================================================================================================
! . MNDO electron-nuclear and electron-electron interactions (OEI and TEIs, respectively).
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "Integer.h"
# include "MNDOElectronNuclearTEIs.h"
# include "MNDOIntegrals.h"
# include "MNDOParameters.h"
# include "Real.h"
# include "RealArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GetGradientDensityTerms   ( const MNDOParameters        *iData     ,  
                                        const Integer                i0        ,  
                                        const MNDOParameters        *jData     ,  
                                        const Integer                j0        ,  
                                        const SymmetricMatrix       *dTotal    ,  
                                        const SymmetricMatrix       *dSpin     ,  
                                              RealArray1D           *dOneI     ,  
                                              RealArray1D           *dOneJ     ,  
                                              RealArray2D           *dTwoIJ    ) ;
static void GetGradientDensityTermsCI ( const Integer                nActive   ,
                                        const Integer                nCore ,
                                        const MNDOParameters        *iData     ,
                                        const Integer                i0        ,
                                        const MNDOParameters        *jData     ,
                                        const Integer                j0        ,
                                        const DoubleSymmetricMatrix *twoPDM    ,
                                        const RealArray2D           *orbitals  ,
                                        const SymmetricMatrix       *dCore     ,
                                        const SymmetricMatrix       *dHF       ,
                                        const SymmetricMatrix       *dTotal    ,
                                        const SymmetricMatrix       *onePDM    ,
                                        const SymmetricMatrix       *zMatrix   ,
                                              RealArray1D           *tPDM1     ,
                                              RealArray2D           *tPDM2     ,
                                              RealArrayND           *tPDM3     ,
                                              RealArray1D           *dI        ,
                                              RealArray1D           *dJ        ,
                                              RealArray2D           *dIJ       ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . The electron-nuclear and electron-electron interaction gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDO_ElectronNuclearTEIGradients ( const MNDOParametersContainer *parameters   ,
                                        const IntegerArray1D          *basisIndices ,
                                        const Coordinates3            *coordinates3 ,
                                        const SymmetricMatrix         *dTotal       ,
                                        const SymmetricMatrix         *dSpin        ,
                                              Coordinates3            *gradients3   )
{
    if ( ( parameters   != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( dTotal       != NULL ) &&
         ( gradients3   != NULL ) )
    {
        auto Integer         i, i0, j, j0, nI, nJ ;
        auto MNDOParameters *iData , *jData  ;
        auto Real            gX, gY, gZ, *xI, *xJ ;
        auto RealArray1D    *dOneI, *dOneJ ;
        auto RealArray2D    *dTwoIJ ;
        for ( i = 0 ; i < Coordinates3_Rows ( coordinates3 ) ; i++ )
        {
            iData = parameters->entries[i] ;
            i0    = Array1D_Item ( basisIndices, i ) ;
            nI    = ( iData->norbitals * ( iData->norbitals + 1 ) ) / 2 ;
            xI    = Coordinates3_RowPointer ( coordinates3, i ) ;
            dOneI = RealArray1D_AllocateWithExtent ( nI, NULL ) ; RealArray1D_Set ( dOneI, 0.0e+00 ) ;
            for ( j = 0 ; j < i ; j++ )
            {
                jData  = parameters->entries[j] ;
                j0     = Array1D_Item ( basisIndices, j ) ;
                nJ     = ( jData->norbitals * ( jData->norbitals + 1 ) ) / 2 ;
                xJ     = Coordinates3_RowPointer ( coordinates3, j ) ;
                dOneJ  = RealArray1D_AllocateWithExtent  ( nJ    , NULL ) ; RealArray1D_Set ( dOneJ , 0.0e+00 ) ;
                dTwoIJ = RealArray2D_AllocateWithExtents ( nI, nJ, NULL ) ; RealArray2D_Set ( dTwoIJ, 0.0e+00 ) ;
                GetGradientDensityTerms ( iData, i0, jData, j0, dTotal, dSpin, dOneI, dOneJ, dTwoIJ ) ;
                MNDOIntegrals_MolecularFrame2CIntegralsD ( iData, i0, xI, jData, j0, xJ, dOneI, dOneJ, dTwoIJ, &gX, &gY, &gZ ) ;
                RealArray1D_Deallocate ( &dOneJ  ) ;
                RealArray2D_Deallocate ( &dTwoIJ ) ;
                Coordinates3_IncrementRow ( gradients3, i, gX, gY, gZ ) ;
                Coordinates3_DecrementRow ( gradients3, j, gX, gY, gZ ) ;
            }
            RealArray1D_Deallocate ( &dOneI ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The electron-nuclear and electron-electron interaction gradients for a CI calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDO_ElectronNuclearTEIGradientsCI ( const Integer                  nActive      ,
                                          const Integer                  nCore        ,
                                          const Integer                  nOrbitals    ,
                                          const MNDOParametersContainer *parameters   ,
                                          const IntegerArray1D          *basisIndices ,
                                          const Coordinates3            *coordinates3 ,
                                          const DoubleSymmetricMatrix   *twoPDM       ,
                                          const RealArray2D             *orbitals     , /* . Full set. */
                                          const SymmetricMatrix         *dCore        ,
                                          const SymmetricMatrix         *dHF          , /* . From HF calculation. */
                                          const SymmetricMatrix         *dTotalZ      , /* . From CI and CPHF calculation. */
                                          const SymmetricMatrix         *onePDM       ,
                                          const SymmetricMatrix         *zMatrix      ,
                                                Coordinates3            *gradients3   )
{
    if ( ( parameters   != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( twoPDM       != NULL ) &&
         ( orbitals     != NULL ) &&
         ( dCore        != NULL ) && /* ??? */
         ( dHF          != NULL ) &&
         ( dTotalZ      != NULL ) &&
         ( onePDM       != NULL ) &&
         ( zMatrix      != NULL ) &&
         ( gradients3   != NULL ) )
    {
        auto Integer      extents[3] ;
        auto RealArray1D *tPDM1 ;
        auto RealArray2D *tPDM2 ;
        auto RealArrayND *tPDM3 ;
        extents[0] = nActive ; extents[1] = nActive ; extents[2] = nActive ;
        tPDM1 = RealArray1D_AllocateWithExtent          ( nActive  ,          NULL ) ;
        tPDM2 = RealArray2D_AllocateWithExtents          ( nOrbitals, nActive, NULL ) ;
        tPDM3 = RealArrayND_AllocateWithShape ( 3        , extents, NULL ) ;
        if ( ( tPDM1 != NULL ) && ( tPDM2 != NULL ) && ( tPDM3 != NULL ) )
        {
            auto Integer         i, i0, j, j0, nI, nJ ;
            auto MNDOParameters *iData , *jData  ;
            auto Real            gX, gY, gZ, *xI, *xJ ;
            auto RealArray1D    *dOneI, *dOneJ ;
            auto RealArray2D    *dTwoIJ ;
            for ( i = 0 ; i < Coordinates3_Rows ( coordinates3 ) ; i++ )
            {
                iData = parameters->entries[i] ;
                i0    = Array1D_Item ( basisIndices, i ) ;
                nI    = ( iData->norbitals * ( iData->norbitals + 1 ) ) / 2 ;
                xI    = Coordinates3_RowPointer ( coordinates3, i ) ;
                dOneI = RealArray1D_AllocateWithExtent ( nI, NULL ) ; RealArray1D_Set ( dOneI, 0.0e+00 ) ;
                for ( j = 0 ; j < i ; j++ )
                {
                    jData  = parameters->entries[j] ;
                    j0     = Array1D_Item ( basisIndices, j ) ;
                    nJ     = ( jData->norbitals * ( jData->norbitals + 1 ) ) / 2 ;
                    xJ     = Coordinates3_RowPointer ( coordinates3, j ) ;
                    dOneJ  = RealArray1D_AllocateWithExtent ( nJ    , NULL ) ; RealArray1D_Set ( dOneJ , 0.0e+00 ) ;
                    dTwoIJ = RealArray2D_AllocateWithExtents ( nI, nJ, NULL ) ; RealArray2D_Set ( dTwoIJ, 0.0e+00 ) ;
                    GetGradientDensityTermsCI ( nActive  ,
                                                nCore    ,
                                                iData    ,
                                                i0       ,
                                                jData    ,
                                                j0       ,
                                                twoPDM   ,
                                                orbitals ,
                                                dCore    ,
                                                dHF      ,
                                                dTotalZ  ,
                                                onePDM   ,
                                                zMatrix  ,
                                                tPDM1    ,
                                                tPDM2    ,
                                                tPDM3    ,
                                                dOneI    ,
                                                dOneJ    ,
                                                dTwoIJ   ) ;
                    MNDOIntegrals_MolecularFrame2CIntegralsD ( iData, i0, xI, jData, j0, xJ, dOneI, dOneJ, dTwoIJ, &gX, &gY, &gZ ) ;
                    RealArray1D_Deallocate ( &dOneJ  ) ;
                    RealArray2D_Deallocate ( &dTwoIJ ) ;
                    Coordinates3_IncrementRow ( gradients3, i, gX, gY, gZ ) ;
                    Coordinates3_DecrementRow ( gradients3, j, gX, gY, gZ ) ;
                }
                RealArray1D_Deallocate ( &dOneI ) ;
            }
        }
        RealArray1D_Deallocate ( &tPDM1 ) ;
        RealArray2D_Deallocate ( &tPDM2 ) ;
        RealArrayND_Deallocate ( &tPDM3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The electron-nuclear and electron-electron interaction integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define MNDO_BLOCKSIZE 1024
# define MNDO_UNDERFLOW 1.0e-12
void MNDO_ElectronNuclearTEIIntegrals ( const MNDOParametersContainer *parameters           ,
                                        const IntegerArray1D          *basisIndices         ,
                                        const Coordinates3            *coordinates3         ,
                                              SymmetricMatrix         *oneElectronMatrix    ,
                                              BlockStorage           **twoElectronIntegrals )
{
    if ( ( parameters        != NULL ) &&
         ( basisIndices      != NULL ) &&
         ( coordinates3      != NULL ) &&
         ( oneElectronMatrix != NULL ) )
    {
        auto BlockStorage   *teis = NULL ;
        auto Integer         i, i0, j, j0, nI, nJ, u, v, w ;
        auto MNDOParameters *iData , *jData  ;
        auto Real           *xI, *xJ ;
        auto RealArray1D    *e1b, *e2a ;
        teis = BlockStorage_Allocate ( NULL ) ;
        teis->blockSize      = MNDO_BLOCKSIZE ;
        teis->checkUnderFlow = True ;
        teis->nIndices16     = 4 ;
        teis->nReal          = 1 ;
        teis->underFlow      = MNDO_UNDERFLOW ;
        for ( i = 0 ; i < Coordinates3_Rows ( coordinates3 ) ; i++ )
        {
            iData = parameters->entries[i] ;
            i0    = Array1D_Item ( basisIndices, i ) ;
            nI    = iData->norbitals ;
            xI    = Coordinates3_RowPointer ( coordinates3, i ) ;
            MNDOIntegrals_AddInOneCenterTEIs ( iData, i0, teis ) ;
            for ( j = 0 ; j < i ; j++ )
            {
                jData = parameters->entries[j] ;
                j0    = Array1D_Item ( basisIndices, j ) ;
                nJ    = jData->norbitals ;
                xJ    = Coordinates3_RowPointer ( coordinates3, j ) ;
                e1b   = RealArray1D_AllocateWithExtent ( ( nI * ( nI + 1 ) ) / 2, NULL ) ;
                e2a   = RealArray1D_AllocateWithExtent ( ( nJ * ( nJ + 1 ) ) / 2, NULL ) ;
                MNDOIntegrals_MolecularFrame2CIntegrals ( iData, i0, xI, jData, j0, xJ, e1b, e2a, teis ) ;
                for ( u = i0, w = 0 ; u < ( i0+nI ) ; u++ )
                {
                    for ( v = i0 ; v <= u ; v++, w++ ) SymmetricMatrix_Item ( oneElectronMatrix, u, v ) += e1b->data[w] ;
                }
                for ( u = j0, w = 0 ; u < ( j0+nJ ) ; u++ )
                {
                    for ( v = j0 ; v <= u ; v++, w++ ) SymmetricMatrix_Item ( oneElectronMatrix, u, v ) += e2a->data[w] ;
                }
                RealArray1D_Deallocate ( &e1b ) ;
                RealArray1D_Deallocate ( &e2a ) ;
            }
        }
        (*twoElectronIntegrals) = teis ;
    }
}
# undef MNDO_BLOCKSIZE
# undef MNDO_UNDERFLOW

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the densities for the gradient terms between atoms i and j.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GetGradientDensityTerms ( const MNDOParameters  *iData  ,
                                      const Integer          i0     ,
                                      const MNDOParameters  *jData  ,
                                      const Integer          j0     ,
                                      const SymmetricMatrix *dTotal ,
                                      const SymmetricMatrix *dSpin  ,
                                            RealArray1D     *dI     ,
                                            RealArray1D     *dJ     ,
                                            RealArray2D     *dIJ    )
{
    Integer      i, ij, j, k, kk, kl, l, ll, m, mk, ml, mn, n, nk, nl ;
    Real         aa, bb, f ;
    if ( ( dI != NULL ) && ( dJ != NULL ) && ( dIJ != NULL ) )
    {
        /* . One-center terms. */
        for ( i = ij = 0 ; i < iData->norbitals ; i++, ij++ )
        {
            mn = ( ( i + i0 ) * ( i + i0 + 1 ) ) / 2 + i0 ;
            for ( j = 0 ; j < i ; ij++, j++, mn++ ) Array1D_Item ( dI , ij ) = 2.0e+00 * dTotal->data[mn] ;
            Array1D_Item ( dI , ij ) = dTotal->data[mn] ;
        }
        for ( i = ij = 0 ; i < jData->norbitals ; i++, ij++ )
        {
            mn = ( ( i + j0 ) * ( i + j0 + 1 ) ) / 2 + j0 ;
            for ( j = 0 ; j < i ; ij++, j++, mn++ ) Array1D_Item ( dJ, ij ) = 2.0e+00 * dTotal->data[mn] ;
            Array1D_Item ( dJ, ij ) = dTotal->data[mn] ;
        }
        /* . Two-center exchange. */
        for ( k = i0 ; k < i0 + iData->norbitals ; ++k )
        {
	    aa = 1.0e+00 ;
	    kk = k * ( k + 1 ) / 2 ;
	    for ( l = k ; l < i0 + iData->norbitals ; ++l )
            {
	        ll = l * ( l + 1 ) / 2;
	        kl = ( k - i0 ) + ( ( l - i0 ) * ( l - i0 + 1 ) ) / 2 ;
	        for ( m = j0 ; m < j0 + jData->norbitals ; ++m )
                {
		    bb = 1.0e+00 ;
		    for ( n = m ; n < j0 + jData->norbitals ; ++n )
                    {
		        mn = ( m - j0 ) + ( ( n - j0 ) * ( n - j0 + 1 ) ) / 2 ;
		        mk = m + kk ;
		        nk = n + kk ;
		        ml = m + ll ;
		        nl = n + ll ;
                        f  = dTotal->data[mk] * dTotal->data[nl] + dTotal->data[nk] * dTotal->data[ml] ;
                        if ( dSpin != NULL ) f += ( dSpin->data[mk] * dSpin->data[nl] + dSpin->data[nk] * dSpin->data[ml] ) ;
                        Array2D_Item ( dIJ, kl, mn ) = 0.25e+00 * aa * bb * f ;
		        bb = 2.0e+00 ;
		    }
	        }
	        aa = 2.0e+00 ;
	    }
        }
        /* . Two-center Coulomb. */
        for ( i = 0 ; i < View1D_Extent ( dI ) ; i++ )
        {
            f = Array1D_Item ( dI, i ) ;
            for ( j = 0 ; j < View1D_Extent ( dJ ) ; j++ ) Array2D_Item ( dIJ, i, j ) -= ( f * Array1D_Item ( dJ, j ) ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the densities for the gradient terms between atoms i and j for a CI calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GetGradientDensityTermsCI ( const Integer                nActive  ,
                                        const Integer                nCore    ,
                                        const MNDOParameters        *iData    ,
                                        const Integer                i0       ,
                                        const MNDOParameters        *jData    ,
                                        const Integer                j0       ,
                                        const DoubleSymmetricMatrix *twoPDM   ,
                                        const RealArray2D           *orbitals ,
                                        const SymmetricMatrix       *dCore    ,
                                        const SymmetricMatrix       *dHF      ,
                                        const SymmetricMatrix       *dTotalZ  ,
                                        const SymmetricMatrix       *onePDM   ,
                                        const SymmetricMatrix       *zMatrix  ,
                                              RealArray1D           *tPDM1    ,
                                              RealArray2D           *tPDM2    ,
                                              RealArrayND           *tPDM3    ,
                                              RealArray1D           *dI       ,
                                              RealArray1D           *dJ       ,
                                              RealArray2D           *dIJ      )
{
    if ( ( dI != NULL ) && ( dJ != NULL ) && ( dIJ != NULL ) )
    {
        auto Integer  i, ij, j, k, kk, kl, kl0, l, ll, m, mk, ml, mn, mn0, n, nk, nl, p, q, r, s ;
        auto Real     aa, bb, f, f1, f2, f3, f4 ;
        /* . One-center terms. */
        for ( i = ij = 0 ; i < iData->norbitals ; i++, ij++ )
        {
            mn = ( ( i + i0 ) * ( i + i0 + 1 ) ) / 2 + i0 ;
            for ( j = 0 ; j < i ; ij++, j++, mn++ ) Array1D_Item ( dI, ij ) = 2.0e+00 * dTotalZ->data[mn] ;
            Array1D_Item ( dI, ij ) = dTotalZ->data[mn] ;
        }
        for ( i = ij = 0 ; i < jData->norbitals ; i++, ij++ )
        {
            mn = ( ( i + j0 ) * ( i + j0 + 1 ) ) / 2 + j0 ;
            for ( j = 0 ; j < i ; ij++, j++, mn++ ) Array1D_Item ( dJ, ij ) = 2.0e+00 * dTotalZ->data[mn] ;
            Array1D_Item ( dJ, ij ) = dTotalZ->data[mn] ;
        }
        /* . Two-center terms. */
        for ( k = i0 ; k < i0 + iData->norbitals ; ++k )
        {
	    aa = 1.0e+00 ;
	    kk = k * ( k + 1 ) / 2 ;
            /* . First orbital transformation. */
            for ( s = 0 ; s < nActive ; s++ )
            {
                for ( r = 0 ; r < nActive ; r++ )
                {
                    for ( q = 0 ; q < nActive ; q++ )
                    {
                        for ( f = 0.0e+00, p = 0 ; p < nActive ; p++ ) f += ( Array2D_Item ( orbitals, k, p+nCore ) * DoubleSymmetricMatrix_GetItem ( twoPDM, p, q, r, s, NULL ) ) ;
                        ArrayND_Item3D ( tPDM3, q, r, s ) = f ;
                    }
                }
            }
	    for ( l = k ; l < i0 + iData->norbitals ; ++l )
            {
	        ll = l * ( l + 1 ) / 2;
                /* . Second orbital transformation. */
                for ( s = 0 ; s < nActive ; s++ )
                {
                    for ( r = 0 ; r < nActive ; r++ )
                    {
                        for ( f = 0.0e+00, q = 0 ; q < nActive ; q++ ) f += ( Array2D_Item ( orbitals, l, q+nCore ) * ArrayND_Item3D ( tPDM3, q, r, s ) ) ;
                        Array2D_Item ( tPDM2, r, s ) = f ;
                    }
                }
	        kl0 = ( k - i0 ) + ( ( l - i0 ) * ( l - i0 + 1 ) ) / 2 ;
                kl  = k + ll ;
	        for ( m = j0 ; m < j0 + jData->norbitals ; ++m )
                {
		    bb = 1.0e+00 ;
                    /* . Third orbital transformation. */
                    for ( s = 0 ; s < nActive ; s++ )
                    {
                        for ( f = 0.0e+00, r = 0 ; r < nActive ; r++ ) f += ( Array2D_Item ( orbitals, m, r+nCore ) * Array2D_Item ( tPDM2, r, s ) ) ;
                        Array1D_Item ( tPDM1, s ) = f ;
                    }
		    for ( n = m ; n < j0 + jData->norbitals ; ++n )
                    {
		        mn0 = ( m - j0 ) + ( ( n - j0 ) * ( n - j0 + 1 ) ) / 2 ;
                        mn  = m + ( n * ( n + 1 ) ) / 2 ;
		        mk  = m + kk ;
		        nk  = n + kk ;
		        ml  = m + ll ;
		        nl  = n + ll ;
                        /* . dCore/dCore term. */
                        f1 = ( dCore->data[kl] * dCore->data[mn] ) - 0.25e+00 * ( dCore->data[mk] * dCore->data[nl] + dCore->data[nk] * dCore->data[ml] ) ;
                        /* . OnePDM/dCore term. */
                        f2 = 0.5e+00   * ( onePDM->data[kl] * dCore ->data[mn] + dCore ->data[kl] * onePDM->data[mn] ) -
                             0.125e+00 * ( onePDM->data[mk] * dCore ->data[nl] + onePDM->data[nk] * dCore ->data[ml] +
                                           dCore ->data[mk] * onePDM->data[nl] + dCore ->data[nk] * onePDM->data[ml] ) ;
                        /* . TwoPDM term. */
                        for ( f3 = 0.0e+00, s = 0 ; s < nActive ; s++ ) f3 += ( Array2D_Item ( orbitals, n, s+nCore ) * Array1D_Item ( tPDM1, s ) ) ;
                        /* . Zmatrix term. */
                        f4 = 0.5e+00   * ( dHF    ->data[kl] * zMatrix->data[mn] + zMatrix->data[kl] * dHF    ->data[mn] ) -
                             0.125e+00 * ( dHF    ->data[mk] * zMatrix->data[nl] + dHF    ->data[nk] * zMatrix->data[ml] +
                                           zMatrix->data[mk] * dHF    ->data[nl] + zMatrix->data[nk] * dHF    ->data[ml] ) ;
                        /* . Total contribution. */
                        Array2D_Item ( dIJ, kl0, mn0 ) = - aa * bb * ( f1 + 2.0e+00 * ( f2 + f3 + f4 ) ) ;
		        bb = 2.0e+00 ;
		    }
	        }
	        aa = 2.0e+00 ;
	    }
        }
    }
}
