/*==================================================================================================================================
! . Procedures for calculating the integrals in a MNDO method.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "BlockStorage.h"
# include "Memory.h"
# include "MNDODefinitions.h"
# include "MNDOIntegralDefinitions.h"
# include "MNDOIntegrals.h"
# include "MNDOIntegralUtilities.h"
# include "MNDOParameters.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add in the one-center TEIs.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The maximum number of unique integrals - 1 (s), 16 (sp), 155 (spd). */
void MNDOIntegrals_AddInOneCenterTEIs ( const MNDOParameters *self                 ,
                                        const Integer         i0                   ,
                                              BlockStorage   *twoElectronIntegrals )
{
    if ( self->nocteis > 0 )
    {
        auto Integer    i ;
        auto Cardinal16 indices[4*N1CTEIS] ;
        for ( i = 0 ; i < 4*self->nocteis ; i++ ) indices[i] = self->octeiindices[i] + i0 ;
        BlockStorage_AddData ( twoElectronIntegrals, self->nocteis, self->octeivalues, indices, NULL, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the integrals in the molecular frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOIntegrals_MolecularFrame2CIntegrals ( const MNDOParameters *iData                ,
                                               const Integer         i0                   ,
                                               const Real           *xI                   ,
                                               const MNDOParameters *jData                ,
                                               const Integer         j0                   ,
                                               const Real           *xJ                   ,
                                                     RealArray1D    *mfcore1b             ,
                                                     RealArray1D    *mfcore2a             ,
                                                     BlockStorage   *twoElectronIntegrals )
{
    if ( ( iData != NULL ) && ( jData != NULL ) && ( mfcore1b != NULL ) && ( mfcore2a != NULL ) && ( twoElectronIntegrals != NULL ) )
    {
        auto Cardinal16   indices[4*N2CTEIS] ;
        auto Integer      i, j, k, l, n, ni, nj ;
        auto Real         r, x, y, z ;
        auto RealArray2D *hfteis = NULL, *iTransformation = NULL, *jTransformation = NULL, *lfteis = NULL, *mfteis = NULL ;
        auto RealArray1D *lfcore1b = NULL, *lfcore2a = NULL ;
        /* . Get the transformation matrices. */
        MNDOIntegralUtilities_GetDisplacement ( xI, xJ, r, x, y, z ) ;
        MNDOIntegralUtilities_GetTransformationMatrices ( iData->norbitals, jData->norbitals, r, x, y, z, &iTransformation, &jTransformation, NULL, NULL, NULL, NULL, NULL, NULL ) ;
        /* . Allocate space for the integrals. */
        ni = ( iData->norbitals * ( iData->norbitals + 1 ) ) / 2 ;
        nj = ( jData->norbitals * ( jData->norbitals + 1 ) ) / 2 ;
        lfteis = RealArray2D_AllocateWithExtents ( ni, nj, NULL ) ;
        if ( iTransformation == NULL ) lfcore1b = mfcore1b ;
        else                           lfcore1b = RealArray1D_AllocateWithExtent ( ni, NULL ) ;
        if ( jTransformation == NULL ) lfcore2a = mfcore2a ;
        else                           lfcore2a = RealArray1D_AllocateWithExtent ( nj, NULL ) ;
        /* . Get the integrals in the local frame. */
        MNDOIntegralUtilities_LocalFrame2CTEIs ( iData, jData, r, lfteis, lfcore1b, lfcore2a, NULL, NULL, NULL ) ;
        /* . Transform from the local to molecular frames. */
        /* . OEIs then TEIs. */
        if ( iTransformation == NULL )
        {
            lfcore1b = NULL   ;
            hfteis   = lfteis ;
            lfteis   = NULL   ;
        }
        else
        {
            RealArray2D_VectorMultiply ( False, 1.0e+00, iTransformation, lfcore1b, 0.0e+00, mfcore1b, NULL ) ;
            hfteis   = RealArray2D_AllocateWithExtents ( ni, nj, NULL ) ;
            RealArray2D_MatrixMultiply ( False, False, 1.0e+00, iTransformation, lfteis, 0.0e+00, hfteis, NULL ) ;
        }

        if ( jTransformation == NULL )
        {
            lfcore2a = NULL     ;
            mfteis   = hfteis   ;
            hfteis   = NULL     ;
        }
        else
        {
            RealArray2D_VectorMultiply ( False, 1.0e+00, jTransformation, lfcore2a, 0.0e+00, mfcore2a, NULL ) ;
            mfteis   = RealArray2D_AllocateWithExtents ( ni, nj, NULL ) ;
            RealArray2D_MatrixMultiply ( False, True, 1.0e+00, hfteis, jTransformation, 0.0e+00, mfteis, NULL ) ;
        }
        /* . Determine the TEI indices. */
        /* . There is no restriction on i, j, k and l as their order is checked when building the Fock matrices. */
        n = 0 ;
        for ( i = 0 ; i < iData->norbitals ; i++ )
        {
            for ( j = 0 ; j <= i ; j++ )
            {
                for ( k = 0 ; k < jData->norbitals ; k++ )
                {
                    for ( l = 0 ; l <= k ; l++ )
                    {
                        indices[n  ] = i + i0 ;
                        indices[n+1] = j + i0 ;
                        indices[n+2] = k + j0 ;
                        indices[n+3] = l + j0 ;
                        n += 4 ;
                    }
                }
            }
        }
        /* . Save the data. */
        BlockStorage_AddData ( twoElectronIntegrals, ( ni * nj ), mfteis->data, indices, NULL, NULL ) ;
        /* . Finish up. */
        RealArray2D_Deallocate ( &hfteis   ) ;
        RealArray2D_Deallocate ( &lfteis   ) ;
        RealArray2D_Deallocate ( &mfteis   ) ;
        RealArray1D_Deallocate ( &lfcore1b ) ;
        RealArray1D_Deallocate ( &lfcore2a ) ;
        /* . The transformations can be identical so be careful about deallocation. */
        if ( iTransformation != jTransformation ) RealArray2D_Deallocate ( &iTransformation ) ;
        RealArray2D_Deallocate ( &jTransformation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the derivatives in the molecular frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOIntegrals_MolecularFrame2CIntegralsD ( const MNDOParameters *iData  ,
                                                const Integer         i0     ,
                                                const Real           *xI     ,
                                                const MNDOParameters *jData  ,
                                                const Integer         j0     ,
                                                const Real           *xJ     ,
                                                const RealArray1D    *dOneI  ,
                                                const RealArray1D    *dOneJ  ,
                                                const RealArray2D    *dTwoIJ ,
                                                      Real           *gX     ,
                                                      Real           *gY     ,
                                                      Real           *gZ     )
{
    Boolean      doI, doJ ;
    Integer      i, ix, j, ni, nj ;
    Real         doei0, doei0f, doei1, dR[3], dtei0, dtei0f, dtei1, eng[3], factor, r, x, y, z ;
    RealArray1D *dlfcore1b = NULL, *dlfcore2a = NULL, *lfcore1b = NULL, *lfcore2a = NULL, *temporaryi = NULL, *temporaryj = NULL ;
    RealArray2D *dhfteis = NULL, *dlfteis = NULL, *lfteis = NULL, *LTj = NULL, *temporaryij = NULL, *TiL = NULL ;
    RealArray2D *iTransformation = NULL, *iTransformationD[3], *iTransformationX = NULL, *iTransformationY = NULL, *iTransformationZ = NULL,
                *jTransformation = NULL, *jTransformationD[3], *jTransformationX = NULL, *jTransformationY = NULL, *jTransformationZ = NULL ;
    /* . Initialization. */
    for ( i = 0 ; i < 3; i++ ) eng[i] = 0.0e+00;
    /* . Get the transformation matrices. */
    MNDOIntegralUtilities_GetDisplacement ( xI, xJ, r, x, y, z ) ;
    MNDOIntegralUtilities_GetTransformationMatrices ( iData->norbitals  ,
                                                      jData->norbitals  ,
                                                      r, x, y, z ,
                                                      &iTransformation  ,
                                                      &jTransformation  ,
                                                      &iTransformationX ,
                                                      &iTransformationY ,
                                                      &iTransformationZ ,
                                                      &jTransformationX ,
                                                      &jTransformationY ,
                                                      &jTransformationZ ) ;
    /* . Allocate space. */
    ni = ( iData->norbitals * ( iData->norbitals + 1 ) ) / 2 ;
    nj = ( jData->norbitals * ( jData->norbitals + 1 ) ) / 2 ;
    dlfteis   = RealArray2D_AllocateWithExtents ( ni, nj, NULL ) ; RealArray2D_Set ( dlfteis   , 0.0e+00 ) ;
    lfteis    = RealArray2D_AllocateWithExtents ( ni, nj, NULL ) ; RealArray2D_Set ( lfteis    , 0.0e+00 ) ;
    dlfcore1b = RealArray1D_AllocateWithExtent  ( ni    , NULL ) ; RealArray1D_Set ( dlfcore1b , 0.0e+00 ) ;
    dlfcore2a = RealArray1D_AllocateWithExtent  ( nj    , NULL ) ; RealArray1D_Set ( dlfcore2a , 0.0e+00 ) ;
    lfcore1b  = RealArray1D_AllocateWithExtent  ( ni    , NULL ) ; RealArray1D_Set ( lfcore1b  , 0.0e+00 ) ;
    lfcore2a  = RealArray1D_AllocateWithExtent  ( nj    , NULL ) ; RealArray1D_Set ( lfcore2a  , 0.0e+00 ) ;
    /* . Define some transformation factors. */
    dR[0] = x / r ; dR[1] = y / r ; dR[2] = z / r ;
    iTransformationD[0] = iTransformationX ; iTransformationD[1] = iTransformationY ; iTransformationD[2] = iTransformationZ ;
    jTransformationD[0] = jTransformationX ; jTransformationD[1] = jTransformationY ; jTransformationD[2] = jTransformationZ ;
    /* . Compute the integrals and derivatives in the local frame. */
    MNDOIntegralUtilities_LocalFrame2CTEIs ( iData, jData, r, lfteis, lfcore1b, lfcore2a, dlfteis, dlfcore1b, dlfcore2a ) ;
    /* . Set some flags. */
    doI = ( iTransformation != NULL ) ;
    doJ = ( jTransformation != NULL ) ;
    /* . Local frame terms (only depend on r). */
    /* . OEIs. */
    doei0f = 0.0e+00 ;
    if ( doI )
    {
        temporaryi = RealArray1D_AllocateWithExtent ( ni, NULL ) ;
        RealArray2D_VectorMultiply ( False, 1.0e+00, iTransformation, dlfcore1b, 0.0e+00, temporaryi, NULL ) ;
    }
    else
    {
        temporaryi = dlfcore1b ;
        dlfcore1b  = NULL      ;
    }
    factor = RealArray1D_Dot ( dOneI, temporaryi, NULL ) ; doei0f += factor ;
    if ( doJ )
    {
        temporaryj = RealArray1D_AllocateWithExtent ( nj, NULL ) ;
        RealArray2D_VectorMultiply ( False, 1.0e+00, jTransformation, dlfcore2a, 0.0e+00, temporaryj, NULL ) ;
    }
    else
    {
        temporaryj = dlfcore2a ;
        dlfcore2a  = NULL      ;
    }
    factor = RealArray1D_Dot ( dOneJ, temporaryj, NULL ) ; doei0f += factor ;
    /* . TEIs. */
    dtei0f = 0.0e+00 ;
    if ( doI )
    {
        dhfteis = RealArray2D_AllocateWithExtents ( ni, nj, NULL ) ;
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, iTransformation, dlfteis, 0.0e+00, dhfteis, NULL ) ;
    }
    else
    {
        dhfteis = dlfteis ;
        dlfteis = NULL    ;
    }
    if ( doJ )
    {
        temporaryij = RealArray2D_AllocateWithExtents ( ni, nj, NULL ) ;
        RealArray2D_MatrixMultiply ( False, True, 1.0e+00, dhfteis, jTransformation, 0.0e+00, temporaryij, NULL ) ;
    }
    else
    {
        temporaryij = dhfteis ;
        dhfteis     = NULL   ;
    }
    /* . Coulomb and exchange terms. */
    for ( i = 0 ; i < ni ; i++ )
    {
        for ( j = 0 ; j < nj ; j++ ) dtei0f -= Array2D_Item ( dTwoIJ, i, j ) * Array2D_Item ( temporaryij, i, j ) ;
    }
    /* . Determine some intermediate matrices for later. */
    if ( doI )
    {
        LTj = RealArray2D_AllocateWithExtents ( ni, nj, NULL ) ;
        if ( doJ ) RealArray2D_MatrixMultiply ( False, True, 1.0e+00, lfteis, jTransformation, 0.0e+00, LTj, NULL ) ;
        else       RealArray2D_CopyTo ( lfteis, LTj, NULL ) ;
    }
    if ( doJ )
    {
        TiL = RealArray2D_AllocateWithExtents ( ni, nj, NULL ) ;
        if ( doI ) RealArray2D_MatrixMultiply ( False, False, 1.0e+00, iTransformation, lfteis, 0.0e+00, TiL, NULL ) ;
        else       RealArray2D_CopyTo ( lfteis, TiL, NULL ) ;
    }
    /* . Electronic terms. */
    /* . Loop over the Cartesian components. */
    for ( ix = 0 ; ix < 3 ; ++ix )
    {
        /* . OEI. */
        doei0 = - dR[ix] * doei0f ;
        doei1 = 0.0e+00 ;
        if ( doI )
        {
            RealArray2D_VectorMultiply ( False, 1.0e+00, iTransformationD[ix], lfcore1b, 0.0e+00, temporaryi, NULL ) ;
            factor = RealArray1D_Dot ( dOneI, temporaryi, NULL ) ; doei1 -= factor ;
        }
        if ( doJ )
        {
            RealArray2D_VectorMultiply ( False, 1.0e+00, jTransformationD[ix], lfcore2a, 0.0e+00, temporaryj, NULL ) ;
            factor = RealArray1D_Dot ( dOneJ, temporaryj, NULL ) ; doei1 -= factor ;
        }

        /* . TEI. */
        dtei0 = - dR[ix] * dtei0f ;
        dtei1 = 0.0e+00 ;
        if ( doI || doJ )
        {
            /* . Get the total integrals. */
            factor = 0.0e+00 ;
            if ( doI ) { RealArray2D_MatrixMultiply ( False, False, 1.0e+00, iTransformationD[ix], LTj, 0.0e+00, temporaryij, NULL ) ; factor = 1.0e+00 ; }
            if ( doJ ) { RealArray2D_MatrixMultiply ( False, True,   1.0e+00, TiL, jTransformationD[ix], factor,  temporaryij, NULL ) ; }

            /* . Coulomb and exchange terms. */
            for ( i = 0 ; i < ni ; i++ )
            {
                for ( j = 0 ; j < nj ; j++ ) dtei1 += Array2D_Item ( dTwoIJ, i, j ) * Array2D_Item ( temporaryij, i, j ) ;
            }
        }

        /* . Save the gradient terms. */
        eng[ix] = doei0 + doei1 + dtei0 + dtei1 ;
    }
    /* . Clear up. */
    RealArray2D_Deallocate ( &dhfteis     ) ;
    RealArray2D_Deallocate ( &dlfteis     ) ;
    RealArray2D_Deallocate ( &lfteis      ) ;
    RealArray2D_Deallocate ( &LTj         ) ;
    RealArray2D_Deallocate ( &TiL         ) ;
    RealArray2D_Deallocate ( &temporaryij ) ;
    RealArray1D_Deallocate ( &dlfcore1b   ) ;
    RealArray1D_Deallocate ( &dlfcore2a   ) ;
    RealArray1D_Deallocate ( &lfcore1b    ) ;
    RealArray1D_Deallocate ( &lfcore2a    ) ;
    RealArray1D_Deallocate ( &temporaryi  ) ;
    RealArray1D_Deallocate ( &temporaryj  ) ;
    /* . The transformations can be identical so be careful about deallocation. */
    if ( iTransformation  != jTransformation  ) RealArray2D_Deallocate ( &iTransformation  ) ;
    if ( iTransformationX != jTransformationX ) RealArray2D_Deallocate ( &iTransformationX ) ;
    if ( iTransformationY != jTransformationY ) RealArray2D_Deallocate ( &iTransformationY ) ;
    if ( iTransformationZ != jTransformationZ ) RealArray2D_Deallocate ( &iTransformationZ ) ;
    RealArray2D_Deallocate ( &jTransformation  ) ;
    RealArray2D_Deallocate ( &jTransformationX ) ;
    RealArray2D_Deallocate ( &jTransformationY ) ;
    RealArray2D_Deallocate ( &jTransformationZ ) ;
    /* . Finish up. */
    (*gX) = eng[0] ;
    (*gY) = eng[1] ;
    (*gZ) = eng[2] ;
}
