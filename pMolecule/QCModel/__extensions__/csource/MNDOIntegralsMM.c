/*==================================================================================================================================
! . MNDO integrals for QC/MM interactions.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "Boolean.h"
# include "MNDODefinitions.h"
# include "MNDOIntegralDefinitions.h"
# include "MNDOIntegralsMM.h"
# include "MNDOIntegralUtilities.h"
# include "NumericalMacros.h"
# include "RealArray2D.h"
# include "Units.h"

# include "Status.h"

/* . All quantities are in atomic units. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real MNDOIntegralsMM_LocalFrame2COEI ( const ChargeInteractionFunction Evaluate ,
                                              const MNDOParameters           *qData    ,
                                              const Integer                   ij       ,
                                              const Integer                   i        ,
                                              const Integer                   j        ,
                                              const Real                      r        ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . MM parameters. */
# define _AlpMM   ( 5.0e+00 / Units_Length_Angstroms_To_Bohrs ) /* . From original QC/MM paper. */
# define _AlpQC   ( 3.0e+00 / Units_Length_Angstroms_To_Bohrs ) /* . Default when alp = 0 (e.g. for PM6). */
# define _GPhotMM   1.0e+00
# define _PO8MM     0.0e+00 /* . Same as Rho0MM. */
# define _Rho0MM    0.0e+00

/*----------------------------------------------------------------------------------------------------------------------------------
! . The QC/MM core-charge interaction.
! . Signed and unsigned portions are returned separately.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . PM6 diatomic parameters have not been optimized for these interactions. */
# define _MinimumR 0.1e+00
void MNDOIntegralsMM_CoreCharge ( const MNDOParameters *qData  ,
                                  const Real            qM     ,
                                  const Real            R      ,
                                        Real           *fCore0 ,
                                        Real           *fCore1 ,
                                        Real           *gCore0 ,
                                        Real           *gCore1 )
{
    Integer i ;
    Real    a, alpQC, anam1, d, dGam, exI, exJ, f0, f1, g0, g1, gam, r, scale, zQ, zQabs ;
    /* . Basic integral. */
    gam   = 1.0e+00 / sqrt ( R*R + pow ( ( qData->po[8] + _PO8MM ), 2 ) ) ;
    dGam  = -R * gam * gam * gam ;
    /* . Standard core terms. */
    alpQC = qData->alp ;
    if ( fabs ( alpQC ) < 1.0e-01 ) alpQC = _AlpQC ;
    exI   = exp ( - alpQC * R ) ;
    exJ   = exp ( -_AlpMM * R ) ;
    scale = exI + exJ ;
    zQ    = qData->zcore * qM ;
    zQabs = fabs ( zQ )       ;
    f0    = zQ    *    gam         ;
    f1    = zQabs *    gam * scale ;
    g0    = zQ    *   dGam         ;
    g1    = zQabs * ( dGam * scale - gam * ( alpQC * exI + _AlpMM * exJ ) ) ;
    /* . Compute the AM1/PM3-specific terms. */
    anam1 = 0.0e+00 ;
    scale = 0.0e+00 ;
    for ( i = 0 ; i < qData->nam1pm3g ; i++ )
    {
        r = ( R < _MinimumR ) ? _MinimumR : R ;
        d = r - qData->fn3[i] ;
        a = qData->fn2[i] * d * d ;
        if ( a <= EXPONENT_TOLERANCE )
        {
            exI    = qData->fn1[i] * exp( -a ) / r ;
            anam1 += ( 1.0e+00 / r + 2.0e+00 * qData->fn2[i] * d ) * exI ;
            scale += exI ;
        }
    }
    g0 -= anam1 * qData->gphot * _GPhotMM * zQ ;
    f0 += scale * qData->gphot * _GPhotMM * zQ ;
    /* . Finish up. */
    (*fCore0) = f0 ;
    (*fCore1) = f1 ;
    (*gCore0) = g0 ;
    (*gCore1) = g1 ;
}
# undef _MinimumR

/*----------------------------------------------------------------------------------------------------------------------------------
! . The QC/MM core-charge and electron-charge integrals and derivatives from splines.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Very little checking. */
void MNDOIntegralsMM_FromSpline ( const MNDOParameters *qData      ,
                                  const CubicSpline    *qSpline    ,
                                  const Real            qM         ,
                                  const Real            R          ,
                                        Real           *fCore      ,
                                        Real           *gCore      ,
                                        RealArray1D    *integrals  ,
                                        RealArray1D    *gIntegrals )
{
    Boolean doGradients = ( gIntegrals != NULL ) ;
    Integer am = 0, c, i, ij, j, l, n, u ;
    Real    d, f, f0, f1, g, g0, g1, s, t ;
    /* . Evaluate the interval given R. */
    CubicSpline_EvaluateLUDST  ( qSpline, R, &l, &u, &d, &s, &t ) ;
    /* . Core terms - signed and unsigned. */
    CubicSpline_FastEvaluateFGN ( qSpline, 0, l, u, d, s, t, f0, g0 ) ;
    CubicSpline_FastEvaluateFGN ( qSpline, 1, l, u, d, s, t, f1, g1 ) ;
    (*fCore) = qM * f0 + fabs ( qM ) * f1 ;
    (*gCore) = qM * g0 + fabs ( qM ) * g1 ;
    /* . Get the highest AM. */
    switch ( qData->norbitals )
    {
        case 1: am = 0 ; break ;
        case 4: am = 1 ; break ;
        case 9: am = 2 ; break ;
    }
    /* . Unique integrals. */
    for ( c = n = 0 ; c < _NCUniqueSPD[am] ; c++, n += 2 )
    {
        i  = _CUniqueSPD[n  ] ;
        j  = _CUniqueSPD[n+1] ;
        ij = ( i * ( i + 1 ) ) / 2 + j ;
        CubicSpline_FastEvaluateFGN ( qSpline, c+2, l, u, d, s, t, f, g ) ;
        Array1D_Item ( integrals, ij ) = f ;
        if ( doGradients ) Array1D_Item ( gIntegrals, ij ) = g ;
    }
    /* . Integrals related by symmetry. */
    for ( c = n = 0 ; c < NCPOSITIVE[am] ; c++, n+= 2 )
    {
        Array1D_Item ( integrals, CPOSITIVE[n] ) = Array1D_Item ( integrals, CPOSITIVE[n+1] ) ;
        if ( doGradients ) Array1D_Item ( gIntegrals, CPOSITIVE[n] ) = Array1D_Item ( gIntegrals, CPOSITIVE[n+1] ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The QC/MM electron-charge integrals and derivatives as a function of R in the local frame for unit MM charge.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOIntegralsMM_LocalFrame ( const MNDOParameters *qData      ,
                                  const Real            R          ,
                                        RealArray1D    *integrals  ,
                                        RealArray1D    *gIntegrals )
{
    Boolean  doGradients = ( gIntegrals != NULL ) ;
    Integer  c, i, iAM = 0, ij, j, t ;
    /* . Unique non-zero integrals in the local frame and then those related by symmetry. */
    /* . Get the highest AM. */
    switch ( qData->norbitals )
    {
        case 1: iAM = 0 ; break ;
        case 4: iAM = 1 ; break ;
        case 9: iAM = 2 ; break ;
    }
    /* . sp integrals. */
    MNDOIntegralUtilities_LocalFrame2COEIsSP ( qData, _PO8MM, R, False, integrals, gIntegrals ) ;
    /* . Other integrals. */
    for ( c = t = 0 ; c < NCUNIQUE[iAM] ; c++, t += 2 )
    {
        i  = CUNIQUE[t  ] ;
        j  = CUNIQUE[t+1] ;
        ij = ( i * ( i + 1 ) ) / 2 + j ;
        Array1D_Item ( integrals, ij ) = - MNDOIntegralsMM_LocalFrame2COEI ( &MNDOIntegralUtilities_2CChargeInteraction, qData, ij, ORBITALAM[i], ORBITALAM[j], R ) ;
    }
    for ( c = t = 0 ; c < NCPOSITIVE[iAM] ; c++, t+= 2 ) Array1D_Item ( integrals, CPOSITIVE[t] ) = Array1D_Item ( integrals, CPOSITIVE[t+1] ) ;
    /* . Derivatives. */
    if ( doGradients )
    {
        for ( c = t = 0 ; c < NCUNIQUE[iAM] ; c++, t += 2 )
        {
            i  = CUNIQUE[t  ] ;
            j  = CUNIQUE[t+1] ;
            ij = ( i * ( i + 1 ) ) / 2 + j ;
            Array1D_Item ( gIntegrals, ij ) = - MNDOIntegralsMM_LocalFrame2COEI ( &MNDOIntegralUtilities_2CChargeInteractionD, qData, ij, ORBITALAM[i], ORBITALAM[j], R ) ;
        }
        for ( c = t = 0 ; c < NCPOSITIVE[iAM] ; c++, t+= 2 ) Array1D_Item ( gIntegrals, CPOSITIVE[t] ) = Array1D_Item ( gIntegrals, CPOSITIVE[t+1] ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The QC/MM integrals and derivatives in the molecular frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOIntegralsMM_MolecularFrame ( const Integer      nQ          ,
                                      const Real         r           ,
                                      const Real         x           ,
                                      const Real         y           ,
                                      const Real         z           ,
                                      const RealArray1D *iLocal      ,
                                      const RealArray1D *gLocal      ,
                                            RealArray1D *iMolecular  ,
                                            RealArray1D *gMolecularX ,
                                            RealArray1D *gMolecularY ,
                                            RealArray1D *gMolecularZ )
{
    Boolean doGradients ;
    Real    dR[3] ;
    /* . Initialization. */
    doGradients = ( gLocal      != NULL ) &&
                  ( gMolecularX != NULL ) &&
                  ( gMolecularY != NULL ) &&
                  ( gMolecularZ != NULL ) ;
    dR[0] = - x / r ; dR[1] = - y / r ; dR[2] = - z / r ;
    /* . No transformation matrices - s functions. */
    if ( nQ == 1 )
    {
        Array1D_Item ( iMolecular, 0 ) = Array1D_Item ( iLocal, 0 ) ;
        if ( doGradients )
        {
            Array1D_Item ( gMolecularX, 0 ) = dR[0] * Array1D_Item ( gLocal, 0 ) ;
            Array1D_Item ( gMolecularY, 0 ) = dR[1] * Array1D_Item ( gLocal, 0 ) ;
            Array1D_Item ( gMolecularZ, 0 ) = dR[2] * Array1D_Item ( gLocal, 0 ) ;
        }
    }
    /* . Transformation matrices. */
    else
    {
        auto RealArray2D *iT = NULL, *iTx = NULL, *iTy = NULL, *iTz = NULL, **pTx = NULL, **pTy = NULL, **pTz = NULL ;
        if ( doGradients ) { pTx = &iTx ; pTy = &iTy ; pTz = &iTz ; }
        MNDOIntegralUtilities_GetTransformationMatrices ( nQ, 0, r, -x, -y, -z, &iT, NULL, pTx, pTy, pTz, NULL, NULL, NULL ) ;
        RealArray2D_VectorMultiply ( False, 1.0e+00, iT, iLocal, 0.0e+00, iMolecular, NULL ) ;
        if ( doGradients )
        {
            RealArray2D_VectorMultiply ( False, dR[0]  , iT , gLocal, 0.0e+00, gMolecularX, NULL ) ;
            RealArray2D_VectorMultiply ( False, 1.0e+00, iTx, iLocal, 1.0e+00, gMolecularX, NULL ) ;
            RealArray2D_VectorMultiply ( False, dR[1]  , iT , gLocal, 0.0e+00, gMolecularY, NULL ) ;
            RealArray2D_VectorMultiply ( False, 1.0e+00, iTy, iLocal, 1.0e+00, gMolecularY, NULL ) ;
            RealArray2D_VectorMultiply ( False, dR[2]  , iT , gLocal, 0.0e+00, gMolecularZ, NULL ) ;
            RealArray2D_VectorMultiply ( False, 1.0e+00, iTz, iLocal, 1.0e+00, gMolecularZ, NULL ) ;
        }
        RealArray2D_Deallocate ( &iT  ) ;
        RealArray2D_Deallocate ( &iTx ) ;
        RealArray2D_Deallocate ( &iTy ) ;
        RealArray2D_Deallocate ( &iTz ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Values of the core-charge and core-electron integrals and derivatives for unit MM charge.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Very little checking and not terribly efficient so as to reuse old code. */
void MNDOIntegralsMM_Values ( const MNDOParameters *qData ,
                              const Real            R     ,
                              const Integer         index ,
                                    Real           *f     ,
                                    Real           *g     )
{
    Boolean doG = ( g != NULL ) ;
    /* . Core-charge (index = 0,1). */
    if ( ( index == 0 ) || ( index == 1 ) )
    {
        auto Real f0, f1, g0, g1 ;
        MNDOIntegralsMM_CoreCharge ( qData, 1.0e+00, R, &f0, &f1, &g0, &g1 ) ;
        if ( index == 0 ) { if ( doG ) (*g) = g0 ; else (*f) = f0 ; }
        else              { if ( doG ) (*g) = g1 ; else (*f) = f1 ; }
    }
    /* . Electron-charge. */
    /* . sp (index = 2,3,4,5). */
    else if ( ( index >= 2 ) && ( index <= 5 ) )
    {
        auto RealArray1D *gIntegrals, *integrals ;
        gIntegrals = RealArray1D_AllocateWithExtent ( 10, NULL ) ;
        integrals  = RealArray1D_AllocateWithExtent ( 10, NULL ) ;
        MNDOIntegralUtilities_LocalFrame2COEIsSP ( qData, _PO8MM, R, False, integrals, gIntegrals ) ;
        if ( doG )
        {
            switch ( index )
            {
                case 2: (*g) = Array1D_Item ( gIntegrals, SS   ) ; break ;
                case 3: (*g) = Array1D_Item ( gIntegrals, PZS  ) ; break ;
                case 4: (*g) = Array1D_Item ( gIntegrals, PZPZ ) ; break ;
                case 5: (*g) = Array1D_Item ( gIntegrals, PXPX ) ; break ;
            }
        }
        else
        {
            switch ( index )
            {
                case 2: (*f) = Array1D_Item ( integrals, SS   ) ; break ;
                case 3: (*f) = Array1D_Item ( integrals, PZS  ) ; break ;
                case 4: (*f) = Array1D_Item ( integrals, PZPZ ) ; break ;
                case 5: (*f) = Array1D_Item ( integrals, PXPX ) ; break ;
            }
        }
        RealArray1D_Deallocate ( &gIntegrals ) ;
        RealArray1D_Deallocate ( & integrals ) ;
    }
    /* . d (index = 6,7,8,9,10,11). */
    else
    {
        auto Integer i, ij, j, t ;
        t  = 2 * ( index - 6 ) ;
        i  = CUNIQUE[t  ] ;
        j  = CUNIQUE[t+1] ;
        ij = ( i * ( i + 1 ) ) / 2 + j ;
        if ( doG ) (*g) = - MNDOIntegralsMM_LocalFrame2COEI ( &MNDOIntegralUtilities_2CChargeInteractionD, qData, ij, ORBITALAM[i], ORBITALAM[j], R ) ;
        else       (*f) = - MNDOIntegralsMM_LocalFrame2COEI ( &MNDOIntegralUtilities_2CChargeInteraction , qData, ij, ORBITALAM[i], ORBITALAM[j], R ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a two-center OEI or its derivative in the local frame (d orbitals only).
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real MNDOIntegralsMM_LocalFrame2COEI ( const ChargeInteractionFunction Evaluate ,
                                              const MNDOParameters           *qData    ,
                                              const Integer                   ij       ,
                                              const Integer                   i        ,
                                              const Integer                   j        ,
                                              const Real                      r        )
{
    Real integral = 0.0e+00 ;
    if ( ( qData != NULL ) && ( NCHTERMS[ij] > 0 ) )
    {
        auto Integer lij, lmin, lm1, lm2, l1, l1Max, l1Min, l1Offset, l1Terms[MAXCHTERMS], m, mTerms[MAXCHTERMS], nTerms, t ;
        auto Real    add, chijkl[MAXCHTERMS], dij = 0.0e+00, pij = 0.0e+00 ;
        /* . Determine the number of terms. */
        l1Min = i - j ;
        l1Max = Minimum ( i + j, 2 ) ;
        lij   = ( i * ( i + 1 ) ) / 2 + j ;
        for ( l1 = l1Min, nTerms = 0 ; l1 <= l1Max ; l1++ )
        {
            l1Offset = ij * CHINCREMENT1 + l1 * CHINCREMENT2 + CHINCREMENT3 ;
            lmin     = Minimum ( l1, 0 ) ;
            for ( m = -lmin ; m <= lmin ; m++ )
            {
                lm1 = CHINDICES[l1Offset+m    ] ;
                lm2 = CHINDICES[CHINCREMENT3+m] ;
                if ( ( lm1 > 0 ) && ( lm2 > 0 ) )
                {
                    l1Terms[nTerms] = l1 ;
                    mTerms [nTerms] = abs ( m ) ;
                    chijkl [nTerms] = CHTERMS[lm1-1] * CHTERMS[lm2-1] ;
                    nTerms ++ ;
                }
            }
        }
        /* . Calculate the terms. */
        for ( t = 0 ; t < nTerms ; t++ )
        {
            l1 = l1Terms[t] ;
            m  = mTerms [t] ;
            if ( l1 == 0 )
            {
                dij = 0.0e+00 ;
                switch ( i )
                {
                    case 0:
                        pij = qData->po[0] ;
                        break ;
                    case 1:
                        pij = qData->po[6] ;
                        break ;
                    case 2:
                        pij = qData->po[7] ;
                        break ;
                }
            }
            else
            {
                dij = qData->ddp[lij] ;
                pij = qData->po [lij] ;
            }
            add = pow ( ( pij + _PO8MM ), 2 ) ;
            integral += chijkl[t] * Evaluate ( r, l1, 0, m, dij, 0.0e+00, add ) ;
        }
    }
    return integral ;
}


