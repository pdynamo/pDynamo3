/*==================================================================================================================================
! . Integrals - 4 basis, 2 electrons, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "Boolean.h"
# include "GaussianBasisIntegrals_b4e2n0.h"
# include "GaussianBasisSubsidiary.h"
# include "Integer.h"
# include "NumericalMacros.h"
# include "Real.h"
# include "RysQuadrature.h"

/*# define _PrintIntegrals*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the two-electron integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _IntegralSize ( ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP1 ) * MAXAMP1 * ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP1 ) * MAXAMP1 * _MAXRYS )
void GaussianBasisIntegrals_TEIs ( const GaussianBasis *iBasis     ,
                                   const Real          *rI         ,
                                   const GaussianBasis *jBasis     ,
                                   const Real          *rJ         ,
                                   const Real          *rIJ        ,
                                   const Real           rIJ2       ,
                                   const GaussianBasis *kBasis     ,
                                   const Real          *rK         ,
                                   const GaussianBasis *lBasis     ,
                                   const Real          *rL         ,
                                   const Real          *rKL        ,
                                   const Real           rKL2       ,
                                   const Boolean        jLessThanL ,
                                         Block         *block      )
{
    /* . Pair data. */
    Boolean        iIsJ, iIsK, jIsL, kIsL ;
    /* . Shell loops. */
    Boolean        iAndJ  , ijAndKL, kAndL, qIJ0, qIJ1, qKL0, qKL1 ;
    Integer        iAMMax , iAMMaxt, iCBFInd, iShell,         nCFuncI ,
                   jAMMax , jAMMaxt, jCBFInd, jShell, jUpper, nCFuncJ ,
                   kAMMax , kAMMaxt, kCBFInd, kShell, kUpper, nCFuncK ,
                   lAMMax , lAMMaxt, lCBFInd, lShell, lUpper, nCFuncL ,
                   mAMMax , nAMMax ,
                   nRoots , strideI, strideIt, strideJ, strideJt ,
                   strideK, strideKt, strideL, strideLt, strideM ;
    Real           xIJt, xKLt, yIJt, yKLt, zIJt, zKLt ;
    Real           g[MAXCBF*MAXCBF*MAXCBF*MAXCBF] ;
    const Real    *rC, *rD ;
    /* . Primitive loops. */
    Integer        iP, jP, kP, lP ;
    Real           aa, aaInv, aAndB, ab, aI, aJ, aK, aL, ar2I, ar2K, arg, argIJ, bb, bbInv,
                   axac, ayac, azac, axad, ayad, azad, bxbc, bybc, bzbc, bxbd, bybd, bzbd,
                   c1x , c2x , c3x , c4x , c1y , c2y , c3y , c4y , c1z , c2z , c3z , c4z,
                   expFac, rho, xAB , yAB , zAB ;
    Real           arI[3], arK[3], rA[3], rB[3], xInt[_IntegralSize], yInt[_IntegralSize], zInt[_IntegralSize] ;
    RysQuadrature  roots ;
# ifdef _PrintIntegrals
auto Integer ntotal = 0 ;
# endif
    /* . Initialization. */
    block->count = 0 ;
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    iIsK = ( iBasis == kBasis ) && ( rI == rK ) ;
    jIsL = ( jBasis == lBasis ) && ( rJ == rL ) ;
    kIsL = ( kBasis == lBasis ) && ( rK == rL ) ;
    /* . Quadruple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        iAMMax  = iBasis->shells[iShell].type->angularmomentum_high ;
        iCBFInd = iBasis->shells[iShell].type->cbfindex ;
        nCFuncI = iBasis->shells[iShell].type->ncbf     ;
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nshells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jAMMax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jCBFInd = jBasis->shells[jShell].type->cbfindex ;
            nCFuncJ = jBasis->shells[jShell].type->ncbf     ;
            nAMMax  = iAMMax + jAMMax ;
            if ( iAMMax >= jAMMax )
            {
               iAMMaxt = iAMMax ;
               jAMMaxt = jAMMax ;
               xIJt    = rIJ[0] ;
               yIJt    = rIJ[1] ;
               zIJt    = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iAMMaxt = jAMMax ;
               jAMMaxt = iAMMax ;
               xIJt    = - rIJ[0] ;
               yIJt    = - rIJ[1] ;
               zIJt    = - rIJ[2] ;
               rC      = rJ ;
            }
            iAndJ = iIsJ && ( iShell == jShell ) ;
            qIJ0  = ( nAMMax == 0 ) ;
            qIJ1  = ( nAMMax <= 1 ) ;
            if ( iIsK ) kUpper = iShell + 1 ;
            else        kUpper = kBasis->nshells ;
            for ( kShell = 0 ; kShell < kUpper ; kShell++ )
            {
                kAMMax  = kBasis->shells[kShell].type->angularmomentum_high ;
                kCBFInd = kBasis->shells[kShell].type->cbfindex ;
                nCFuncK = kBasis->shells[kShell].type->ncbf     ;
                if ( iIsK && ( iShell == kShell ) )
                {
                          if ( jIsL       ) { lUpper = jShell + 1      ; }
                     else if ( jLessThanL ) { lUpper = -1              ; }
                     else                   { lUpper = lBasis->nshells ; }
                }
                else
                {    if ( kIsL ) { lUpper = kShell + 1      ; }
                     else        { lUpper = lBasis->nshells ; }
                }
                for ( lShell = 0 ; lShell < lUpper ; lShell++ )
                {
                    auto Integer  i ;
                    lAMMax  = lBasis->shells[lShell].type->angularmomentum_high ;
                    lCBFInd = lBasis->shells[lShell].type->cbfindex ;
                    nCFuncL = lBasis->shells[lShell].type->ncbf     ;
                    mAMMax  = kAMMax + lAMMax ;
                    if ( kAMMax >= lAMMax )
                    {
                       kAMMaxt = kAMMax ;
                       lAMMaxt = lAMMax ;
                       xKLt    = rKL[0] ;
                       yKLt    = rKL[1] ;
                       zKLt    = rKL[2] ;
                       rD      = rK ;
                    }
                    else
                    {
                       kAMMaxt = lAMMax ;
                       lAMMaxt = kAMMax ;
                       xKLt    = - rKL[0] ;
                       yKLt    = - rKL[1] ;
                       zKLt    = - rKL[2] ;
                       rD      = rL ;
                    }
                    kAndL   = kIsL && ( kShell == lShell ) ;
                    ijAndKL = iIsK && ( iShell == kShell ) && jIsL && ( jShell == lShell ) ;
                    qKL0    = ( mAMMax == 0 ) ;
                    qKL1    = ( mAMMax <= 1 ) ;
                    /* . Integral initialization. */
                    nRoots  = ( mAMMax + nAMMax ) / 2 + 1 ;
                    for ( i = 0 ; i < ( nCFuncI * nCFuncJ * nCFuncK * nCFuncL ) ; i++ ) g[i] = 0.0e+00 ;
                    strideL = 1 ;
                    strideK = ( lAMMaxt + 1 ) * strideL ;
                    strideJ = ( mAMMax  + 1 ) * strideK ;
                    strideI = ( jAMMaxt + 1 ) * strideJ ;
                    strideM = ( nAMMax  + 1 ) * strideI ;
                    if ( iAMMax >= jAMMax ) { strideIt = strideI ; strideJt = strideJ ; }
                    else                    { strideIt = strideJ ; strideJt = strideI ; }
                    if ( kAMMax >= lAMMax ) { strideKt = strideK ; strideLt = strideL ; }
                    else                    { strideKt = strideL ; strideLt = strideK ; }
    /*------------------------------------------------------------------------------------------------------------------------------
    ! . Quadruple loop over primitives.
    !-----------------------------------------------------------------------------------------------------------------------------*/
    for ( iP = 0 ; iP < iBasis->shells[iShell].nprimitives ; iP++ )
    {
	aI   = iBasis->shells[iShell].primitives[iP].exponent ;
	ar2I = aI * rIJ2 ;
	for ( i = 0 ; i < 3 ; i++ ) arI[i] = aI * rI[i] ;
        for ( jP = 0 ; jP < jBasis->shells[jShell].nprimitives ; jP++ )
        {
	    aJ    = jBasis->shells[jShell].primitives[jP].exponent ;
	    aa    = aI + aJ ;
	    aaInv = 1.0e+00 / aa ;
	    argIJ = aJ * ar2I * aaInv ;
	    if ( argIJ > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
	    for ( i = 0 ; i < 3 ; i++ ) rA[i] = ( arI[i] + aJ * rJ[i] ) * aaInv ;
            axad = aa * ( rA[0] - rD[0] ) ;
            ayad = aa * ( rA[1] - rD[1] ) ;
            azad = aa * ( rA[2] - rD[2] ) ;
            axac = aa * ( rA[0] - rC[0] ) ;
            ayac = aa * ( rA[1] - rC[1] ) ;
            azac = aa * ( rA[2] - rC[2] ) ;
            for ( kP = 0 ; kP < kBasis->shells[kShell].nprimitives ; kP++ )
            {
	        aK   = kBasis->shells[kShell].primitives[kP].exponent ;
	        ar2K = aK * rKL2 ;
	        for ( i = 0 ; i < 3 ; i++ ) arK[i] = aK * rK[i] ;
                for ( lP = 0 ; lP < lBasis->shells[lShell].nprimitives ; lP++ )
                {
                    auto Integer  m ;
	            aL     = lBasis->shells[lShell].primitives[lP].exponent ;
	            bb     = aK + aL ;
	            bbInv  = 1.0e+00 / bb ;
	            arg    = argIJ + aL * ar2K * bbInv ;
	            if ( arg > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    ab     = aa * bb ;
                    aAndB  = aa + bb ;
                    rho    = ab / aAndB ;
                    expFac = exp ( - arg ) * PI252 / ( ab * sqrt ( aAndB ) ) ;
	            for ( i = 0 ; i < 3 ; i++ ) rB[i] = ( arK[i] + aL * rL[i] ) * bbInv ;
                    bxbd = bb * ( rB[0] - rD[0] ) ;
                    bybd = bb * ( rB[1] - rD[1] ) ;
                    bzbd = bb * ( rB[2] - rD[2] ) ;
                    bxbc = bb * ( rB[0] - rC[0] ) ;
                    bybc = bb * ( rB[1] - rC[1] ) ;
                    bzbc = bb * ( rB[2] - rC[2] ) ;
                    c1x  = bxbd + axad ;
                    c2x  = aa * bxbd   ;
                    c3x  = bxbc + axac ;
                    c4x  = bb * axac   ;
                    c1y  = bybd + ayad ;
                    c2y  = aa * bybd   ;
                    c3y  = bybc + ayac ;
                    c4y  = bb * ayac   ;
                    c1z  = bzbd + azad ;
                    c2z  = aa * bzbd   ;
                    c3z  = bzbc + azac ;
                    c4z  = bb * azac   ;
                    xAB  = rA[0] - rB[0] ;
                    yAB  = rA[1] - rB[1] ;
                    zAB  = rA[2] - rB[2] ;
                    RysQuadrature_Roots ( &roots, nRoots, rho * ( xAB*xAB + yAB*yAB + zAB*zAB ) ) ;
                    for ( m = 0 ; m < nRoots ; m++ )
                    {
                        auto Real b00, b10, bp01, f00, fac, fac2, u2, xc00, xcp00, yc00, ycp00, zc00, zcp00 ;
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] * expFac ;
                        fac   = 1.0e+00 / ( ab + u2 * aAndB ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( aa   + u2 ) * fac2 ;
                        b00   =          u2   * fac2 ;
                        b10   = ( bb   + u2 ) * fac2 ;
                        xcp00 = ( u2 * c1x + c2x ) * fac ;
                        ycp00 = ( u2 * c1y + c2y ) * fac ;
                        zcp00 = ( u2 * c1z + c2z ) * fac ;
                        xc00  = ( u2 * c3x + c4x ) * fac ;
                        yc00  = ( u2 * c3y + c4y ) * fac ;
                        zc00  = ( u2 * c3z + c4z ) * fac ;
                        Subsidiary_Integral_Nuclear4C ( iAMMaxt ,
                                                        jAMMaxt ,
                                                        nAMMax  ,
                                                        kAMMaxt ,
                                                        lAMMaxt ,
                                                        mAMMax  ,
                                                        qIJ0    ,
                                                        qIJ1    ,
                                                        qKL0    ,
                                                        qKL1    ,
                                                        b00     ,
                                                        b10     ,
                                                        bp01    ,
                                                        xIJt    ,
                                                        yIJt    ,
                                                        zIJt    ,
                                                        xKLt    ,
                                                        yKLt    ,
                                                        zKLt    ,
                                                        f00     ,
                                                        xc00    ,
                                                        xcp00   ,
                                                        yc00    ,
                                                        ycp00   ,
                                                        zc00    ,
                                                        zcp00   ,
                                                        strideI ,
                                                        strideJ ,
                                                        strideK ,
                                                        &xInt[m*strideM] ,
                                                        &yInt[m*strideM] ,
                                                        &zInt[m*strideM] ) ;
                    }
                    /* . Assemble the integrals. */
                    {
                        auto Integer  i, ix, iy, iz, j, jix, jiy, jiz, k, kjix, kjiy, kjiz, l, lkjix, lkjiy, lkjiz, m, n ;
                        auto Real    fac, ti, tij, tijk, tijkl ;
                        for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                        {
   	                    ix = CBFPOWX[i+iCBFInd]*strideIt ;
	                    iy = CBFPOWY[i+iCBFInd]*strideIt ;
	                    iz = CBFPOWZ[i+iCBFInd]*strideIt ;
                            ti = iBasis->shells[iShell].primitives[iP].ccbf[i] ;
                            for ( j = 0 ; j < nCFuncJ ; j++ )
                            {
	                        jix = CBFPOWX[j+jCBFInd]*strideJt + ix ;
	                        jiy = CBFPOWY[j+jCBFInd]*strideJt + iy ;
	                        jiz = CBFPOWZ[j+jCBFInd]*strideJt + iz ;
                                tij = ti * jBasis->shells[jShell].primitives[jP].ccbf[j] ;
                                for ( k = 0 ; k < nCFuncK ; k++ )
                                {
                                    kjix = CBFPOWX[k+kCBFInd]*strideKt + jix ;
                                    kjiy = CBFPOWY[k+kCBFInd]*strideKt + jiy ;
                                    kjiz = CBFPOWZ[k+kCBFInd]*strideKt + jiz ;
                                    tijk = tij * kBasis->shells[kShell].primitives[kP].ccbf[k] ;
                                    for ( l = 0 ; l < nCFuncL ; l++, n++ )
                                    {
                                        lkjix = CBFPOWX[l+lCBFInd]*strideLt + kjix ;
                                        lkjiy = CBFPOWY[l+lCBFInd]*strideLt + kjiy ;
                                        lkjiz = CBFPOWZ[l+lCBFInd]*strideLt + kjiz ;
                                        for ( m = 0, fac = 0.0e+00 ; m < nRoots ; m++ ) fac += xInt[lkjix+m*strideM] * yInt[lkjiy+m*strideM] * zInt[lkjiz+m*strideM] ;
                                        tijkl = tijk * lBasis->shells[lShell].primitives[lP].ccbf[l] ;
                                        g[n] += tijkl * fac ;
                                    }
                                }
                            }
                        }
                    }
                } /* . lP. */
            } /* . kP. */
        } /* . jP. */
    } /* . iP. */
    /*----------------------------------------------------------------------------------------------------------------------------*/
                    /* . Save the integrals. */
                    {
                        auto Boolean     skip, skipIJ ;
                        auto Integer     i, ii, ij, j, jj, k, kk, kl, l, ll, m, n ;
                        auto Cardinal16 *indices16 = block->indices16 ;
                        auto Real       *integrals = block->data      ;
                        m = block->count ;
                        for ( i = 0, ij = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
                        {
                            ii = iBasis->shells[iShell].nstartw + i ;
                            for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; ij++, j++ )
                            {
                                jj     = jBasis->shells[jShell].nstartw + j ;
                                skipIJ = iAndJ && ( j > i ) ; /* . Do not skip here as need full loops for n counter. */
                                for ( k = 0, kl = 0 ; k < kBasis->shells[kShell].nbasisw ; k++ )
                                {
                                    kk = kBasis->shells[kShell].nstartw + k ;
                                    for ( l = 0 ; l < lBasis->shells[lShell].nbasisw ; kl++, l++, n++ )
                                    {
                                        ll   = lBasis->shells[lShell].nstartw + l ;
                                        skip = skipIJ || ( ijAndKL && ( ij < kl ) ) || ( kAndL && ( l > k ) ) ;
# ifdef _PrintIntegrals
ntotal += 1 ;
printf ( "%10d %6d %6d %6d %6d %25.15f | %6u %6u | %6d %6d %6d %6d | %6d %6d %6d %6d\n", ntotal, ii, jj, kk, ll, g[n], skipIJ, skip, i, j, k, l, iShell, jShell, kShell, lShell ) ;
fflush ( stdout ) ;
# endif
                                        if ( ! skip )
                                        {
                                            auto Integer  m4 = 4 * m ;
                                            indices16[m4  ] = ii ;
                                            indices16[m4+1] = jj ;
                                            indices16[m4+2] = kk ;
                                            indices16[m4+3] = ll ;
                                            integrals[m]    = g[n] ;
                                            m++ ;
                                        }
                                    }
                                }
                            }
                        }
                        block->count = m ;
                    }
                } /* . lShell. */
            } /* . kShell. */
        } /* . jShell. */
    } /* . iShell. */
}
# undef _IntegralSize

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the two-electron integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _DensityTolerance 1.0e-12
# define _IntegralSize0    ( ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP2 ) * MAXAMP2 * ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP2 ) * MAXAMP1 * _MAXRYS ) /* NI, NJ and NK increased by 1. */
# define _IntegralSize1    ( MAXAMP2 * MAXAMP2 * MAXAMP2 * MAXAMP1 * _MAXRYS )
# define _IntegralSizeS    ( MAXCBF*MAXCBF*MAXCBF*MAXCBF )
void GaussianBasisIntegrals_TEIsD ( const GaussianBasis *iBasis     ,
                                    const Real          *rI         ,
                                    const GaussianBasis *jBasis     ,
                                    const Real          *rJ         ,
                                    const Real          *rIJ        ,
                                    const Real           rIJ2       ,
                                    const GaussianBasis *kBasis     ,
                                    const Real          *rK         ,
                                    const GaussianBasis *lBasis     ,
                                    const Real          *rL         ,
                                    const Real          *rKL        ,
                                    const Real           rKL2       ,
                                    const Boolean        jLessThanL ,
                                          Block         *block      )
{
    /* . Pair data. */
    Boolean        iIsJ, iIsK, jIsL, kIsL ;
    /* . Shell loops. */
    Boolean        iAndJ  , ijAndKL, kAndL, qKL1 ;
    Integer        iAMMax , iAMMaxt, iCBFInd, iShell,         nCFuncI ,
                   jAMMax , jAMMaxt, jCBFInd, jShell, jUpper, nCFuncJ ,
                   kAMMax , kAMMaxt, kCBFInd, kShell, kUpper, nCFuncK ,
                   lAMMax , lAMMaxt, lCBFInd, lShell, lUpper, nCFuncL ,
                   mAMMax , nAMMax , nRoots ,
                   dStrideI, dStrideJ, dStrideK, dStrideL, dStrideM,
                   strideI, strideIt, strideJ, strideJt,
                   strideK, strideKt, strideL, strideLt, strideM ;
    Real           xIJt, xKLt, yIJt, yKLt, zIJt, zKLt ;
    Real           gIx[_IntegralSizeS], gIy[_IntegralSizeS], gIz[_IntegralSizeS] ,
                   gJx[_IntegralSizeS], gJy[_IntegralSizeS], gJz[_IntegralSizeS] ,
                   gKx[_IntegralSizeS], gKy[_IntegralSizeS], gKz[_IntegralSizeS] ;
    const Real    *rC, *rD ;
    /* . Primitive loops. */
    Integer        iP, jP, kP, lP ;
    Real           aa, aaInv, aAndB, ab, aI, aJ, aK, aL, ar2I, ar2K, arg, argIJ, bb, bbInv,
                   axac, ayac, azac, axad, ayad, azad, bxbc, bybc, bzbc, bxbd, bybd, bzbd,
                   c1x , c2x , c3x , c4x , c1y , c2y , c3y , c4y , c1z , c2z , c3z , c4z,
                   expFac, rho, xAB , yAB , zAB ;
    Real           arI[3], arK[3], rA[3], rB[3],
                   xDI [_IntegralSize1], yDI [_IntegralSize1], zDI [_IntegralSize1] ,
                   xDJ [_IntegralSize1], yDJ [_IntegralSize1], zDJ [_IntegralSize1] ,
                   xDK [_IntegralSize1], yDK [_IntegralSize1], zDK [_IntegralSize1] ,
                   xInt[_IntegralSize0], yInt[_IntegralSize0], zInt[_IntegralSize0] ;
    RysQuadrature  roots ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    iIsK = ( iBasis == kBasis ) && ( rI == rK ) ;
    jIsL = ( jBasis == lBasis ) && ( rJ == rL ) ;
    kIsL = ( kBasis == lBasis ) && ( rK == rL ) ;
    /* . Quadruple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        iAMMax  = iBasis->shells[iShell].type->angularmomentum_high ;
        iCBFInd = iBasis->shells[iShell].type->cbfindex ;
        nCFuncI = iBasis->shells[iShell].type->ncbf     ;
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nshells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jAMMax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jCBFInd = jBasis->shells[jShell].type->cbfindex ;
            nCFuncJ = jBasis->shells[jShell].type->ncbf     ;
            nAMMax  = iAMMax + jAMMax + 2 ; /* . I and J increased by 1. */
            if ( iAMMax >= jAMMax )
            {
               iAMMaxt = iAMMax + 1 ; /* . AMMax normal, AMMaxt increased by 1. */
               jAMMaxt = jAMMax + 1 ;
               xIJt    = rIJ[0] ;
               yIJt    = rIJ[1] ;
               zIJt    = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iAMMaxt = jAMMax + 1 ;
               jAMMaxt = iAMMax + 1 ;
               xIJt    = - rIJ[0] ;
               yIJt    = - rIJ[1] ;
               zIJt    = - rIJ[2] ;
               rC      = rJ ;
            }
            iAndJ = iIsJ && ( iShell == jShell ) ;
            if ( iIsK ) kUpper = iShell + 1 ;
            else        kUpper = kBasis->nshells ;
            for ( kShell = 0 ; kShell < kUpper ; kShell++ )
            {
                kAMMax  = kBasis->shells[kShell].type->angularmomentum_high ;
                kCBFInd = kBasis->shells[kShell].type->cbfindex ;
                nCFuncK = kBasis->shells[kShell].type->ncbf     ;
                if ( iIsK && ( iShell == kShell ) )
                {
                          if ( jIsL       ) { lUpper = jShell + 1      ; }
                     else if ( jLessThanL ) { lUpper = -1              ; }
                     else                   { lUpper = lBasis->nshells ; }
                 }
                else
                {    if ( kIsL ) { lUpper = kShell + 1      ; }
                     else        { lUpper = lBasis->nshells ; }
                }
                for ( lShell = 0 ; lShell < lUpper ; lShell++ )
                {
                    auto Integer  i ;
                    lAMMax  = lBasis->shells[lShell].type->angularmomentum_high ;
                    lCBFInd = lBasis->shells[lShell].type->cbfindex ;
                    nCFuncL = lBasis->shells[lShell].type->ncbf     ;
                    mAMMax  = kAMMax + lAMMax + 1 ; /* . K increased by 1. */
                    if ( kAMMax + 1 >= lAMMax )
                    {
                       kAMMaxt = kAMMax + 1 ; /* . AMMax normal, AMMaxt increased by 1. */
                       lAMMaxt = lAMMax ;
                       xKLt    = rKL[0] ;
                       yKLt    = rKL[1] ;
                       zKLt    = rKL[2] ;
                       rD      = rK ;
                    }
                    else
                    {
                       kAMMaxt = lAMMax ;
                       lAMMaxt = kAMMax + 1 ;
                       xKLt    = - rKL[0] ;
                       yKLt    = - rKL[1] ;
                       zKLt    = - rKL[2] ;
                       rD      = rL ;
                    }
                    kAndL    = kIsL && ( kShell == lShell ) ;
                    ijAndKL  = iIsK && ( iShell == kShell ) && jIsL && ( jShell == lShell ) ;
                    qKL1     = ( mAMMax <= 1 ) ;
                    /* . Integral initialization. */
                    dStrideL = 1                         ; strideL = 1                         ;
                    dStrideK = ( lAMMax + 1 ) * dStrideL ; strideK = ( lAMMaxt + 1 ) * strideL ;
                    dStrideJ = ( kAMMax + 1 ) * dStrideK ; strideJ = ( mAMMax  + 1 ) * strideK ;
                    dStrideI = ( jAMMax + 1 ) * dStrideJ ; strideI = ( jAMMaxt + 1 ) * strideJ ;
                    dStrideM = ( iAMMax + 1 ) * dStrideI ; strideM = ( nAMMax  + 1 ) * strideI ;
                    if ( iAMMax     >= jAMMax ) { strideIt = strideI ; strideJt = strideJ ; }
                    else                        { strideIt = strideJ ; strideJt = strideI ; }
                    if ( kAMMax + 1 >= lAMMax ) { strideKt = strideK ; strideLt = strideL ; }
                    else                        { strideKt = strideL ; strideLt = strideK ; }
                    nRoots  = ( mAMMax + nAMMax ) / 2 + 1 ;
                    for ( i = 0 ; i < ( nCFuncI * nCFuncJ * nCFuncK * nCFuncL ) ; i++ )
                    {
                        gIx[i] = gIy[i] = gIz[i] = gJx[i] = gJy[i] = gJz[i] = gKx[i] = gKy[i] = gKz[i] = 0.0e+00 ;
                    }
    /*------------------------------------------------------------------------------------------------------------------------------
    ! . Quadruple loop over primitives.
    !-----------------------------------------------------------------------------------------------------------------------------*/
    for ( iP = 0 ; iP < iBasis->shells[iShell].nprimitives ; iP++ )
    {
        aI   = iBasis->shells[iShell].primitives[iP].exponent ;
        ar2I = aI * rIJ2 ;
        for ( i = 0 ; i < 3 ; i++ ) arI[i] = aI * rI[i] ;
        for ( jP = 0 ; jP < jBasis->shells[jShell].nprimitives ; jP++ )
        {
            aJ    = jBasis->shells[jShell].primitives[jP].exponent ;
            aa    = aI + aJ ;
            aaInv = 1.0e+00 / aa ;
            argIJ = aJ * ar2I * aaInv ;
            if ( argIJ > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
            for ( i = 0 ; i < 3 ; i++ ) rA[i] = ( arI[i] + aJ * rJ[i] ) * aaInv ;
            axad = aa * ( rA[0] - rD[0] ) ;
            ayad = aa * ( rA[1] - rD[1] ) ;
            azad = aa * ( rA[2] - rD[2] ) ;
            axac = aa * ( rA[0] - rC[0] ) ;
            ayac = aa * ( rA[1] - rC[1] ) ;
            azac = aa * ( rA[2] - rC[2] ) ;
            for ( kP = 0 ; kP < kBasis->shells[kShell].nprimitives ; kP++ )
            {
                aK   = kBasis->shells[kShell].primitives[kP].exponent ;
                ar2K = aK * rKL2 ;
                for ( i = 0 ; i < 3 ; i++ ) arK[i] = aK * rK[i] ;
                for ( lP = 0 ; lP < lBasis->shells[lShell].nprimitives ; lP++ )
                {
                    auto Integer  m ;
                    aL     = lBasis->shells[lShell].primitives[lP].exponent ;
                    bb     = aK + aL ;
                    bbInv  = 1.0e+00 / bb ;
                    arg    = argIJ + aL * ar2K * bbInv ;
                    if ( arg > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    ab     = aa * bb ;
                    aAndB  = aa + bb ;
                    rho    = ab / aAndB ;
                    expFac = exp ( - arg ) * PI252 / ( ab * sqrt ( aAndB ) ) ;
                    for ( i = 0 ; i < 3 ; i++ ) rB[i] = ( arK[i] + aL * rL[i] ) * bbInv ;
                    bxbd = bb * ( rB[0] - rD[0] ) ;
                    bybd = bb * ( rB[1] - rD[1] ) ;
                    bzbd = bb * ( rB[2] - rD[2] ) ;
                    bxbc = bb * ( rB[0] - rC[0] ) ;
                    bybc = bb * ( rB[1] - rC[1] ) ;
                    bzbc = bb * ( rB[2] - rC[2] ) ;
                    c1x  = bxbd + axad ;
                    c2x  = aa * bxbd   ;
                    c3x  = bxbc + axac ;
                    c4x  = bb * axac   ;
                    c1y  = bybd + ayad ;
                    c2y  = aa * bybd   ;
                    c3y  = bybc + ayac ;
                    c4y  = bb * ayac   ;
                    c1z  = bzbd + azad ;
                    c2z  = aa * bzbd   ;
                    c3z  = bzbc + azac ;
                    c4z  = bb * azac   ;
                    xAB  = rA[0] - rB[0] ;
                    yAB  = rA[1] - rB[1] ;
                    zAB  = rA[2] - rB[2] ;
                    RysQuadrature_Roots ( &roots, nRoots, rho * ( xAB*xAB + yAB*yAB + zAB*zAB ) ) ;
                    for ( m = 0 ; m < nRoots ; m++ )
                    {
                        auto Real b00, b10, bp01, f00, fac, fac2, u2, xc00, xcp00, yc00, ycp00, zc00, zcp00 ;
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] * expFac ;
                        fac   = 1.0e+00 / ( ab + u2 * aAndB ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( aa   + u2 ) * fac2 ;
                        b00   =          u2   * fac2 ;
                        b10   = ( bb   + u2 ) * fac2 ;
                        xcp00 = ( u2 * c1x + c2x ) * fac ;
                        ycp00 = ( u2 * c1y + c2y ) * fac ;
                        zcp00 = ( u2 * c1z + c2z ) * fac ;
                        xc00  = ( u2 * c3x + c4x ) * fac ;
                        yc00  = ( u2 * c3y + c4y ) * fac ;
                        zc00  = ( u2 * c3z + c4z ) * fac ;
                        Subsidiary_Integral_Nuclear4C ( iAMMaxt ,
                                                        jAMMaxt ,
                                                        nAMMax  ,
                                                        kAMMaxt ,
                                                        lAMMaxt ,
                                                        mAMMax  ,
                                                        False   ,
                                                        False   ,
                                                        False   ,
                                                        qKL1    ,
                                                        b00     ,
                                                        b10     ,
                                                        bp01    ,
                                                        xIJt    ,
                                                        yIJt    ,
                                                        zIJt    ,
                                                        xKLt    ,
                                                        yKLt    ,
                                                        zKLt    ,
                                                        f00     ,
                                                        xc00    ,
                                                        xcp00   ,
                                                        yc00    ,
                                                        ycp00   ,
                                                        zc00    ,
                                                        zcp00   ,
                                                        strideI ,
                                                        strideJ ,
                                                        strideK ,
                                                        &xInt[m*strideM] ,
                                                        &yInt[m*strideM] ,
                                                        &zInt[m*strideM] ) ;
                        Subsidiary_Integral_Derivative4 ( iAMMax   ,
                                                          jAMMax   ,
                                                          kAMMax   ,
                                                          lAMMax   ,
                                                          strideIt ,
                                                          strideJt ,
                                                          strideKt ,
                                                          strideLt ,
                                                          dStrideI ,
                                                          dStrideJ ,
                                                          dStrideK ,
                                                          dStrideL ,
                                                          aI       ,
                                                          aJ       ,
                                                          aK       ,
                                                          &xInt[m* strideM] ,
                                                          &yInt[m* strideM] ,
                                                          &zInt[m* strideM] ,
                                                          &xDI [m*dStrideM] ,
                                                          &yDI [m*dStrideM] ,
                                                          &zDI [m*dStrideM] ,
                                                          &xDJ [m*dStrideM] ,
                                                          &yDJ [m*dStrideM] ,
                                                          &zDJ [m*dStrideM] ,
                                                          &xDK [m*dStrideM] ,
                                                          &yDK [m*dStrideM] ,
                                                          &zDK [m*dStrideM] ) ;
                    }
                    /* . Assemble the integrals. */
                    {
                        auto Integer  i, ix, ixd, iy, iyd, iz, izd, j, jix, jixd, jiy, jiyd, jiz, jizd,
                                     k, kjix, kjixd, kjiy, kjiyd, kjiz, kjizd, l, lkjix, lkjixd, lkjiy, lkjiyd, lkjiz, lkjizd, m, n ;
                        auto Real    facIx, facIy, facIz, facJx, facJy, facJz, facKx, facKy, facKz, ti, tij, tijk, tijkl ;
                        for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                        {
                            ix  = CBFPOWX[i+iCBFInd]*strideIt ;
                            iy  = CBFPOWY[i+iCBFInd]*strideIt ;
                            iz  = CBFPOWZ[i+iCBFInd]*strideIt ;
                            ixd = CBFPOWX[i+iCBFInd]*dStrideI ;
                            iyd = CBFPOWY[i+iCBFInd]*dStrideI ;
                            izd = CBFPOWZ[i+iCBFInd]*dStrideI ;
                            ti  = iBasis->shells[iShell].primitives[iP].ccbf[i] ;
                            for ( j = 0 ; j < nCFuncJ ; j++ )
                            {
                                jix  = CBFPOWX[j+jCBFInd]*strideJt + ix  ;
                                jiy  = CBFPOWY[j+jCBFInd]*strideJt + iy  ;
                                jiz  = CBFPOWZ[j+jCBFInd]*strideJt + iz  ;
                                jixd = CBFPOWX[j+jCBFInd]*dStrideJ + ixd ;
                                jiyd = CBFPOWY[j+jCBFInd]*dStrideJ + iyd ;
                                jizd = CBFPOWZ[j+jCBFInd]*dStrideJ + izd ;
                                tij  = ti * jBasis->shells[jShell].primitives[jP].ccbf[j] ;
                                for ( k = 0 ; k < nCFuncK ; k++ )
                                {
                                    kjix  = CBFPOWX[k+kCBFInd]*strideKt + jix  ;
                                    kjiy  = CBFPOWY[k+kCBFInd]*strideKt + jiy  ;
                                    kjiz  = CBFPOWZ[k+kCBFInd]*strideKt + jiz  ;
                                    kjixd = CBFPOWX[k+kCBFInd]*dStrideK + jixd ;
                                    kjiyd = CBFPOWY[k+kCBFInd]*dStrideK + jiyd ;
                                    kjizd = CBFPOWZ[k+kCBFInd]*dStrideK + jizd ;
                                    tijk  = tij * kBasis->shells[kShell].primitives[kP].ccbf[k] ;
                                    for ( l = 0 ; l < nCFuncL ; l++, n++ )
                                    {
                                        lkjix  = CBFPOWX[l+lCBFInd]*strideLt + kjix  ;
                                        lkjiy  = CBFPOWY[l+lCBFInd]*strideLt + kjiy  ;
                                        lkjiz  = CBFPOWZ[l+lCBFInd]*strideLt + kjiz  ;
                                        lkjixd = CBFPOWX[l+lCBFInd]*dStrideL + kjixd ;
                                        lkjiyd = CBFPOWY[l+lCBFInd]*dStrideL + kjiyd ;
                                        lkjizd = CBFPOWZ[l+lCBFInd]*dStrideL + kjizd ;
                                        for ( m = 0, facIx = facIy = facIz = facJx = facJy = facJz = facKx = facKy = facKz = 0.0e+00 ; m < nRoots ; m++ )
                                        {
                                            facIx += xDI [lkjixd+m*dStrideM] * yInt[lkjiy+ m* strideM] * zInt[lkjiz+ m* strideM] ;
                                            facIy += xInt[lkjix+ m* strideM] * yDI [lkjiyd+m*dStrideM] * zInt[lkjiz+ m* strideM] ;
                                            facIz += xInt[lkjix+ m* strideM] * yInt[lkjiy+ m* strideM] * zDI [lkjizd+m*dStrideM] ;
                                            facJx += xDJ [lkjixd+m*dStrideM] * yInt[lkjiy+ m* strideM] * zInt[lkjiz+ m* strideM] ;
                                            facJy += xInt[lkjix+ m* strideM] * yDJ [lkjiyd+m*dStrideM] * zInt[lkjiz+ m* strideM] ;
                                            facJz += xInt[lkjix+ m* strideM] * yInt[lkjiy+ m* strideM] * zDJ [lkjizd+m*dStrideM] ;
                                            facKx += xDK [lkjixd+m*dStrideM] * yInt[lkjiy+ m* strideM] * zInt[lkjiz+ m* strideM] ;
                                            facKy += xInt[lkjix+ m* strideM] * yDK [lkjiyd+m*dStrideM] * zInt[lkjiz+ m* strideM] ;
                                            facKz += xInt[lkjix+ m* strideM] * yInt[lkjiy+ m* strideM] * zDK [lkjizd+m*dStrideM] ;
                                        }
                                        tijkl   = tijk * lBasis->shells[lShell].primitives[lP].ccbf[l] ;
                                        gIx[n] += tijkl * facIx ;
                                        gIy[n] += tijkl * facIy ;
                                        gIz[n] += tijkl * facIz ;
                                        gJx[n] += tijkl * facJx ;
                                        gJy[n] += tijkl * facJy ;
                                        gJz[n] += tijkl * facJz ;
                                        gKx[n] += tijkl * facKx ;
                                        gKy[n] += tijkl * facKy ;
                                        gKz[n] += tijkl * facKz ;
                                    }
                                }
                            }
                        }
                    }
                } /* . lP. */
            } /* . kP. */
        } /* . jP. */
    } /* . iP. */
    /*----------------------------------------------------------------------------------------------------------------------------*/
                    /* . Save the integrals. */
                    {
                        auto Boolean     skip, skipIJ ;
                        auto Integer     i, ii, ij, j, jj, k, kk, kl, l, ll, m, n ;
                        auto Cardinal16 *indices16 = block->indices16 ;
                        auto Real       *integrals = block->data      ;
                        m = block->count ;
                        for ( i = 0, ij = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
                        {
                            ii = iBasis->shells[iShell].nstartw + i ;
                            for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; ij++, j++ )
                            {
                                jj     = jBasis->shells[jShell].nstartw + j ;
                                skipIJ = iAndJ && ( j > i ) ; /* . Do not skip here as need full loops for n counter. */
                                for ( k = 0, kl = 0 ; k < kBasis->shells[kShell].nbasisw ; k++ )
                                {
                                    kk = kBasis->shells[kShell].nstartw + k ;
                                    for ( l = 0 ; l < lBasis->shells[lShell].nbasisw ; kl++, l++, n++ )
                                    {
                                        ll   = lBasis->shells[lShell].nstartw + l ;
                                        skip = skipIJ || ( ijAndKL && ( ij < kl ) ) || ( kAndL && ( l > k ) ) ;
                                        if ( ! skip )
                                        {
                                            auto Integer  m4 = 4 * m, m9 = 9 * m ;
                                            indices16[m4  ] = ii ;
                                            indices16[m4+1] = jj ;
                                            indices16[m4+2] = kk ;
                                            indices16[m4+3] = ll ;
                                            integrals[m9  ] = gIx[n] ;
                                            integrals[m9+1] = gIy[n] ;
                                            integrals[m9+2] = gIz[n] ;
                                            integrals[m9+3] = gJx[n] ;
                                            integrals[m9+4] = gJy[n] ;
                                            integrals[m9+5] = gJz[n] ;
                                            integrals[m9+6] = gKx[n] ;
                                            integrals[m9+7] = gKy[n] ;
                                            integrals[m9+8] = gKz[n] ;
                                            m++ ;
                                        }
                                    }
                                }
                            }
                        }
                        block->count = m ;
                    }
                } /* . lShell. */
            } /* . kShell. */
        } /* . jShell. */
    } /* . iShell. */
}
# undef _DensityTolerance
# undef _IntegralSize0
# undef _IntegralSize1
# undef _IntegralSizeS
