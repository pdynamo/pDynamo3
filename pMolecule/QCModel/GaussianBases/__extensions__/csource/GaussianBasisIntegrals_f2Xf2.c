/*==================================================================================================================================
! . Integrals - 4 basis, 2 electrons, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "Boolean.h"
# include "GaussianBasisIntegrals_f2Xf2.h"
# include "GaussianBasisSubsidiary.h"
# include "GaussianBasisTransform.h"
# include "Integer.h"
# include "NumericalMacros.h"
# include "Real.h"
# include "RysQuadrature.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the anti-Coulomb two-electron integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s4 and real 3 * s4 where s4 = ( maximum shell size )^4. */
# define MAXAMP21 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP1 )
# define MAXAMP23 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP3 )
void GaussianBasisIntegrals_f2Af2i ( const GaussianBasis *iBasis     ,
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
                                     const Integer        s4         ,
                                           Integer       *iWork      ,
                                           Real          *rWork      ,
                                           Block         *block      )
{
    /* . Pair data. */
    Boolean        iIsJ, iIsK, jIsL, kIsL ;
    /* . Shell loops. */
    Boolean        iAndJ  , ijAndKL, kAndL ;
    Integer        iAMMax , iAMMaxT, iShell,         nCFuncI ,
                   jAMMax , jAMMaxT, jShell, jUpper, nCFuncJ ,
                   kAMMax , kAMMaxT, kShell, kUpper, nCFuncK ,
                   lAMMax , lAMMaxT, lShell, lUpper, nCFuncL ,
                   mAMMax , nAMMax , nRoots ,
                   uStrideI, uStrideIt, uStrideJ, uStrideJt , uStrideKL           , uStrideM ,
                   vStrideI, vStrideJ , vStrideK, vStrideKt , vStrideL , vStrideLt, vStrideM ;
    Real           xIJt, xKLt, yIJt, yKLt, zIJt, zKLt ;
    Real          *pG, *values = NULL, *work = NULL ;
    Real          *g, *gT ;
    const Real    *rC, *rD ;
    /* . Primitive loops. */
    Integer        iP, jP, kP, lP ;
    Real           aa, aaInv, aAndB, ab, aI, aJ, aK, aL, ar2I, ar2K, arg, argIJ, bb, bbInv,
                   axac, ayac, azac, axad, ayad, azad, bxbc, bybc, bzbc, bxbd, bybd, bzbd,
                   c1x , c2x , c3x , c4x , c1y , c2y , c3y , c4y , c1z , c2z , c3z , c4z,
                   expFac, rho, xAB , xCD , yAB , yCD , zAB , zCD ;
    Real           arI[3], arK[3], rA[3], rB[3],
                   Gx[MAXAMP23*MAXAMP23              ] ,
                   Gy[MAXAMP23*MAXAMP23              ] ,
                   Gz[MAXAMP23*MAXAMP23              ] ,
                   Hx[MAXAMP21*MAXAMP21              ] ,
                   Hy[MAXAMP21*MAXAMP21              ] ,
                   Hz[MAXAMP21*MAXAMP21              ] ,
                   Sx[MAXAMP21*MAXAMP1*MAXAMP1       ] ,
                   Sy[MAXAMP21*MAXAMP1*MAXAMP1       ] ,
                   Sz[MAXAMP21*MAXAMP1*MAXAMP1       ] ,
                   Tx[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ,
                   Ty[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ,
                   Tz[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ,
                   Ux[MAXAMP21*MAXAMP1*MAXAMP1       ] ,
                   Uy[MAXAMP21*MAXAMP1*MAXAMP1       ] ,
                   Uz[MAXAMP21*MAXAMP1*MAXAMP1       ] ,
                   Vx[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ,
                   Vy[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ,
                   Vz[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ;
    /* . Function loops. */
    Integer        f, i, ix, iy, iz, j, jix, jiy, jiz, k, kjix, kjiy, kjiz, l, m, nCFunc ;
    Integer       *Ix, *Iy, *Iz  ;
    Real           tI, tIJ, tIJK ;
    Real          *Cijkl ;
    RysQuadrature  roots ;
# ifdef _PrintIntegrals_
auto Integer ntotal = 0 ;
# endif
    /* . Initialization. */
    block->count = 0 ;
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    iIsK = ( iBasis == kBasis ) && ( rI == rK ) ;
    jIsL = ( jBasis == lBasis ) && ( rJ == rL ) ;
    kIsL = ( kBasis == lBasis ) && ( rK == rL ) ;
    /* . Set pointers. */
    Cijkl = &rWork[0] ; g  = &rWork[s4] ; gT = &rWork[2*s4] ;
    Ix    = &iWork[0] ; Iy = &iWork[s4] ; Iz = &iWork[2*s4] ;
    /* . Quadruple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iAMMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jAMMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            nAMMax  = iAMMax + jAMMax ;
            if ( iAMMax >= jAMMax )
            {
               iAMMaxT = iAMMax ;
               jAMMaxT = jAMMax ;
               xIJt    = rIJ[0] ;
               yIJt    = rIJ[1] ;
               zIJt    = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iAMMaxT = jAMMax ;
               jAMMaxT = iAMMax ;
               xIJt    = - rIJ[0] ;
               yIJt    = - rIJ[1] ;
               zIJt    = - rIJ[2] ;
               rC      = rJ ;
            }
            iAndJ = iIsJ && ( iShell == jShell ) ;
            if ( iIsK ) kUpper = iShell + 1 ;
            else        kUpper = kBasis->nShells ;
            for ( kShell = 0 ; kShell < kUpper ; kShell++ )
            {
                kAMMax  = kBasis->shells[kShell].lHigh ;
                nCFuncK = kBasis->shells[kShell].nCBF     ;
                if ( iIsK && ( iShell == kShell ) )
                {
                          if ( jIsL       ) { lUpper = jShell + 1      ; }
                     else if ( jLessThanL ) { lUpper = -1              ; }
                     else                   { lUpper = lBasis->nShells ; }
                }
                else
                {    if ( kIsL ) { lUpper = kShell + 1      ; }
                     else        { lUpper = lBasis->nShells ; }
                }
                for ( lShell = 0 ; lShell < lUpper ; lShell++ )
                {
                    lAMMax  = lBasis->shells[lShell].lHigh ;
                    nCFuncL = lBasis->shells[lShell].nCBF     ;
                    mAMMax  = kAMMax + lAMMax ;
                    kAndL   = kIsL && ( kShell == lShell ) ;
                    ijAndKL = iIsK && ( iShell == kShell ) && jIsL && ( jShell == lShell ) ;
                    if ( kAMMax >= lAMMax )
                    {
                       kAMMaxT = kAMMax ;
                       lAMMaxT = lAMMax ;
                       xKLt    = rKL[0] ;
                       yKLt    = rKL[1] ;
                       zKLt    = rKL[2] ;
                       rD      = rK ;
                    }
                    else
                    {
                       kAMMaxT = lAMMax ;
                       lAMMaxT = kAMMax ;
                       xKLt    = - rKL[0] ;
                       yKLt    = - rKL[1] ;
                       zKLt    = - rKL[2] ;
                       rD      = rL ;
                    }
                    /* . Displacement. */
                    xCD = ( rC[0] - rD[0] ) ;
                    yCD = ( rC[1] - rD[1] ) ;
                    zCD = ( rC[2] - rD[2] ) ;
                    /* . Roots. */
                    nRoots  = ( mAMMax + nAMMax + 4 ) / 2 + 1 ; /* . +4 */
                    /* . Strides. */
                    uStrideKL = 1 ;
                    uStrideJ  = ( mAMMax + 1 ) * uStrideKL ;
                    uStrideI  = ( jAMMax + 1 ) * uStrideJ  ;
                    uStrideM  = ( iAMMax + 1 ) * uStrideI  ;
                    vStrideL  = 1 ;
                    vStrideK  = ( lAMMax + 1 ) * vStrideL  ;
                    vStrideJ  = ( kAMMax + 1 ) * vStrideK  ;
                    vStrideI  = ( jAMMax + 1 ) * vStrideJ  ;
                    vStrideM  = ( iAMMax + 1 ) * vStrideI  ;
                    if ( iAMMax >= jAMMax ) { uStrideIt = uStrideI ; uStrideJt = uStrideJ ; }
                    else                    { uStrideIt = uStrideJ ; uStrideJt = uStrideI ; }
                    if ( kAMMax >= lAMMax ) { vStrideKt = vStrideK ; vStrideLt = vStrideL ; }
                    else                    { vStrideKt = vStrideL ; vStrideLt = vStrideK ; }
                    /* . Index arrays. */
                    for ( i = 0, f = 0 ; i < nCFuncI ; i++ )
                    {
   	                ix = iBasis->shells[iShell].cbfPowX[i]*vStrideI ;
	                iy = iBasis->shells[iShell].cbfPowY[i]*vStrideI ;
	                iz = iBasis->shells[iShell].cbfPowZ[i]*vStrideI ;
                        for ( j = 0 ; j < nCFuncJ ; j++ )
                        {
	                    jix = jBasis->shells[jShell].cbfPowX[j]*vStrideJ + ix ;
	                    jiy = jBasis->shells[jShell].cbfPowY[j]*vStrideJ + iy ;
	                    jiz = jBasis->shells[jShell].cbfPowZ[j]*vStrideJ + iz ;
                            for ( k = 0 ; k < nCFuncK ; k++ )
                            {
                                kjix = kBasis->shells[kShell].cbfPowX[k]*vStrideK + jix ;
                                kjiy = kBasis->shells[kShell].cbfPowY[k]*vStrideK + jiy ;
                                kjiz = kBasis->shells[kShell].cbfPowZ[k]*vStrideK + jiz ;
                                for ( l = 0 ; l < nCFuncL ; l++, f++ )
                                {
                                    Ix[f] = lBasis->shells[lShell].cbfPowX[l] + kjix ;
                                    Iy[f] = lBasis->shells[lShell].cbfPowY[l] + kjiy ;
                                    Iz[f] = lBasis->shells[lShell].cbfPowZ[l] + kjiz ;
                                }
                            }
                        }
                    }
                    /* . Function initialization. */
                    nCFunc = nCFuncI * nCFuncJ * nCFuncK * nCFuncL ;
                    for ( f = 0 ; f < nCFunc ; f++ ) g[f] = 0.0e+00 ;
    /*------------------------------------------------------------------------------------------------------------------------------
    ! . Quadruple loop over primitives.
    !-----------------------------------------------------------------------------------------------------------------------------*/
    for ( iP = 0 ; iP < iBasis->shells[iShell].nPrimitives ; iP++ )
    {
	aI   = iBasis->shells[iShell].primitives[iP].exponent ;
	ar2I = aI * rIJ2 ;
	for ( i = 0 ; i < 3 ; i++ ) arI[i] = aI * rI[i] ;
        for ( jP = 0 ; jP < jBasis->shells[jShell].nPrimitives ; jP++ )
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
            for ( kP = 0 ; kP < kBasis->shells[kShell].nPrimitives ; kP++ )
            {
	        aK   = kBasis->shells[kShell].primitives[kP].exponent ;
	        ar2K = aK * rKL2 ;
	        for ( i = 0 ; i < 3 ; i++ ) arK[i] = aK * rK[i] ;
                for ( lP = 0 ; lP < lBasis->shells[lShell].nPrimitives ; lP++ )
                {
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
                    /* . Get coefficient array. */
                    for ( i = 0, f = 0 ; i < nCFuncI ; i++ )
                    {
                        tI = iBasis->shells[iShell].primitives[iP].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++ )
                        {
                            tIJ = tI * jBasis->shells[jShell].primitives[jP].cCBF[j] ;
                            for ( k = 0 ; k < nCFuncK ; k++ )
                            {
                                tIJK = tIJ * kBasis->shells[kShell].primitives[kP].cCBF[k] ;
                                for ( l = 0 ; l < nCFuncL ; l++, f++ ) Cijkl[f] = tIJK * lBasis->shells[lShell].primitives[lP].cCBF[l] ;
                            }
                        }
                    }
                    /* . Loop over Rys roots. */
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
                        GaussianBasisSubsidiary_f1Cg1  ( nAMMax+2  ,
                                                         mAMMax+2  ,
                                                         b00       ,
                                                         b10       ,
                                                         bp01      ,
                                                         f00       ,
                                                         xc00      ,
                                                         xcp00     ,
                                                         yc00      ,
                                                         ycp00     ,
                                                         zc00      ,
                                                         zcp00     ,
                                                         mAMMax+3  , /* . gStrideIJ. */
                                                         Gx        ,
                                                         Gy        ,
                                                         Gz        ) ;
                        GaussianBasisSubsidiary_f1Ag1  ( nAMMax    ,
                                                         mAMMax    ,
                                                         mAMMax+3  , /* . gStrideIJ. */
                                                         Gx        ,
                                                         Gy        ,
                                                         Gz        ,
                                                         xCD       ,
                                                         yCD       ,
                                                         zCD       ,
                                                         mAMMax+1  , /* . hStrideIJ. */
                                                         Hx        ,
                                                         Hy        ,
                                                         Hz        ) ;
                        /* . Reset S, T, U and V. */
                        for ( i = 0 ; i < uStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = Ux[i] = Uy[i] = Uz[i] = 0.0e+00 ;
                        for ( i = 0 ; i < vStrideM ; i++ ) Tx[i] = Ty[i] = Tz[i] = Vx[i] = Vy[i] = Vz[i] = 0.0e+00 ;
                        GaussianBasisSubsidiary_f1Xg2i ( iAMMaxT   ,
                                                         jAMMaxT   ,
                                                         mAMMax    ,
                                                         mAMMax+3  , /* . gStrideIJ. */
                                                         1         , /* . gStrideKL. */
                                                         Gx        ,
                                                         Gy        ,
                                                         Gz        ,
                                                         xIJt      ,
                                                         yIJt      ,
                                                         zIJt      ,
                                                         uStrideIt ,
                                                         uStrideJt ,
                                                         1         , /* . uStrideKL. */
                                                         Sx        ,
                                                         Sy        ,
                                                         Sz        ) ;
                        GaussianBasisSubsidiary_f1Xg2i ( kAMMaxT   ,
                                                         lAMMaxT   ,
                                                         (iAMMaxT+1)*(jAMMaxT+1)-1, /* . IJ loop range upper limit [0,u]. */
                                                         1         , /* . uStrideKL. */
                                                         uStrideJ  , /* . uStrideIJ. */
                                                         Sx        ,
                                                         Sy        ,
                                                         Sz        ,
                                                         xKLt      ,
                                                         yKLt      ,
                                                         zKLt      ,
                                                         vStrideKt ,
                                                         vStrideLt ,
                                                         vStrideJ  , /* . vStrideIJ. */
                                                         Tx        ,
                                                         Ty        ,
                                                         Tz        ) ;
                        GaussianBasisSubsidiary_f1Xg2i ( iAMMaxT   ,
                                                         jAMMaxT   ,
                                                         mAMMax    ,
                                                         mAMMax+1  , /* . hStrideIJ. */
                                                         1         , /* . hStrideKL. */
                                                         Hx        ,
                                                         Hy        ,
                                                         Hz        ,
                                                         xIJt      ,
                                                         yIJt      ,
                                                         zIJt      ,
                                                         uStrideIt ,
                                                         uStrideJt ,
                                                         1         , /* . uStrideKL. */
                                                         Ux        ,
                                                         Uy        ,
                                                         Uz        ) ;
                        GaussianBasisSubsidiary_f1Xg2i ( kAMMaxT   ,
                                                         lAMMaxT   ,
                                                         (iAMMaxT+1)*(jAMMaxT+1)-1, /* . IJ loop range upper limit [0,u]. */
                                                         1         , /* . uStrideKL. */
                                                         uStrideJ  , /* . uStrideIJ. */
                                                         Ux        ,
                                                         Uy        ,
                                                         Uz        ,
                                                         xKLt      ,
                                                         yKLt      ,
                                                         zKLt      ,
                                                         vStrideKt ,
                                                         vStrideLt ,
                                                         vStrideJ  , /* . vStrideIJ. */
                                                         Vx        ,
                                                         Vy        ,
                                                         Vz        ) ;
                        /* . Assemble the integrals. */
                        for ( f = 0 ; f < nCFunc ; f++ )
                        {
                            g[f] += Cijkl[f] * ( Vx[Ix[f]] * Ty[Iy[f]] * Tz[Iz[f]] + 
                                                 Tx[Ix[f]] * Vy[Iy[f]] * Tz[Iz[f]] + 
                                                 Tx[Ix[f]] * Ty[Iy[f]] * Vz[Iz[f]] ) ;
/* . Coulomb TEIs as test.
                            g[f] += ( Cijkl[f] * Tx[Ix[f]] * Ty[Iy[f]] * Tz[Iz[f]] ) ;
*/
                        }
                    } /* . nRoots. */
                } /* . lP. */
            } /* . kP. */
        } /* . jP. */
    } /* . iP. */
    /*----------------------------------------------------------------------------------------------------------------------------*/
    /* . Transform and save the integrals. */
    {
        auto Boolean     skip, skipIJ ;
        auto Integer     i, ii, ij, j, jj, k, kk, kl, l, ll, m, n ;
        auto Cardinal16 *indices16 = block->indices16 ;
        auto Real       *integrals = block->data      ;
        values = g  ;
        work   = gT ;
        GaussianBasisTransform4 ( nCFuncI                    ,
                                  nCFuncJ                    ,
                                  nCFuncK                    ,
                                  nCFuncL                    ,
                                  iBasis->shells[iShell].c2s ,
                                  jBasis->shells[jShell].c2s ,
                                  kBasis->shells[kShell].c2s ,
                                  lBasis->shells[lShell].c2s ,
                                  &values                    ,
                                  &work                      ) ;
        pG = values ;
        m  = block->count ;
        for ( i = 0, ij = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
        {
            ii = iBasis->shells[iShell].nStart + i ;
            for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; ij++, j++ )
            {
                jj     = jBasis->shells[jShell].nStart + j ;
                skipIJ = iAndJ && ( j > i ) ; /* . Do not skip here as need full loops for n counter. */
                for ( k = 0, kl = 0 ; k < kBasis->shells[kShell].nBasis ; k++ )
                {
                    kk = kBasis->shells[kShell].nStart + k ;
                    for ( l = 0 ; l < lBasis->shells[lShell].nBasis ; kl++, l++, n++ )
                    {
                        ll   = lBasis->shells[lShell].nStart + l ;
                        skip = skipIJ || ( ijAndKL && ( ij < kl ) ) || ( kAndL && ( l > k ) ) ;
# ifdef _PrintIntegrals_
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
                            integrals[m]    = -pG[n] ; /* . -r12 operator. */
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
# undef MAXAMP21
# undef MAXAMP23

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the Coulomb two-electron integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s4 and real 3 * s4 where s4 = ( maximum shell size )^4. */
# define MAXAMP21 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP1 )
void GaussianBasisIntegrals_f2Cf2i ( const GaussianBasis *iBasis     ,
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
                                     const Integer        s4         ,
                                           Integer       *iWork      ,
                                           Real          *rWork      ,
                                           Block         *block      )
{
    /* . Pair data. */
    Boolean        iIsJ, iIsK, jIsL, kIsL ;
    /* . Shell loops. */
    Boolean        iAndJ  , ijAndKL, kAndL ;
    Integer        iAMMax , iAMMaxT, iShell,         nCFuncI ,
                   jAMMax , jAMMaxT, jShell, jUpper, nCFuncJ ,
                   kAMMax , kAMMaxT, kShell, kUpper, nCFuncK ,
                   lAMMax , lAMMaxT, lShell, lUpper, nCFuncL ,
                   mAMMax , nAMMax , nRoots ,
                   sStrideI, sStrideIt, sStrideJ, sStrideJt , sStrideKL           , sStrideM ,
                   tStrideI, tStrideJ , tStrideK, tStrideKt , tStrideL , tStrideLt, tStrideM ;
    Real           xIJt, xKLt, yIJt, yKLt, zIJt, zKLt ;
    Real          *pG, *values = NULL, *work = NULL ;
    Real          *g, *gT ;
    const Real    *rC, *rD ;
    /* . Primitive loops. */
    Integer        iP, jP, kP, lP ;
    Real           aa, aaInv, aAndB, ab, aI, aJ, aK, aL, ar2I, ar2K, arg, argIJ, bb, bbInv,
                   axac, ayac, azac, axad, ayad, azad, bxbc, bybc, bzbc, bxbd, bybd, bzbd,
                   c1x , c2x , c3x , c4x , c1y , c2y , c3y , c4y , c1z , c2z , c3z , c4z,
                   expFac, rho, xAB , yAB , zAB ;
    Real           arI[3], arK[3], rA[3], rB[3],
                   Gx[MAXAMP21*MAXAMP21              ] ,
                   Gy[MAXAMP21*MAXAMP21              ] ,
                   Gz[MAXAMP21*MAXAMP21              ] ,
                   Sx[MAXAMP21*MAXAMP1*MAXAMP1       ] ,
                   Sy[MAXAMP21*MAXAMP1*MAXAMP1       ] ,
                   Sz[MAXAMP21*MAXAMP1*MAXAMP1       ] ,
                   Tx[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ,
                   Ty[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ,
                   Tz[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ;
    /* . Function loops. */
    Integer        f, i, ix, iy, iz, j, jix, jiy, jiz, k, kjix, kjiy, kjiz, l, m, nCFunc ;
    Integer       *Ix, *Iy, *Iz  ;
    Real           tI, tIJ, tIJK ;
    Real          *Cijkl ;
    RysQuadrature  roots ;
# ifdef _PrintIntegrals_
auto Integer ntotal = 0 ;
# endif
    /* . Initialization. */
    block->count = 0 ;
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    iIsK = ( iBasis == kBasis ) && ( rI == rK ) ;
    jIsL = ( jBasis == lBasis ) && ( rJ == rL ) ;
    kIsL = ( kBasis == lBasis ) && ( rK == rL ) ;
    /* . Set pointers. */
    Cijkl = &rWork[0] ; g  = &rWork[s4] ; gT = &rWork[2*s4] ;
    Ix    = &iWork[0] ; Iy = &iWork[s4] ; Iz = &iWork[2*s4] ;
    /* . Quadruple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iAMMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jAMMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            nAMMax  = iAMMax + jAMMax ;
            if ( iAMMax >= jAMMax )
            {
               iAMMaxT = iAMMax ;
               jAMMaxT = jAMMax ;
               xIJt    = rIJ[0] ;
               yIJt    = rIJ[1] ;
               zIJt    = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iAMMaxT = jAMMax ;
               jAMMaxT = iAMMax ;
               xIJt    = - rIJ[0] ;
               yIJt    = - rIJ[1] ;
               zIJt    = - rIJ[2] ;
               rC      = rJ ;
            }
            iAndJ = iIsJ && ( iShell == jShell ) ;
            if ( iIsK ) kUpper = iShell + 1 ;
            else        kUpper = kBasis->nShells ;
            for ( kShell = 0 ; kShell < kUpper ; kShell++ )
            {
                kAMMax  = kBasis->shells[kShell].lHigh ;
                nCFuncK = kBasis->shells[kShell].nCBF     ;
                if ( iIsK && ( iShell == kShell ) )
                {
                          if ( jIsL       ) { lUpper = jShell + 1      ; }
                     else if ( jLessThanL ) { lUpper = -1              ; }
                     else                   { lUpper = lBasis->nShells ; }
                }
                else
                {    if ( kIsL ) { lUpper = kShell + 1      ; }
                     else        { lUpper = lBasis->nShells ; }
                }
                for ( lShell = 0 ; lShell < lUpper ; lShell++ )
                {
                    lAMMax  = lBasis->shells[lShell].lHigh ;
                    nCFuncL = lBasis->shells[lShell].nCBF     ;
                    mAMMax  = kAMMax + lAMMax ;
                    if ( kAMMax >= lAMMax )
                    {
                       kAMMaxT = kAMMax ;
                       lAMMaxT = lAMMax ;
                       xKLt    = rKL[0] ;
                       yKLt    = rKL[1] ;
                       zKLt    = rKL[2] ;
                       rD      = rK ;
                    }
                    else
                    {
                       kAMMaxT = lAMMax ;
                       lAMMaxT = kAMMax ;
                       xKLt    = - rKL[0] ;
                       yKLt    = - rKL[1] ;
                       zKLt    = - rKL[2] ;
                       rD      = rL ;
                    }
                    kAndL   = kIsL && ( kShell == lShell ) ;
                    ijAndKL = iIsK && ( iShell == kShell ) && jIsL && ( jShell == lShell ) ;
                    /* . Roots. */
                    nRoots  = ( mAMMax + nAMMax ) / 2 + 1 ;
                    /* . Strides. */
                    sStrideKL = 1 ;
                    sStrideJ  = ( mAMMax + 1 ) * sStrideKL ;
                    sStrideI  = ( jAMMax + 1 ) * sStrideJ  ;
                    sStrideM  = ( iAMMax + 1 ) * sStrideI  ;
                    tStrideL  = 1 ;
                    tStrideK  = ( lAMMax + 1 ) * tStrideL  ;
                    tStrideJ  = ( kAMMax + 1 ) * tStrideK  ;
                    tStrideI  = ( jAMMax + 1 ) * tStrideJ  ;
                    tStrideM  = ( iAMMax + 1 ) * tStrideI  ;
                    if ( iAMMax >= jAMMax ) { sStrideIt = sStrideI ; sStrideJt = sStrideJ ; }
                    else                    { sStrideIt = sStrideJ ; sStrideJt = sStrideI ; }
                    if ( kAMMax >= lAMMax ) { tStrideKt = tStrideK ; tStrideLt = tStrideL ; }
                    else                    { tStrideKt = tStrideL ; tStrideLt = tStrideK ; }
                    /* . Index arrays. */
                    for ( i = 0, f = 0 ; i < nCFuncI ; i++ )
                    {
   	                ix = iBasis->shells[iShell].cbfPowX[i]*tStrideI ;
	                iy = iBasis->shells[iShell].cbfPowY[i]*tStrideI ;
	                iz = iBasis->shells[iShell].cbfPowZ[i]*tStrideI ;
                        for ( j = 0 ; j < nCFuncJ ; j++ )
                        {
	                    jix = jBasis->shells[jShell].cbfPowX[j]*tStrideJ + ix ;
	                    jiy = jBasis->shells[jShell].cbfPowY[j]*tStrideJ + iy ;
	                    jiz = jBasis->shells[jShell].cbfPowZ[j]*tStrideJ + iz ;
                            for ( k = 0 ; k < nCFuncK ; k++ )
                            {
                                kjix = kBasis->shells[kShell].cbfPowX[k]*tStrideK + jix ;
                                kjiy = kBasis->shells[kShell].cbfPowY[k]*tStrideK + jiy ;
                                kjiz = kBasis->shells[kShell].cbfPowZ[k]*tStrideK + jiz ;
                                for ( l = 0 ; l < nCFuncL ; l++, f++ )
                                {
                                    Ix[f] = lBasis->shells[lShell].cbfPowX[l] + kjix ;
                                    Iy[f] = lBasis->shells[lShell].cbfPowY[l] + kjiy ;
                                    Iz[f] = lBasis->shells[lShell].cbfPowZ[l] + kjiz ;
                                }
                            }
                        }
                    }
                    /* . Function initialization. */
                    nCFunc = nCFuncI * nCFuncJ * nCFuncK * nCFuncL ;
                    for ( f = 0 ; f < nCFunc ; f++ ) g[f] = 0.0e+00 ;
    /*------------------------------------------------------------------------------------------------------------------------------
    ! . Quadruple loop over primitives.
    !-----------------------------------------------------------------------------------------------------------------------------*/
    for ( iP = 0 ; iP < iBasis->shells[iShell].nPrimitives ; iP++ )
    {
	aI   = iBasis->shells[iShell].primitives[iP].exponent ;
	ar2I = aI * rIJ2 ;
	for ( i = 0 ; i < 3 ; i++ ) arI[i] = aI * rI[i] ;
        for ( jP = 0 ; jP < jBasis->shells[jShell].nPrimitives ; jP++ )
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
            for ( kP = 0 ; kP < kBasis->shells[kShell].nPrimitives ; kP++ )
            {
	        aK   = kBasis->shells[kShell].primitives[kP].exponent ;
	        ar2K = aK * rKL2 ;
	        for ( i = 0 ; i < 3 ; i++ ) arK[i] = aK * rK[i] ;
                for ( lP = 0 ; lP < lBasis->shells[lShell].nPrimitives ; lP++ )
                {
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
                    /* . Get coefficient array. */
                    for ( i = 0, f = 0 ; i < nCFuncI ; i++ )
                    {
                        tI = iBasis->shells[iShell].primitives[iP].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++ )
                        {
                            tIJ = tI * jBasis->shells[jShell].primitives[jP].cCBF[j] ;
                            for ( k = 0 ; k < nCFuncK ; k++ )
                            {
                                tIJK = tIJ * kBasis->shells[kShell].primitives[kP].cCBF[k] ;
                                for ( l = 0 ; l < nCFuncL ; l++, f++ ) Cijkl[f] = tIJK * lBasis->shells[lShell].primitives[lP].cCBF[l] ;
                            }
                        }
                    }
                    /* . Loop over Rys roots. */
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
                        GaussianBasisSubsidiary_f1Cg1  ( nAMMax    ,
                                                         mAMMax    ,
                                                         b00       ,
                                                         b10       ,
                                                         bp01      ,
                                                         f00       ,
                                                         xc00      ,
                                                         xcp00     ,
                                                         yc00      ,
                                                         ycp00     ,
                                                         zc00      ,
                                                         zcp00     ,
                                                         mAMMax+1  , /* . gStrideIJ. */
                                                         Gx        ,
                                                         Gy        ,
                                                         Gz        ) ;
                        /* . Reset S and T. */
                        for ( i = 0 ; i < sStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = 0.0e+00 ;
                        for ( i = 0 ; i < tStrideM ; i++ ) Tx[i] = Ty[i] = Tz[i] = 0.0e+00 ;
                        GaussianBasisSubsidiary_f1Xg2i ( iAMMaxT   ,
                                                         jAMMaxT   ,
                                                         mAMMax    ,
                                                         mAMMax+1  , /* . gStrideIJ. */
                                                         1         , /* . gStrideKL. */
                                                         Gx        ,
                                                         Gy        ,
                                                         Gz        ,
                                                         xIJt      ,
                                                         yIJt      ,
                                                         zIJt      ,
                                                         sStrideIt ,
                                                         sStrideJt ,
                                                         1         , /* . sStrideKL. */
                                                         Sx        ,
                                                         Sy        ,
                                                         Sz        ) ;
                        GaussianBasisSubsidiary_f1Xg2i ( kAMMaxT   ,
                                                         lAMMaxT   ,
                                                         (iAMMaxT+1)*(jAMMaxT+1)-1, /* . IJ loop range upper limit [0,u]. */
                                                         1         , /* . sStrideKL. */
                                                         sStrideJ  , /* . sStrideIJ. */
                                                         Sx        ,
                                                         Sy        ,
                                                         Sz        ,
                                                         xKLt      ,
                                                         yKLt      ,
                                                         zKLt      ,
                                                         tStrideKt ,
                                                         tStrideLt ,
                                                         tStrideJ  , /* . tStrideIJ. */
                                                         Tx        ,
                                                         Ty        ,
                                                         Tz        ) ;
                        /* . Assemble the integrals. */
                        for ( f = 0 ; f < nCFunc ; f++ ) g[f] += ( Cijkl[f] * Tx[Ix[f]] * Ty[Iy[f]] * Tz[Iz[f]] ) ;
                    } /* . nRoots. */
                } /* . lP. */
            } /* . kP. */
        } /* . jP. */
    } /* . iP. */
    /*----------------------------------------------------------------------------------------------------------------------------*/
    /* . Transform and save the integrals. */
    {
        auto Boolean     skip, skipIJ ;
        auto Integer     i, ii, ij, j, jj, k, kk, kl, l, ll, m, n ;
        auto Cardinal16 *indices16 = block->indices16 ;
        auto Real       *integrals = block->data      ;
        values = g  ;
        work   = gT ;
        GaussianBasisTransform4 ( nCFuncI                    ,
                                  nCFuncJ                    ,
                                  nCFuncK                    ,
                                  nCFuncL                    ,
                                  iBasis->shells[iShell].c2s ,
                                  jBasis->shells[jShell].c2s ,
                                  kBasis->shells[kShell].c2s ,
                                  lBasis->shells[lShell].c2s ,
                                  &values                    ,
                                  &work                      ) ;
        pG = values ;
        m  = block->count ;
        for ( i = 0, ij = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
        {
            ii = iBasis->shells[iShell].nStart + i ;
            for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; ij++, j++ )
            {
                jj     = jBasis->shells[jShell].nStart + j ;
                skipIJ = iAndJ && ( j > i ) ; /* . Do not skip here as need full loops for n counter. */
                for ( k = 0, kl = 0 ; k < kBasis->shells[kShell].nBasis ; k++ )
                {
                    kk = kBasis->shells[kShell].nStart + k ;
                    for ( l = 0 ; l < lBasis->shells[lShell].nBasis ; kl++, l++, n++ )
                    {
                        ll   = lBasis->shells[lShell].nStart + l ;
                        skip = skipIJ || ( ijAndKL && ( ij < kl ) ) || ( kAndL && ( l > k ) ) ;
# ifdef _PrintIntegrals_
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
                            integrals[m]    = pG[n] ;
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
# undef MAXAMP21

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the two-electron integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 6 * s4 and real 11 * s4 where s4 = ( maximum shell size )^4. */
# define _DensityTolerance 1.0e-12
# define  MAXAMP22         ( 2 * MAXAMP1     )
# define  MAXAMP23         ( 2 * MAXAMP1 + 1 )
# define _IntegralSizeD    ( MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1 )
# define _IntegralSizeG    ( MAXAMP23*MAXAMP22        )
# define _IntegralSizeS    ( MAXAMP22*MAXAMP2*MAXAMP2 )
# define _IntegralSizeT    ( MAXAMP2*MAXAMP2*MAXAMP2*MAXAMP1 )
void GaussianBasisIntegrals_f2Cf2r1 ( const GaussianBasis *iBasis     ,
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
                                      const Integer        s4         ,
                                            Integer       *iWork      ,
                                            Real          *rWork      ,
                                            Block         *block      )
{
    /* . Pair data. */
    Boolean        iIsJ, iIsK, jIsL, kIsL ;
    /* . Shell loops. */
    Boolean        iAndJ  , ijAndKL, kAndL ;
    Integer        iAMMax , iAMMaxT, iShell,         nCFuncI ,
                   jAMMax , jAMMaxT, jShell, jUpper, nCFuncJ ,
                   kAMMax , kAMMaxT, kShell, kUpper, nCFuncK ,
                   lAMMax , lAMMaxT, lShell, lUpper, nCFuncL ,
                   mAMMax , nAMMax , nRoots ,
                   dStrideI, dStrideJ , dStrideK, dStrideL  ,
                   sStrideI, sStrideIt, sStrideJ, sStrideJt , sStrideKL           , sStrideM ,
                   tStrideI, tStrideJ , tStrideK, tStrideKt , tStrideL , tStrideLt, tStrideM ;
    Real           xIJt, xKLt, yIJt, yKLt, zIJt, zKLt ;
    Real          *pGIx, *pGIy, *pGIz, *pGJx, *pGJy, *pGJz, *pGKx, *pGKy, *pGKz,
                  *values = NULL, *work = NULL ;
    Real          *gIx, *gIy, *gIz ,
                  *gJx, *gJy, *gJz ,
                  *gKx, *gKy, *gKz ,
                  *gT ;
    const Real    *rC, *rD ;
    /* . Primitive loops. */
    Integer        iP, jP, kP, lP ;
    Real           aa, aaInv, aAndB, ab, aI, aJ, aK, aL, ar2I, ar2K, arg, argIJ, bb, bbInv,
                   axac, ayac, azac, axad, ayad, azad, bxbc, bybc, bzbc, bxbd, bybd, bzbd,
                   c1x , c2x , c3x , c4x , c1y , c2y , c3y , c4y , c1z , c2z , c3z , c4z,
                   expFac, rho, xAB , yAB , zAB ;
    Real           arI[3], arK[3], rA[3], rB[3],
                   Gx [_IntegralSizeG] , Gy [_IntegralSizeG] , Gz [_IntegralSizeG] ,
                   Sx [_IntegralSizeS] , Sy [_IntegralSizeS] , Sz [_IntegralSizeS] ,
                   Tx [_IntegralSizeT] , Ty [_IntegralSizeT] , Tz [_IntegralSizeT] ,
                   xDI[_IntegralSizeD] , yDI[_IntegralSizeD] , zDI[_IntegralSizeD] ,
                   xDJ[_IntegralSizeD] , yDJ[_IntegralSizeD] , zDJ[_IntegralSizeD] ,
                   xDK[_IntegralSizeD] , yDK[_IntegralSizeD] , zDK[_IntegralSizeD] ;
    /* . Function loops. */
    Integer        f, i, ix, ixd, iy, iyd, iz, izd, j, jix, jixd, jiy, jiyd, jiz, jizd,
                      k, kjix, kjixd, kjiy, kjiyd, kjiz, kjizd, l, m, nCFunc ;
    Integer       *Ix, *Ixd ,
                  *Iy, *Iyd ,
                  *Iz, *Izd ;
    Real           tI, tIJ, tIJK ;
    Real          *Cijkl ;
    RysQuadrature  roots ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    iIsK = ( iBasis == kBasis ) && ( rI == rK ) ;
    jIsL = ( jBasis == lBasis ) && ( rJ == rL ) ;
    kIsL = ( kBasis == lBasis ) && ( rK == rL ) ;
    /* . Set pointers. */
    Cijkl = &rWork[    0] ;
    gIx   = &rWork[   s4] ; gIy = &rWork[2*s4] ; gIz = &rWork[3*s4] ;
    gJx   = &rWork[ 4*s4] ; gJy = &rWork[5*s4] ; gJz = &rWork[6*s4] ;
    gKx   = &rWork[ 7*s4] ; gKy = &rWork[8*s4] ; gKz = &rWork[9*s4] ;
    gT    = &rWork[10*s4] ;
    Ix    = &iWork[    0] ; Iy  = &iWork[  s4] ; Iz  = &iWork[2*s4] ;
    Ixd   = &iWork[ 3*s4] ; Iyd = &iWork[4*s4] ; Izd = &iWork[5*s4] ;
    /* . Quadruple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iAMMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jAMMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            nAMMax  = iAMMax + jAMMax + 2 ; /* . I and J increased by 1. */
            if ( iAMMax >= jAMMax )
            {
               iAMMaxT = iAMMax + 1 ; /* . AMMax normal, AMMaxT increased by 1. */
               jAMMaxT = jAMMax + 1 ;
               xIJt    = rIJ[0] ;
               yIJt    = rIJ[1] ;
               zIJt    = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iAMMaxT = jAMMax + 1 ;
               jAMMaxT = iAMMax + 1 ;
               xIJt    = - rIJ[0] ;
               yIJt    = - rIJ[1] ;
               zIJt    = - rIJ[2] ;
               rC      = rJ ;
            }
            iAndJ = iIsJ && ( iShell == jShell ) ;
            if ( iIsK ) kUpper = iShell + 1 ;
            else        kUpper = kBasis->nShells ;
            for ( kShell = 0 ; kShell < kUpper ; kShell++ )
            {
                kAMMax  = kBasis->shells[kShell].lHigh ;
                nCFuncK = kBasis->shells[kShell].nCBF     ;
                if ( iIsK && ( iShell == kShell ) )
                {
                          if ( jIsL       ) { lUpper = jShell + 1      ; }
                     else if ( jLessThanL ) { lUpper = -1              ; }
                     else                   { lUpper = lBasis->nShells ; }
                 }
                else
                {    if ( kIsL ) { lUpper = kShell + 1      ; }
                     else        { lUpper = lBasis->nShells ; }
                }
                for ( lShell = 0 ; lShell < lUpper ; lShell++ )
                {
                    lAMMax  = lBasis->shells[lShell].lHigh ;
                    nCFuncL = lBasis->shells[lShell].nCBF     ;
                    mAMMax  = kAMMax + lAMMax + 1 ; /* . K increased by 1. */
                    if ( ( kAMMax + 1 ) >= lAMMax )
                    {
                       kAMMaxT = kAMMax + 1 ; /* . AMMax normal, AMMaxT increased by 1. */
                       lAMMaxT = lAMMax ;
                       xKLt    = rKL[0] ;
                       yKLt    = rKL[1] ;
                       zKLt    = rKL[2] ;
                       rD      = rK ;
                    }
                    else
                    {
                       kAMMaxT = lAMMax ;
                       lAMMaxT = kAMMax + 1 ;
                       xKLt    = - rKL[0] ;
                       yKLt    = - rKL[1] ;
                       zKLt    = - rKL[2] ;
                       rD      = rL ;
                    }
                    kAndL   = kIsL && ( kShell == lShell ) ;
                    ijAndKL = iIsK && ( iShell == kShell ) && jIsL && ( jShell == lShell ) ;
                    /* . Integral initialization. */
                    nRoots  = ( mAMMax + nAMMax ) / 2 + 1 ;
                    /* . Strides. */
                    sStrideKL = 1 ;
                    sStrideJ  = ( mAMMax + 1 ) * sStrideKL ;
                    sStrideI  = ( jAMMax + 2 ) * sStrideJ  ;
                    sStrideM  = ( iAMMax + 2 ) * sStrideI  ;
                    dStrideL  = 1                         ; tStrideL = 1                         ;
                    dStrideK  = ( lAMMax + 1 ) * dStrideL ; tStrideK = ( lAMMax + 1 ) * tStrideL ;
                    dStrideJ  = ( kAMMax + 1 ) * dStrideK ; tStrideJ = ( kAMMax + 2 ) * tStrideK ;
                    dStrideI  = ( jAMMax + 1 ) * dStrideJ ; tStrideI = ( jAMMax + 2 ) * tStrideJ ;
                                                            tStrideM = ( iAMMax + 2 ) * tStrideI ;
                    if ( iAMMax     >= jAMMax ) { sStrideIt = sStrideI ; sStrideJt = sStrideJ ; }
                    else                        { sStrideIt = sStrideJ ; sStrideJt = sStrideI ; }
                    if ( kAMMax + 1 >= lAMMax ) { tStrideKt = tStrideK ; tStrideLt = tStrideL ; }
                    else                        { tStrideKt = tStrideL ; tStrideLt = tStrideK ; }
                    /* . Get index arrays. */
                    for ( i = 0, f = 0 ; i < nCFuncI ; i++ )
                    {
                        ix  = iBasis->shells[iShell].cbfPowX[i]*tStrideI ;
                        iy  = iBasis->shells[iShell].cbfPowY[i]*tStrideI ;
                        iz  = iBasis->shells[iShell].cbfPowZ[i]*tStrideI ;
                        ixd = iBasis->shells[iShell].cbfPowX[i]*dStrideI ;
                        iyd = iBasis->shells[iShell].cbfPowY[i]*dStrideI ;
                        izd = iBasis->shells[iShell].cbfPowZ[i]*dStrideI ;
                        for ( j = 0 ; j < nCFuncJ ; j++ )
                        {
                            jix  = jBasis->shells[jShell].cbfPowX[j]*tStrideJ + ix  ;
                            jiy  = jBasis->shells[jShell].cbfPowY[j]*tStrideJ + iy  ;
                            jiz  = jBasis->shells[jShell].cbfPowZ[j]*tStrideJ + iz  ;
                            jixd = jBasis->shells[jShell].cbfPowX[j]*dStrideJ + ixd ;
                            jiyd = jBasis->shells[jShell].cbfPowY[j]*dStrideJ + iyd ;
                            jizd = jBasis->shells[jShell].cbfPowZ[j]*dStrideJ + izd ;
                            for ( k = 0 ; k < nCFuncK ; k++ )
                            {
                                kjix  = kBasis->shells[kShell].cbfPowX[k]*tStrideK + jix  ;
                                kjiy  = kBasis->shells[kShell].cbfPowY[k]*tStrideK + jiy  ;
                                kjiz  = kBasis->shells[kShell].cbfPowZ[k]*tStrideK + jiz  ;
                                kjixd = kBasis->shells[kShell].cbfPowX[k]*dStrideK + jixd ;
                                kjiyd = kBasis->shells[kShell].cbfPowY[k]*dStrideK + jiyd ;
                                kjizd = kBasis->shells[kShell].cbfPowZ[k]*dStrideK + jizd ;
                                for ( l = 0 ; l < nCFuncL ; l++, f++ )
                                {
                                    Ix [f] = lBasis->shells[lShell].cbfPowX[l] + kjix  ;
                                    Iy [f] = lBasis->shells[lShell].cbfPowY[l] + kjiy  ;
                                    Iz [f] = lBasis->shells[lShell].cbfPowZ[l] + kjiz  ;
                                    Ixd[f] = lBasis->shells[lShell].cbfPowX[l] + kjixd ;
                                    Iyd[f] = lBasis->shells[lShell].cbfPowY[l] + kjiyd ;
                                    Izd[f] = lBasis->shells[lShell].cbfPowZ[l] + kjizd ;
                                }
                            }
                        }
                    }
                    /* . Function initialization. */
                    nCFunc = nCFuncI * nCFuncJ * nCFuncK * nCFuncL ;
                    for ( f = 0 ; f < nCFunc ; f++ )
                    {
                        gIx[f] = gIy[f] = gIz[f] = gJx[f] = gJy[f] = gJz[f] = gKx[f] = gKy[f] = gKz[f] = 0.0e+00 ;
                    }
    /*------------------------------------------------------------------------------------------------------------------------------
    ! . Quadruple loop over primitives.
    !-----------------------------------------------------------------------------------------------------------------------------*/
    for ( iP = 0 ; iP < iBasis->shells[iShell].nPrimitives ; iP++ )
    {
        aI   = iBasis->shells[iShell].primitives[iP].exponent ;
        ar2I = aI * rIJ2 ;
        for ( i = 0 ; i < 3 ; i++ ) arI[i] = aI * rI[i] ;
        for ( jP = 0 ; jP < jBasis->shells[jShell].nPrimitives ; jP++ )
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
            for ( kP = 0 ; kP < kBasis->shells[kShell].nPrimitives ; kP++ )
            {
                aK   = kBasis->shells[kShell].primitives[kP].exponent ;
                ar2K = aK * rKL2 ;
                for ( i = 0 ; i < 3 ; i++ ) arK[i] = aK * rK[i] ;
                for ( lP = 0 ; lP < lBasis->shells[lShell].nPrimitives ; lP++ )
                {
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
                    /* . Get coefficient array. */
                    for ( i = 0, f = 0 ; i < nCFuncI ; i++ )
                    {
                        tI = iBasis->shells[iShell].primitives[iP].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++ )
                        {
                            tIJ = tI * jBasis->shells[jShell].primitives[jP].cCBF[j] ;
                            for ( k = 0 ; k < nCFuncK ; k++ )
                            {
                                tIJK = tIJ * kBasis->shells[kShell].primitives[kP].cCBF[k] ;
                                for ( l = 0 ; l < nCFuncL ; l++, f++ ) Cijkl[f] = tIJK * lBasis->shells[lShell].primitives[lP].cCBF[l] ;
                            }
                        }
                    }
                    /* . Loop over Rys roots. */
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
                        GaussianBasisSubsidiary_f1Cg1  ( nAMMax    ,
                                                         mAMMax    ,
                                                         b00       ,
                                                         b10       ,
                                                         bp01      ,
                                                         f00       ,
                                                         xc00      ,
                                                         xcp00     ,
                                                         yc00      ,
                                                         ycp00     ,
                                                         zc00      ,
                                                         zcp00     ,
                                                         mAMMax+1  , /* . gStrideIJ. */
                                                         Gx        ,
                                                         Gy        ,
                                                         Gz        ) ;
                        /* . Reset S and T. */
                        for ( i = 0 ; i < sStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = 0.0e+00 ;
                        for ( i = 0 ; i < tStrideM ; i++ ) Tx[i] = Ty[i] = Tz[i] = 0.0e+00 ;
                        GaussianBasisSubsidiary_f1Xg2i ( iAMMaxT   ,
                                                         jAMMaxT   ,
                                                         mAMMax    ,
                                                         mAMMax+1  , /* . gStrideIJ. */
                                                         1         , /* . gStrideKL. */
                                                         Gx        ,
                                                         Gy        ,
                                                         Gz        ,
                                                         xIJt      ,
                                                         yIJt      ,
                                                         zIJt      ,
                                                         sStrideIt ,
                                                         sStrideJt ,
                                                         1         , /* . sStrideKL. */
                                                         Sx        ,
                                                         Sy        ,
                                                         Sz        ) ;
                        GaussianBasisSubsidiary_f1Xg2i ( kAMMaxT   ,
                                                         lAMMaxT   ,
                                                         (iAMMaxT+1)*(jAMMaxT+1)-1, /* . IJ loop range upper limit [0,u]. */
                                                         1         , /* . sStrideKL. */
                                                         sStrideJ  , /* . sStrideIJ. */
                                                         Sx        ,
                                                         Sy        ,
                                                         Sz        ,
                                                         xKLt      ,
                                                         yKLt      ,
                                                         zKLt      ,
                                                         tStrideKt ,
                                                         tStrideLt ,
                                                         tStrideJ  , /* . tStrideIJ. */
                                                         Tx        ,
                                                         Ty        ,
                                                         Tz        ) ;
                        GaussianBasisSubsidiary_f2Xg2r ( iAMMax    ,
                                                         jAMMax    ,
                                                         kAMMax    ,
                                                         lAMMax    ,
                                                         tStrideI  ,
                                                         tStrideJ  ,
                                                         tStrideK  ,
                                                         tStrideL  ,
                                                         dStrideI  ,
                                                         dStrideJ  ,
                                                         dStrideK  ,
                                                         dStrideL  ,
                                                         aI        ,
                                                         aJ        ,
                                                         aK        ,
                                                         Tx        ,
                                                         Ty        ,
                                                         Tz        ,
                                                         xDI       ,
                                                         yDI       ,
                                                         zDI       ,
                                                         xDJ       ,
                                                         yDJ       ,
                                                         zDJ       ,
                                                         xDK       ,
                                                         yDK       ,
                                                         zDK       ) ;
                        /* . Assemble the integrals. */
                        for ( f = 0 ; f < nCFunc ; f++ )
                        {
                            gIx[f] += ( Cijkl[f] * xDI[Ixd[f]] * Ty [Iy [f]] * Tz [Iz [f]] ) ;
                            gIy[f] += ( Cijkl[f] * Tx [Ix [f]] * yDI[Iyd[f]] * Tz [Iz [f]] ) ;
                            gIz[f] += ( Cijkl[f] * Tx [Ix [f]] * Ty [Iy [f]] * zDI[Izd[f]] ) ;
                            gJx[f] += ( Cijkl[f] * xDJ[Ixd[f]] * Ty [Iy [f]] * Tz [Iz [f]] ) ;
                            gJy[f] += ( Cijkl[f] * Tx [Ix [f]] * yDJ[Iyd[f]] * Tz [Iz [f]] ) ;
                            gJz[f] += ( Cijkl[f] * Tx [Ix [f]] * Ty [Iy [f]] * zDJ[Izd[f]] ) ;
                            gKx[f] += ( Cijkl[f] * xDK[Ixd[f]] * Ty [Iy [f]] * Tz [Iz [f]] ) ;
                            gKy[f] += ( Cijkl[f] * Tx [Ix [f]] * yDK[Iyd[f]] * Tz [Iz [f]] ) ;
                            gKz[f] += ( Cijkl[f] * Tx [Ix [f]] * Ty [Iy [f]] * zDK[Izd[f]] ) ;
                        }
                    } /* . nRoots. */
                } /* . lP. */
            } /* . kP. */
        } /* . jP. */
    } /* . iP. */
    /*----------------------------------------------------------------------------------------------------------------------------*/
    /* . Transform and save the integrals. */
    {
        auto Boolean      skip, skipIJ ;
        auto Integer      i, ii, ij, j, jj, k, kk, kl, l, ll, m, n ;
        auto Cardinal16  *indices16 = block->indices16 ;
        auto Real        *integrals = block->data      ;
        auto RealArray2D *iC2S = iBasis->shells[iShell].c2s ,
                         *jC2S = jBasis->shells[jShell].c2s ,
                         *kC2S = kBasis->shells[kShell].c2s ,
                         *lC2S = lBasis->shells[lShell].c2s ;
        work   = gT  ;
        values = gIx ; GaussianBasisTransform4 ( nCFuncI, nCFuncJ, nCFuncK, nCFuncL, iC2S, jC2S, kC2S, lC2S, &values, &work ) ; pGIx = values ;
        values = gIy ; GaussianBasisTransform4 ( nCFuncI, nCFuncJ, nCFuncK, nCFuncL, iC2S, jC2S, kC2S, lC2S, &values, &work ) ; pGIy = values ;
        values = gIz ; GaussianBasisTransform4 ( nCFuncI, nCFuncJ, nCFuncK, nCFuncL, iC2S, jC2S, kC2S, lC2S, &values, &work ) ; pGIz = values ;
        values = gJx ; GaussianBasisTransform4 ( nCFuncI, nCFuncJ, nCFuncK, nCFuncL, iC2S, jC2S, kC2S, lC2S, &values, &work ) ; pGJx = values ;
        values = gJy ; GaussianBasisTransform4 ( nCFuncI, nCFuncJ, nCFuncK, nCFuncL, iC2S, jC2S, kC2S, lC2S, &values, &work ) ; pGJy = values ;
        values = gJz ; GaussianBasisTransform4 ( nCFuncI, nCFuncJ, nCFuncK, nCFuncL, iC2S, jC2S, kC2S, lC2S, &values, &work ) ; pGJz = values ;
        values = gKx ; GaussianBasisTransform4 ( nCFuncI, nCFuncJ, nCFuncK, nCFuncL, iC2S, jC2S, kC2S, lC2S, &values, &work ) ; pGKx = values ;
        values = gKy ; GaussianBasisTransform4 ( nCFuncI, nCFuncJ, nCFuncK, nCFuncL, iC2S, jC2S, kC2S, lC2S, &values, &work ) ; pGKy = values ;
        values = gKz ; GaussianBasisTransform4 ( nCFuncI, nCFuncJ, nCFuncK, nCFuncL, iC2S, jC2S, kC2S, lC2S, &values, &work ) ; pGKz = values ;
        m         = block->count ;
        for ( i = 0, ij = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
        {
            ii = iBasis->shells[iShell].nStart + i ;
            for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; ij++, j++ )
            {
                jj     = jBasis->shells[jShell].nStart + j ;
                skipIJ = iAndJ && ( j > i ) ; /* . Do not skip here as need full loops for n counter. */
                for ( k = 0, kl = 0 ; k < kBasis->shells[kShell].nBasis ; k++ )
                {
                    kk = kBasis->shells[kShell].nStart + k ;
                    for ( l = 0 ; l < lBasis->shells[lShell].nBasis ; kl++, l++, n++ )
                    {
                        ll   = lBasis->shells[lShell].nStart + l ;
                        skip = skipIJ || ( ijAndKL && ( ij < kl ) ) || ( kAndL && ( l > k ) ) ;
                        if ( ! skip )
                        {
                            auto Integer  m4 = 4 * m, m9 = 9 * m ;
                            indices16[m4  ] = ii ;
                            indices16[m4+1] = jj ;
                            indices16[m4+2] = kk ;
                            indices16[m4+3] = ll ;
                            integrals[m9  ] = pGIx[n] ;
                            integrals[m9+1] = pGIy[n] ;
                            integrals[m9+2] = pGIz[n] ;
                            integrals[m9+3] = pGJx[n] ;
                            integrals[m9+4] = pGJy[n] ;
                            integrals[m9+5] = pGJz[n] ;
                            integrals[m9+6] = pGKx[n] ;
                            integrals[m9+7] = pGKy[n] ;
                            integrals[m9+8] = pGKz[n] ;
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
# undef  MAXAMP22        
# undef  MAXAMP23        
# undef _IntegralSizeC   
# undef _IntegralSizeD   
# undef _IntegralSizeG   
# undef _IntegralSizeS   
# undef _IntegralSizeT

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the one-electron four center overlap integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s4 and real 2 * s4 where s4 = ( maximum shell size )^4. */
void GaussianBasisIntegrals_f2Of2i ( const GaussianBasis *iBasis     ,
                                     const Real          *rI         ,
                                     const GaussianBasis *jBasis     ,
                                     const Real          *rJ         ,
                                     const GaussianBasis *kBasis     ,
                                     const Real          *rK         ,
                                     const GaussianBasis *lBasis     ,
                                     const Real          *rL         ,
                                     const Boolean        jLessThanL ,
                                     const Integer        s4         ,
                                           Integer       *iWork      ,
                                           Real          *rWork      ,
                                           Block         *block      )
{
    /* . Pair data. */
    Boolean  iIsJ, iIsK, jIsL, kIsL ;
    /* . Shell loops. */
    Boolean  iAndJ  , ijAndKL, kAndL ;
    Integer  iAMMax , iShell,         nCFuncI ,
             jAMMax , jShell, jUpper, nCFuncJ ,
             kAMMax , kShell, kUpper, nCFuncK ,
             lAMMax , lShell, lUpper, nCFuncL ,
             strideI, strideJ, strideK ;
    Real     rIJ2, rIK2, rIL2, rJK2, rJL2, rKL2 ;
    Real    *pG, *values = NULL, *work = NULL ;
    Real    *g, *gT ;
    /* . Primitive loops. */
    Integer  iP, jP, kP, lP ;
    Real     aI, aIJ, aIJK, aIJKL, aJ, aK, aL, eIJ, eIJK, eIJKL, expFac ;
    Real     cI[3], cIJ[3], cIJK[3], cIJKL[3],
             Sx[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ,
             Sy[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ,
             Sz[MAXAMP1*MAXAMP1*MAXAMP1*MAXAMP1] ;
    /* . Function loops. */
    Integer  c, f, i, ix, iy, iz, j, jix, jiy, jiz, k, kjix, kjiy, kjiz, l, nCFunc ;
    Integer *Ix, *Iy, *Iz  ;
    Real     tI, tIJ, tIJK ;
# ifdef _PrintIntegrals_
auto Integer ntotal = 0 ;
# endif
    /* . Initialization. */
    block->count = 0 ;
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    iIsK = ( iBasis == kBasis ) && ( rI == rK ) ;
    jIsL = ( jBasis == lBasis ) && ( rJ == rL ) ;
    kIsL = ( kBasis == lBasis ) && ( rK == rL ) ;
    for ( c = 0, rIJ2 = rIK2 = rIL2 = rJK2 = rJL2 = rKL2 = 0.0e+00 ; c < 3 ; c++ )
    {
        rIJ2 += pow ( rI[c] - rJ[c], 2 ) ;
        rIK2 += pow ( rI[c] - rK[c], 2 ) ;
        rIL2 += pow ( rI[c] - rL[c], 2 ) ;
        rJK2 += pow ( rJ[c] - rK[c], 2 ) ;
        rJL2 += pow ( rJ[c] - rL[c], 2 ) ;
        rKL2 += pow ( rK[c] - rL[c], 2 ) ;
    }
    /* . Set pointers. */
    g  = &rWork[0] ; gT = &rWork[s4] ;
    Ix = &iWork[0] ; Iy = &iWork[s4] ; Iz = &iWork[2*s4] ;
    /* . Quadruple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iAMMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jAMMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            iAndJ   = iIsJ && ( iShell == jShell ) ;
            if ( iIsK ) kUpper = iShell + 1 ;
            else        kUpper = kBasis->nShells ;
            for ( kShell = 0 ; kShell < kUpper ; kShell++ )
            {
                kAMMax  = kBasis->shells[kShell].lHigh ;
                nCFuncK = kBasis->shells[kShell].nCBF     ;
                if ( iIsK && ( iShell == kShell ) )
                {
                          if ( jIsL       ) { lUpper = jShell + 1      ; }
                     else if ( jLessThanL ) { lUpper = -1              ; }
                     else                   { lUpper = lBasis->nShells ; }
                }
                else
                {    if ( kIsL ) { lUpper = kShell + 1      ; }
                     else        { lUpper = lBasis->nShells ; }
                }
                for ( lShell = 0 ; lShell < lUpper ; lShell++ )
                {
                    lAMMax  = lBasis->shells[lShell].lHigh ;
                    nCFuncL = lBasis->shells[lShell].nCBF     ;
                    kAndL   = kIsL && ( kShell == lShell ) ;
                    ijAndKL = iIsK && ( iShell == kShell ) && jIsL && ( jShell == lShell ) ;
                    /* . Strides. */
                    /*strideL = 1 ;*/
                    strideK = ( lAMMax + 1 ) ; /* *strideL ;*/
                    strideJ = ( kAMMax + 1 ) * strideK ;
                    strideI = ( jAMMax + 1 ) * strideJ ;
                    /* . Index arrays. */
                    for ( i = 0, f = 0 ; i < nCFuncI ; i++ )
                    {
   	                ix = iBasis->shells[iShell].cbfPowX[i]*strideI ;
	                iy = iBasis->shells[iShell].cbfPowY[i]*strideI ;
	                iz = iBasis->shells[iShell].cbfPowZ[i]*strideI ;
                        for ( j = 0 ; j < nCFuncJ ; j++ )
                        {
	                    jix = jBasis->shells[jShell].cbfPowX[j]*strideJ + ix ;
	                    jiy = jBasis->shells[jShell].cbfPowY[j]*strideJ + iy ;
	                    jiz = jBasis->shells[jShell].cbfPowZ[j]*strideJ + iz ;
                            for ( k = 0 ; k < nCFuncK ; k++ )
                            {
                                kjix = kBasis->shells[kShell].cbfPowX[k]*strideK + jix ;
                                kjiy = kBasis->shells[kShell].cbfPowY[k]*strideK + jiy ;
                                kjiz = kBasis->shells[kShell].cbfPowZ[k]*strideK + jiz ;
                                for ( l = 0 ; l < nCFuncL ; l++, f++ )
                                {
                                    Ix[f] = lBasis->shells[lShell].cbfPowX[l] + kjix ;
                                    Iy[f] = lBasis->shells[lShell].cbfPowY[l] + kjiy ;
                                    Iz[f] = lBasis->shells[lShell].cbfPowZ[l] + kjiz ;
                                }
                            }
                        }
                    }
                    /* . Function initialization. */
                    nCFunc = nCFuncI * nCFuncJ * nCFuncK * nCFuncL ;
                    for ( f = 0 ; f < nCFunc ; f++ ) g[f] = 0.0e+00 ;
    /*------------------------------------------------------------------------------------------------------------------------------
    ! . Quadruple loop over primitives.
    !-----------------------------------------------------------------------------------------------------------------------------*/
    for ( iP = 0 ; iP < iBasis->shells[iShell].nPrimitives ; iP++ )
    {
	aI = iBasis->shells[iShell].primitives[iP].exponent ;
	for ( c = 0 ; c < 3 ; c++ ) cI[c] = aI * rI[c] ;
        for ( jP = 0 ; jP < jBasis->shells[jShell].nPrimitives ; jP++ )
        {
	    aJ  = jBasis->shells[jShell].primitives[jP].exponent ;
            aIJ = aI + aJ ;
	    eIJ = aI * aJ * rIJ2 ;
	    if ( ( eIJ / aIJ ) > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
	    for ( c = 0 ; c < 3 ; c++ ) cIJ[c] = ( cI[c] + aJ * rJ[c] ) ;
            for ( kP = 0 ; kP < kBasis->shells[kShell].nPrimitives ; kP++ )
            {
	        aK   = kBasis->shells[kShell].primitives[kP].exponent ;
                aIJK = aIJ + aK ;
                eIJK = ( eIJ + aI * aK * rIK2 + aJ * aK * rJK2 ) ;
	        if ( ( eIJK / aIJK ) > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
	        for ( c = 0 ; c < 3 ; c++ ) cIJK[c] = ( cIJ[c] + aK * rK[c] ) ;
                for ( lP = 0 ; lP < lBasis->shells[lShell].nPrimitives ; lP++ )
                {
	            aL    = lBasis->shells[lShell].primitives[lP].exponent ;
                    aIJKL = aIJK + aL ;
                    eIJKL = ( eIJK + aI * aL * rIL2 + aJ * aL * rJL2 + aK * aL * rKL2 ) / aIJKL ;
	            if ( eIJKL > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    expFac = exp ( - eIJKL ) ;
	            for ( c = 0 ; c < 3 ; c++ ) cIJKL[c] = ( cIJK[c] + aL * rL[c] ) / aIJKL ;
                    /* . Calculate the overlap integrals. */
                    GaussianBasisSubsidiary_f2Og2 ( Sx, Sy, Sz, aIJKL, cIJKL, rI, rJ, rK, rL, iAMMax, jAMMax, kAMMax, lAMMax ) ;
                    /* . Assemble the integrals. */
                    for ( i = 0, f = 0 ; i < nCFuncI ; i++ )
                    {
                        tI = expFac * iBasis->shells[iShell].primitives[iP].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++ )
                        {
                            tIJ = tI * jBasis->shells[jShell].primitives[jP].cCBF[j] ;
                            for ( k = 0 ; k < nCFuncK ; k++ )
                            {
                                tIJK = tIJ * kBasis->shells[kShell].primitives[kP].cCBF[k] ;
                                for ( l = 0 ; l < nCFuncL ; l++, f++ )
                                {
                                    g[f] += ( tIJK * lBasis->shells[lShell].primitives[lP].cCBF[l] * Sx[Ix[f]] * Sy[Iy[f]] * Sz[Iz[f]] ) ;
                                }
                            }
                        }
                    }
                } /* . lP. */
            } /* . kP. */
        } /* . jP. */
    } /* . iP. */
    /*----------------------------------------------------------------------------------------------------------------------------*/
    /* . Transform and save the integrals. */
    {
        auto Boolean     skip, skipIJ ;
        auto Integer     i, ii, ij, j, jj, k, kk, kl, l, ll, m, n ;
        auto Cardinal16 *indices16 = block->indices16 ;
        auto Real       *integrals = block->data      ;
        values = g  ;
        work   = gT ;
        GaussianBasisTransform4 ( nCFuncI                    ,
                                  nCFuncJ                    ,
                                  nCFuncK                    ,
                                  nCFuncL                    ,
                                  iBasis->shells[iShell].c2s ,
                                  jBasis->shells[jShell].c2s ,
                                  kBasis->shells[kShell].c2s ,
                                  lBasis->shells[lShell].c2s ,
                                  &values                    ,
                                  &work                      ) ;
        pG = values ;
        m  = block->count ;
        for ( i = 0, ij = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
        {
            ii = iBasis->shells[iShell].nStart + i ;
            for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; ij++, j++ )
            {
                jj     = jBasis->shells[jShell].nStart + j ;
                skipIJ = iAndJ && ( j > i ) ; /* . Do not skip here as need full loops for n counter. */
                for ( k = 0, kl = 0 ; k < kBasis->shells[kShell].nBasis ; k++ )
                {
                    kk = kBasis->shells[kShell].nStart + k ;
                    for ( l = 0 ; l < lBasis->shells[lShell].nBasis ; kl++, l++, n++ )
                    {
                        ll   = lBasis->shells[lShell].nStart + l ;
                        skip = skipIJ || ( ijAndKL && ( ij < kl ) ) || ( kAndL && ( l > k ) ) ;
# ifdef _PrintIntegrals_
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
                            integrals[m]    = pG[n] ;
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

