/*==================================================================================================================================
! . Integrals - 2 basis, 2 electrons, 1 nucleus/point.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisIntegrals_f2Cp1.h"
# include "GaussianBasisSubsidiary.h"
# include "GaussianBasisTransform.h"
# include "GaussianNucleus.h"
# include "RysQuadrature.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Selected( selection, i ) ( (selection) == NULL ? True : Block_Item ( selection->flags, i ) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-nuclear/point derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 6 * s2 and real 8 * s2 where s2 = ( maximum shell size )^2. */
# define MAXAMP23 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP3 )
void GaussianBasisIntegrals_f2Cm1R1 ( const GaussianBasis *iBasis     ,
                                      const Real          *rI         ,
                                      const GaussianBasis *jBasis     ,
                                      const Real          *rJ         ,
                                      const RealArray1D   *charges    ,
                                      const RealArray1D   *widthsE    ,
                                      const RealArray1D   *widthsN    ,
                                      const Coordinates3  *rNP        ,
                                      const Selection     *selectionN ,
                                      const RealArray2D   *dOneIJ     ,
                                      const Integer        s2         ,
                                            Integer       *iWork      ,
                                            Real          *rWork      ,
                                            Real          *dRi        ,
                                            Real          *dRj        ,
                                            Coordinates3  *gN         )
{
    Boolean       iIsJ, isDiagonal ;
    Integer       dStrideI, dStrideJ,
                  i, iAMMax, iAMMaxT, ii, ip, iShell, ix, ixd, iy, iyd, iz, izd,
                  j, jAMMax, jAMMaxT, jj, jp, jShell, jUpper,
                  k, m, n, nCFunc, nCFuncI, nCFuncJ, nRoots,
                  sStrideI, sStrideIt, sStrideJ, sStrideJt, sStrideM ;
    Real          aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                  c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dGx, dGy, dGz, dHx, dHy, dHz, dnuc,
                  dxIJt, dyIJt, dzIJt, expfac, expN, fac, facN, fac2, f00, qN, rho, rIJ2, scale, tI, u2,
                  xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real         *pGx, *pGy, *pGz, *pHx, *pHy, *pHz, *values = NULL, *work = NULL ;
    Integer      *Ix, *Ixd, *Iy, *Iyd, *Iz, *Izd ;
    Real          ar[3], ari[3], *Cij, *gT, *gX, *gY, *gZ, *hX, *hY, *hZ,
                  Gx  [MAXAMP23       ], Gy  [MAXAMP23]       , Gz  [MAXAMP23       ] ,
                  Sx  [MAXAMP2*MAXAMP2], Sy  [MAXAMP2*MAXAMP2], Sz  [MAXAMP2*MAXAMP2] ,
                  xidG[MAXAMP1*MAXAMP1], yidG[MAXAMP1*MAXAMP1], zidG[MAXAMP1*MAXAMP1] ,
                  xidH[MAXAMP1*MAXAMP1], yidH[MAXAMP1*MAXAMP1], zidH[MAXAMP1*MAXAMP1] ;
    RealArray2D  *iC2S, *jC2S ;
    const Real   *rC ;
    Real         *rN ;
    RysQuadrature roots ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    for ( i = 0 ; i < 3 ; i++ ) { dRi[i] = dRj[i] = 0.0e+00 ; }
    /* . Set pointers. */
    Cij = &rWork[   0] ; gT  = &rWork[  s2] ;
    gX  = &rWork[2*s2] ; gY  = &rWork[3*s2] ; gZ  = &rWork[4*s2] ;
    hX  = &rWork[5*s2] ; hY  = &rWork[6*s2] ; hZ  = &rWork[7*s2] ;
    Ix  = &iWork[   0] ; Iy  = &iWork[  s2] ; Iz  = &iWork[2*s2] ;
    Ixd = &iWork[3*s2] ; Iyd = &iWork[4*s2] ; Izd = &iWork[5*s2] ;
    /* . Loop over the nuclear densities. */
    for ( k = 0 ; k < Coordinates3_Rows ( rNP ) ; k++ )
    {
        if ( _Selected ( selectionN, k ) )
        {
            expN = _GetWidthE ( widthsE, k ) ;
            facN = _GetWidthN ( widthsN, k ) ;
            qN   = - Array1D_Item ( charges, k ) ;
            rN   = Coordinates3_RowPointer ( rNP, k ) ;
            /* . Initialize some accumulators. */
            dGx = dGy = dGz = dHx = dHy = dHz = 0.0e+00 ;
            /* . Outer loop over shells. */
            for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
            {
                iAMMax  = iBasis->shells[iShell].lHigh ;
                iC2S    = iBasis->shells[iShell].c2s ;
                nCFuncI = iBasis->shells[iShell].nCBF     ;
                /* . Inner loop over shells. */
                if ( iIsJ ) jUpper = iShell + 1 ;
                else        jUpper = jBasis->nShells ;
                for ( jShell = 0 ; jShell < jUpper ; jShell++ )
                {
                    jAMMax  = jBasis->shells[jShell].lHigh ;
                    jC2S    = jBasis->shells[jShell].c2s ;
                    nCFuncJ = jBasis->shells[jShell].nCBF     ;
                    /* . Set the diagonal block flag. */
                    isDiagonal = iIsJ && ( iShell == jShell ) ;
                    /* . Get the number of roots. */
                    nRoots = ( iAMMax + jAMMax + 2 ) / 2 + 1 ;
                    /* . Strides. */
                    dStrideJ =   1 ;
                    dStrideI = ( jAMMax + 1 ) * dStrideJ ;
                    sStrideJ =   1 ;
                    sStrideI = ( jAMMax + 2 ) * sStrideJ ;
                    sStrideM = ( iAMMax + 2 ) * sStrideI ;
                    /* . Select the expansion center for the recurrence relations. */
                    if ( iAMMax >= jAMMax )
                    {
                       iAMMaxT   = iAMMax   ;
                       jAMMaxT   = jAMMax   ;
                       sStrideIt = sStrideI ;
                       sStrideJt = sStrideJ ;
                       dxIJt     = xIJ      ;
                       dyIJt     = yIJ      ;
                       dzIJt     = zIJ      ;
                       rC        = rI       ;
                    }
                    else
                    {
                       iAMMaxT   = jAMMax   ;
                       jAMMaxT   = iAMMax   ;
                       sStrideIt = sStrideJ ;
                       sStrideJt = sStrideI ;
                       dxIJt     = - xIJ    ;
                       dyIJt     = - yIJ    ;
                       dzIJt     = - zIJ    ;
                       rC        = rJ       ;
                    }
                    /* . Index arrays. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
   	                ix  = iBasis->shells[iShell].cbfPowX[i] * sStrideI ;
	                iy  = iBasis->shells[iShell].cbfPowY[i] * sStrideI ;
	                iz  = iBasis->shells[iShell].cbfPowZ[i] * sStrideI ;
   	                ixd = iBasis->shells[iShell].cbfPowX[i] * dStrideI ;
	                iyd = iBasis->shells[iShell].cbfPowY[i] * dStrideI ;
	                izd = iBasis->shells[iShell].cbfPowZ[i] * dStrideI ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                        {
	                    Ix [n] = jBasis->shells[jShell].cbfPowX[j] + ix  ;
	                    Iy [n] = jBasis->shells[jShell].cbfPowY[j] + iy  ;
	                    Iz [n] = jBasis->shells[jShell].cbfPowZ[j] + iz  ;
	                    Ixd[n] = jBasis->shells[jShell].cbfPowX[j] + ixd ;
	                    Iyd[n] = jBasis->shells[jShell].cbfPowY[j] + iyd ;
	                    Izd[n] = jBasis->shells[jShell].cbfPowZ[j] + izd ;
                        }
                    }
                    /* . Initialize the integral blocks. */
                    nCFunc = nCFuncI * nCFuncJ ;
                    for ( n = 0 ; n < nCFunc ; n++ ) gX[n] = gY[n] = gZ[n] = hX[n] = hY[n] = hZ[n] = 0.0e+00 ;
                    /* . Outer loop over primitives. */
                    for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
                    {
                        /* . Get some information for the primitive. */
	                ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	                arri = ai * rIJ2 ;
	                for ( i = 0 ; i < 3 ; i++ ) ari[i] = ai * rI[i] ;
                        /* . Inner loop over primitives. */
                        for ( jp = 0 ; jp < jBasis->shells[jShell].nPrimitives ; jp++ )
                        {
                            /* . Get some information for the primitive. */
	                    aj    = jBasis->shells[jShell].primitives[jp].exponent ;
	                    aa    = ai + aj ;
	                    aainv = 1.0e+00 / aa ;
	                    fac   = aj * arri * aainv ;
	                    if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                            expfac = exp ( - fac ) * PI252 * aainv ;
	                    for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rJ[i] ) * aainv ;
                            /* . Calculate some factors. */
                            ab    = aa * expN ;
                            aandb = aa + expN ;
                            rho   = ab / aandb ;
                            dnuc  = expfac * ( facN * qN ) / ( expN * sqrt ( aandb ) ) ;
                            /* . Calculate the rys polynomial roots. */
                            c1x = ( ar[0] - rN[0] ) ;
                            c1y = ( ar[1] - rN[1] ) ;
                            c1z = ( ar[2] - rN[2] ) ;
                            RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = expN * ( rN[0] - rC[0] ) + axac ;
                            c3y  = expN * ( rN[1] - rC[1] ) + ayac ;
                            c3z  = expN * ( rN[2] - rC[2] ) + azac ;
                            c4x  = expN * axac ;
                            c4y  = expN * ayac ;
                            c4z  = expN * azac ;
                            /* . Coefficient array. */
                            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                            {
                                tI = dnuc * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                                for ( j = 0 ; j < nCFuncJ ; j++, n++ ) Cij[n] = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                            }
                            /* . Loop over the roots and construct the subsidiary integrals. */
                            for ( m = 0 ; m < nRoots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expN + u2 ) * fac2 ;
                                xcp00 = u2 * c1x * fac ;
                                ycp00 = u2 * c1y * fac ;
                                zcp00 = u2 * c1z * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                GaussianBasisSubsidiary_f1Cg1  ( iAMMax+jAMMax+2 ,
                                                                 0               ,
                                                                 b00             ,
                                                                 b10             ,
                                                                 bp01            ,
                                                                 f00             ,
                                                                 xc00            ,
                                                                 xcp00           ,
                                                                 yc00            ,
                                                                 ycp00           ,
                                                                 zc00            ,
                                                                 zcp00           ,
                                                                 1               , /* . gStrideIJ. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ) ;
                                /* . Reset S. */
                                for ( i = 0 ; i < sStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = 0.0e+00 ;
                                GaussianBasisSubsidiary_f1Xg2i ( iAMMaxT+1       ,
                                                                 jAMMaxT+1       ,
                                                                 0               ,
                                                                 1               , /* . gStrideIJ. */
                                                                 1               , /* . gStrideK. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ,
                                                                 dxIJt           ,
                                                                 dyIJt           ,
                                                                 dzIJt           ,
                                                                 sStrideIt       ,
                                                                 sStrideJt       ,
                                                                 1               , /* . sStrideK. */
                                                                 Sx              ,
                                                                 Sy              ,
                                                                 Sz              ) ;
                                GaussianBasisSubsidiary_f1Xg2r ( Sx              ,
                                                                 Sy              ,
                                                                 Sz              ,
                                                                 xidG            ,
                                                                 yidG            ,
                                                                 zidG            ,
                                                                 xidH            ,
                                                                 yidH            ,
                                                                 zidH            ,
                                                                 ai              ,
                                                                 aj              ,
                                                                 iAMMax          ,
                                                                 jAMMax          ,
                                                                 0               ,
                                                                 sStrideJ        ,
                                                                 sStrideI        ,
                                                                 dStrideJ        ,
                                                                 dStrideI        ) ;
                                /* . Assemble the integrals. */
                                for ( n = 0 ; n < nCFunc ; n++ )
                                {
                                    gX[n] += Cij[n] * ( xidG[Ixd[n]] * Sy  [Iy [n]] * Sz  [Iz [n]] ) ;
                                    gY[n] += Cij[n] * ( Sx  [Ix [n]] * yidG[Iyd[n]] * Sz  [Iz [n]] ) ;
                                    gZ[n] += Cij[n] * ( Sx  [Ix [n]] * Sy  [Iy [n]] * zidG[Izd[n]] ) ;
                                    hX[n] += Cij[n] * ( xidH[Ixd[n]] * Sy  [Iy [n]] * Sz  [Iz [n]] ) ;
                                    hY[n] += Cij[n] * ( Sx  [Ix [n]] * yidH[Iyd[n]] * Sz  [Iz [n]] ) ;
                                    hZ[n] += Cij[n] * ( Sx  [Ix [n]] * Sy  [Iy [n]] * zidH[Izd[n]] ) ;
                                }
                            } /* . nRoots. */
                        } /* . jP. */
                    } /* . iP. */
                    /* . Transform the integrals. */
                    work   = gT ;
                    values = gX ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pGx = values ;
                    values = gY ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pGy = values ;
                    values = gZ ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pGz = values ;
                    values = hX ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pHx = values ;
                    values = hY ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pHy = values ;
                    values = hZ ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pHz = values ;
                    /* . Get the scale factor. */
                    if ( isDiagonal ) scale = 1.0e+00 ;
                    else              scale = 2.0e+00 ;
                    /* .  Add in the blocks of integrals to the derivatives. */
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
                    {
                        ii = iBasis->shells[iShell].nStart + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                        {
                            jj   = jBasis->shells[jShell].nStart + j ;
                            fac  = scale * Array2D_Item ( dOneIJ, ii, jj ) ;
                            dGx += fac * pGx[n] ;
                            dGy += fac * pGy[n] ;
                            dGz += fac * pGz[n] ;
                            dHx += fac * pHx[n] ;
                            dHy += fac * pHy[n] ;
                            dHz += fac * pHz[n] ;
	                }
                    }
                } /* . jShell. */
            } /* . iShell. */
            /* . Sum in the contributions to the gradients. */
            dRi[0] += dGx ; dRi[1] += dGy ; dRi[2] += dGz ;
            dRj[0] += dHx ; dRj[1] += dHy ; dRj[2] += dHz ;
            Coordinates3_DecrementRow ( gN, k, (dGx+dHx), (dGy+dHy), (dGz+dHz) ) ;
        } /* . Selected. */
    } /* . k. */
}
# undef MAXAMP23

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-nuclear/point integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s2 and real 4 * s2 where s2 = ( maximum shell size )^2. */
# define MAXAMP21 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP1 )
void GaussianBasisIntegrals_f2Cm1V ( const GaussianBasis *iBasis     ,
                                     const Real          *rI         ,
                                     const GaussianBasis *jBasis     ,
                                     const Real          *rJ         ,
                                     const RealArray1D   *charges    ,
                                     const RealArray1D   *widthsE    ,
                                     const RealArray1D   *widthsN    ,
                                     const Coordinates3  *rNP        ,
                                     const Selection     *selectionN ,
                                     const Integer        s2         ,
                                           Integer       *iWork      ,
                                           Real          *rWork      ,
                                           RealArray2D   *integrals  )
{
    Boolean       iIsJ ;
    Integer       i, iAMMax, iAMMaxT, ii, ip, iShell, ix, iy, iz,
                  j, jAMMax, jAMMaxT, jj, jp, jShell,
                  jUpper, k, m, n, nCFunc, nCFuncI, nCFuncJ, nRoots,
                  sStrideI, sStrideIt, sStrideJ, sStrideJt, sStrideM ;
    Real          aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                  c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dnuc, dxIJt, dyIJt, dzIJt,
                  expfac, expN, fac, facN, fac2, f00, qN, rho, rIJ2, tI, u2,
                  xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real         *pG, *values = NULL, *work = NULL ;
    Integer      *Ix, *Iy, *Iz ;
    Real          ar[3], ari[3], *Cij, *Dij, *g, *gT,
                  Gx [MAXAMP21       ], Gy [MAXAMP21       ], Gz[MAXAMP21       ] ,
                  Sx [MAXAMP1*MAXAMP1], Sy [MAXAMP1*MAXAMP1], Sz[MAXAMP1*MAXAMP1] ;
    const Real   *rC ;
    Real         *rN ;
    RysQuadrature roots ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    /* . Set pointers. */
    Cij = &rWork[0] ; Dij = &rWork[s2] ; g  = &rWork[2*s2] ; gT = &rWork[3*s2] ;
    Ix  = &iWork[0] ; Iy  = &iWork[s2] ; Iz = &iWork[2*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iAMMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            /* . Get information about the shell. */
            jAMMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            /* . Get the number of roots. */
            nRoots = ( iAMMax + jAMMax ) / 2 + 1 ;
            /* . Strides. */
            sStrideJ =   1 ;
            sStrideI = ( jAMMax + 1 ) * sStrideJ ;
            sStrideM = ( iAMMax + 1 ) * sStrideI ;
            /* . Select the expansion center for the recurrence relations. */
            if ( iAMMax >= jAMMax )
            {
               iAMMaxT   = iAMMax   ;
               jAMMaxT   = jAMMax   ;
               sStrideIt = sStrideI ;
               sStrideJt = sStrideJ ;
               dxIJt     = xIJ      ;
               dyIJt     = yIJ      ;
               dzIJt     = zIJ      ;
               rC        = rI       ;
            }
            else
            {
               iAMMaxT   = jAMMax   ;
               jAMMaxT   = iAMMax   ;
               sStrideIt = sStrideJ ;
               sStrideJt = sStrideI ;
               dxIJt     = - xIJ    ;
               dyIJt     = - yIJ    ;
               dzIJt     = - zIJ    ;
               rC        = rJ       ;
            }
            /* . Index arrays. */
            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
            {
   	        ix = iBasis->shells[iShell].cbfPowX[i] * sStrideI ;
	        iy = iBasis->shells[iShell].cbfPowY[i] * sStrideI ;
	        iz = iBasis->shells[iShell].cbfPowZ[i] * sStrideI ;
                for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                {
	            Ix[n] = jBasis->shells[jShell].cbfPowX[j] + ix  ;
	            Iy[n] = jBasis->shells[jShell].cbfPowY[j] + iy  ;
	            Iz[n] = jBasis->shells[jShell].cbfPowZ[j] + iz  ;
                }
            }
            /* . Initialize the integral blocks. */
            nCFunc = nCFuncI * nCFuncJ ;
            for ( n = 0 ; n < nCFunc ; n++ ) g[n] = 0.0e+00 ;
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	        arri = ai * rIJ2 ;
	        for ( i = 0 ; i < 3 ; i++ ) ari[i] = ai * rI[i] ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nPrimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj    = jBasis->shells[jShell].primitives[jp].exponent ;
	            aa    = ai + aj ;
	            aainv = 1.0e+00 / aa ;
	            fac   = aj * arri * aainv ;
	            if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    expfac = exp ( - fac ) * PI252 * aainv ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rJ[i] ) * aainv ;
                    /* . Coefficient array. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
                        tI = iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ ) Cij[n] = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                    }
                    /* . Loop over the nuclear densities. */
                    for ( k = 0 ; k < Coordinates3_Rows ( rNP ) ; k++ )
                    {
                        if ( _Selected ( selectionN, k ) )
                        {
                            expN = _GetWidthE ( widthsE, k ) ;
                            facN = _GetWidthN ( widthsN, k ) ;
                            qN   = - Array1D_Item ( charges, k ) ;
                            rN   = Coordinates3_RowPointer ( rNP, k ) ;
                            /* . Calculate some factors. */
                            ab     = aa * expN ;
                            aandb  = aa + expN ;
                            rho    = ab / aandb ;
                            dnuc   = expfac * ( facN * qN ) / ( expN * sqrt ( aandb ) ) ;
                            /* . Calculate the rys polynomial roots. */
                            c1x = ( ar[0] - rN[0] ) ;
                            c1y = ( ar[1] - rN[1] ) ;
                            c1z = ( ar[2] - rN[2] ) ;
                            RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = expN * ( rN[0] - rC[0] ) + axac ;
                            c3y  = expN * ( rN[1] - rC[1] ) + ayac ;
                            c3z  = expN * ( rN[2] - rC[2] ) + azac ;
                            c4x  = expN * axac ;
                            c4y  = expN * ayac ;
                            c4z  = expN * azac ;
                            /* . Scaled coefficient array. */
                            for ( n = 0 ; n < nCFunc ; n++ ) Dij[n] = dnuc * Cij[n] ;
                            /* . Loop over the roots and construct the subsidiary integrals. */
                            for ( m = 0 ; m < nRoots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expN + u2 ) * fac2 ;
                                xcp00 = u2 * c1x * fac ;
                                ycp00 = u2 * c1y * fac ;
                                zcp00 = u2 * c1z * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                GaussianBasisSubsidiary_f1Cg1  ( iAMMax+jAMMax ,
                                                                 0             ,
                                                                 b00           ,
                                                                 b10           ,
                                                                 bp01          ,
                                                                 f00           ,
                                                                 xc00          ,
                                                                 xcp00         ,
                                                                 yc00          ,
                                                                 ycp00         ,
                                                                 zc00          ,
                                                                 zcp00         ,
                                                                 1             , /* . gStrideIJ. */
                                                                 Gx            ,
                                                                 Gy            ,
                                                                 Gz            ) ;
                                /* . Reset S. */
                                for ( i = 0 ; i < sStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = 0.0e+00 ;
                                GaussianBasisSubsidiary_f1Xg2i ( iAMMaxT       ,
                                                                 jAMMaxT       ,
                                                                 0             ,
                                                                 1             , /* . gStrideIJ. */
                                                                 1             , /* . gStrideF. */
                                                                 Gx            ,
                                                                 Gy            ,
                                                                 Gz            ,
                                                                 dxIJt         ,
                                                                 dyIJt         ,
                                                                 dzIJt         ,
                                                                 sStrideIt     ,
                                                                 sStrideJt     ,
                                                                 1             , /* . sStrideF. */
                                                                 Sx            ,
                                                                 Sy            ,
                                                                 Sz            ) ;
                                /* . Assemble the integrals. */
                                for ( n = 0 ; n < nCFunc ; n++ ) g[n] += ( Dij[n] * Sx[Ix[n]] * Sy[Iy[n]] * Sz[Iz[n]] ) ;
                            } /* . nRoots. */
                        } /* . Selected. */
                    } /* . k. */
                } /* . jp. */
            } /* . ip. */
            /* . Transform the integrals. */
            values = g  ;
            work   = gT ;
            GaussianBasisTransform2 ( nCFuncI                    ,
                                      nCFuncJ                    ,
                                      iBasis->shells[iShell].c2s ,
                                      jBasis->shells[jShell].c2s ,
                                      &values                    ,
                                      &work                      ) ;
            pG = values ;
            /* . Save the integrals. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                ii = iBasis->shells[iShell].nStart + i ;
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    jj = jBasis->shells[jShell].nStart + j ;
                    Array2D_Item ( integrals, ii, jj ) = pG[n] ;
	        }
            }
        } /* . jShell. */
    } /* . iShell. */
}
# undef MAXAMP21

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-nuclear/point potentials.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s2 and real 3 * s2 where s2 = ( maximum shell size )^2. */
# define MAXAMP21 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP1 )
void GaussianBasisIntegrals_f2Cp1V ( const GaussianBasis *iBasis     ,
                                     const Real          *rI         ,
                                     const GaussianBasis *jBasis     ,
                                     const Real          *rJ         ,
                                     const RealArray1D   *widthsE    ,
                                     const RealArray1D   *widthsN    ,
                                     const Coordinates3  *rNP        ,
                                     const Selection     *selectionN ,
                                     const RealArray2D   *dOneIJ     ,
                                     const Integer        s2         ,
                                           Integer       *iWork      ,
                                           Real          *rWork      ,
                                           RealArray1D   *potentials )
{
    Boolean       iIsJ, isDiagonal ;
    Integer       i, iAMMax, iAMMaxT, ii, ip, iShell, ix, iy, iz,
                  j, jAMMax, jAMMaxT, jj, jp, jShell,
                  jUpper, k, m, n, nCFunc, nCFuncI, nCFuncJ, nRoots,
                  sStrideI, sStrideIt, sStrideJ, sStrideJt, sStrideM ;
    Real          aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z,
                  dnuc, dxIJt, dyIJt, dzIJt, expfac, expN, fac, facN, fac2, f00, pot, rho, rIJ2, scale, tI, u2,
                  xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real         *pG, *values = NULL, *work = NULL ;
    Integer      *Ix, *Iy, *Iz ;
    Real          ar[3], ari[3], *Cij, *g, *gT,
                  Gx [MAXAMP21       ], Gy[MAXAMP21       ], Gz[MAXAMP21       ] ,
                  Sx [MAXAMP1*MAXAMP1], Sy[MAXAMP1*MAXAMP1], Sz[MAXAMP1*MAXAMP1] ;
    const Real   *rC ;
    Real         *rN ;
    RysQuadrature roots ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    /* . Set pointers. */
    Cij = &rWork[0] ; g  = &rWork[s2] ; gT = &rWork[2*s2] ;
    Ix  = &iWork[0] ; Iy = &iWork[s2] ; Iz = &iWork[2*s2] ;
    /* . Loop over the points. */
    for ( k = 0 ; k < Coordinates3_Rows ( rNP ) ; k++ )
    {
        if ( _Selected ( selectionN, k ) )
        {
            expN = _GetWidthE ( widthsE, k ) ;
            facN = _GetWidthN ( widthsN, k ) ;
            rN   = Coordinates3_RowPointer ( rNP, k ) ;
            pot  = 0.0e+00 ;
            /* . Outer loop over shells. */
            for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
            {
                iAMMax  = iBasis->shells[iShell].lHigh ;
                nCFuncI = iBasis->shells[iShell].nCBF     ;
                /* . Inner loop over shells. */
                if ( iIsJ ) jUpper = iShell + 1 ;
                else        jUpper = jBasis->nShells ;
                for ( jShell = 0 ; jShell < jUpper ; jShell++ )
                {
                    jAMMax  = jBasis->shells[jShell].lHigh ;
                    nCFuncJ = jBasis->shells[jShell].nCBF     ;
                    /* . Set the diagonal block flag. */
                    isDiagonal = iIsJ && ( iShell == jShell ) ;
                    /* . Get the number of roots. */
                    nRoots = ( iAMMax + jAMMax ) / 2 + 1 ;
                    /* . Strides. */
                    sStrideJ =   1 ;
                    sStrideI = ( jAMMax + 1 ) * sStrideJ ;
                    sStrideM = ( iAMMax + 1 ) * sStrideI ;
                    /* . Select the expansion center for the recurrence relations. */
                    if ( iAMMax >= jAMMax )
                    {
                       iAMMaxT   = iAMMax   ;
                       jAMMaxT   = jAMMax   ;
                       sStrideIt = sStrideI ;
                       sStrideJt = sStrideJ ;
                       dxIJt     = xIJ      ;
                       dyIJt     = yIJ      ;
                       dzIJt     = zIJ      ;
                       rC        = rI       ;
                    }
                    else
                    {
                       iAMMaxT   = jAMMax   ;
                       jAMMaxT   = iAMMax   ;
                       sStrideIt = sStrideJ ;
                       sStrideJt = sStrideI ;
                       dxIJt     = - xIJ    ;
                       dyIJt     = - yIJ    ;
                       dzIJt     = - zIJ    ;
                       rC        = rJ       ;
                    }
                    /* . Index arrays. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
   	                ix = iBasis->shells[iShell].cbfPowX[i] * sStrideI ;
	                iy = iBasis->shells[iShell].cbfPowY[i] * sStrideI ;
	                iz = iBasis->shells[iShell].cbfPowZ[i] * sStrideI ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                        {
	                    Ix[n] = jBasis->shells[jShell].cbfPowX[j] + ix  ;
	                    Iy[n] = jBasis->shells[jShell].cbfPowY[j] + iy  ;
	                    Iz[n] = jBasis->shells[jShell].cbfPowZ[j] + iz  ;
                        }
                    }
                    /* . Initialize the integral blocks. */
                    nCFunc = nCFuncI * nCFuncJ ;
                    for ( n = 0 ; n < nCFunc ; n++ ) g[n] = 0.0e+00 ;
                    /* . Outer loop over primitives. */
                    for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
                    {
                        /* . Get some information for the primitive. */
	                ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	                arri = ai * rIJ2 ;
	                for ( i = 0 ; i < 3 ; i++ ) ari[i] = ai * rI[i] ;
                        /* . Inner loop over primitives. */
                        for ( jp = 0 ; jp < jBasis->shells[jShell].nPrimitives ; jp++ )
                        {
                            /* . Get some information for the primitive. */
	                    aj    = jBasis->shells[jShell].primitives[jp].exponent ;
	                    aa    = ai + aj ;
	                    aainv = 1.0e+00 / aa ;
	                    fac   = aj * arri * aainv ;
	                    if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                            expfac = exp ( - fac ) * PI252 * aainv ;
	                    for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rJ[i] ) * aainv ;
                            /* . Start of point-specific code. */
                            /* . Calculate some factors. */
                            ab     = aa * expN ;
                            aandb  = aa + expN ;
                            rho    = ab / aandb ;
                            dnuc   = expfac * facN / ( expN * sqrt ( aandb ) ) ;
                            /* . Calculate the rys polynomial roots. */
                            c1x = ( ar[0] - rN[0] ) ;
                            c1y = ( ar[1] - rN[1] ) ;
                            c1z = ( ar[2] - rN[2] ) ;
                            RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = expN * ( rN[0] - rC[0] ) + axac ;
                            c3y  = expN * ( rN[1] - rC[1] ) + ayac ;
                            c3z  = expN * ( rN[2] - rC[2] ) + azac ;
                            c4x  = expN * axac ;
                            c4y  = expN * ayac ;
                            c4z  = expN * azac ;
                            /* . Coefficient array. */
                            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                            {
                                tI = dnuc * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                                for ( j = 0 ; j < nCFuncJ ; j++, n++ ) Cij[n] = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                            }
                            /* . Loop over the roots and construct the subsidiary integrals. */
                            for ( m = 0 ; m < nRoots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expN + u2 ) * fac2 ;
                                xcp00 = u2 * c1x * fac ;
                                ycp00 = u2 * c1y * fac ;
                                zcp00 = u2 * c1z * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                GaussianBasisSubsidiary_f1Cg1  ( iAMMax+jAMMax ,
                                                                 0             ,
                                                                 b00           ,
                                                                 b10           ,
                                                                 bp01          ,
                                                                 f00           ,
                                                                 xc00          ,
                                                                 xcp00         ,
                                                                 yc00          ,
                                                                 ycp00         ,
                                                                 zc00          ,
                                                                 zcp00         ,
                                                                 1             , /* . gStrideIJ. */
                                                                 Gx            ,
                                                                 Gy            ,
                                                                 Gz            ) ;
                                /* . Reset S. */
                                for ( i = 0 ; i < sStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = 0.0e+00 ;
                                GaussianBasisSubsidiary_f1Xg2i ( iAMMaxT       ,
                                                                 jAMMaxT       ,
                                                                 0             ,
                                                                 1             , /* . gStrideIJ. */
                                                                 1             , /* . gStrideF. */
                                                                 Gx            ,
                                                                 Gy            ,
                                                                 Gz            ,
                                                                 dxIJt         ,
                                                                 dyIJt         ,
                                                                 dzIJt         ,
                                                                 sStrideIt     ,
                                                                 sStrideJt     ,
                                                                 1             , /* . sStrideF. */
                                                                 Sx            ,
                                                                 Sy            ,
                                                                 Sz            ) ;
                                /* . Assemble the integrals. */
                                for ( n = 0 ; n < nCFunc ; n++ ) g[n] += ( Cij[n] * Sx[Ix[n]] * Sy[Iy[n]] * Sz[Iz[n]] ) ;
                            } /* . nRoots. */
                        } /* . jP. */
                    } /* . iP. */
                    /* . Transform the integrals. */
                    values = g  ;
                    work   = gT ;
                    GaussianBasisTransform2 ( nCFuncI                    ,
                                              nCFuncJ                    ,
                                              iBasis->shells[iShell].c2s ,
                                              jBasis->shells[jShell].c2s ,
                                              &values                    ,
                                              &work                      ) ;
                    pG = values ;
                    /* . Determine scaling. */
                    if ( isDiagonal ) scale = 1.0e+00 ;
                    else              scale = 2.0e+00 ;
                    /* .  Add in the block of integrals to the potential - i is usually greater than j. */
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
                    {
                        ii = iBasis->shells[iShell].nStart + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                        {
                            jj = jBasis->shells[jShell].nStart + j ;
                            pot += ( Array2D_Item ( dOneIJ, ii, jj ) * pG[n] * scale ) ;
	                }
                    }
                } /* . jShell. */
            } /* . iShell. */
            /* . Save the potential. */
            Array1D_Item ( potentials, k ) -= pot ;
        } /* . Selected. */
    } /* . k. */
}
# undef MAXAMP21

# undef _Selected
