/*==================================================================================================================================
! . Integrals - 3 basis, 2 electrons, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisIntegrals_f1Xg2.h"
# include "GaussianBasisSubsidiary.h"
# include "GaussianBasisTransform.h"
# include "RysQuadrature.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Options.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . For testing. */
/*# define _UseLocalAccumulators*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-fit anti-Coulomb integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s3 and real 3 * s3 where s3 = ( maximum shell size )^3. */
# define MAXAMP21 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP1 )
# define MAXAMP23 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP3 )
void GaussianBasisIntegrals_f1Ag2i ( const GaussianBasis *iBasis ,
                                     const Real          *rI     ,
                                     const GaussianBasis *jBasis ,
                                     const Real          *rJ     ,
                                     const Real          *rIJ    ,
                                     const Real           rIJ2   ,
                                     const GaussianBasis *fBasis ,
                                     const Real          *rF     ,
                                     const Integer        s3     ,
                                           Integer       *iWork  ,
                                           Real          *rWork  ,
                                           Block         *block  )
{
    Boolean     iIsJ, isDiagonal ;
    Integer     fammax,          fp, fShell,
                iamMax, iamMaxT, ip, iShell,
                jamMax, jamMaxT, jp, jShell,
                jUpper, nCFuncF, nCFuncI, nCFuncJ, nRoots,
                tStrideI, tStrideIt, tStrideJ, tStrideJt, tStrideM ;
    Real        aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dnuc, dxIJt, dyIJt, dzIJt,
                expfac, expf, fac, fac2, f00, rho, u2,
                xc00, xCF , xcp00, yc00, yCF, ycp00, zc00, zCF, zcp00 ;
    Real       *pG, *values = NULL, *work = NULL ;
    Real        ar[3], arI[3], *g, *gT,
                Gx[MAXAMP23*MAXAMP3], Gy[MAXAMP23*MAXAMP3], Gz[MAXAMP23*MAXAMP3] ,
                Hx[MAXAMP21*MAXAMP1], Hy[MAXAMP21*MAXAMP1], Hz[MAXAMP21*MAXAMP1] ,
                Sx[MAXAMP1*MAXAMP1*MAXAMP1], Sy[MAXAMP1*MAXAMP1*MAXAMP1], Sz[MAXAMP1*MAXAMP1*MAXAMP1] ,
                Tx[MAXAMP1*MAXAMP1*MAXAMP1], Ty[MAXAMP1*MAXAMP1*MAXAMP1], Tz[MAXAMP1*MAXAMP1*MAXAMP1] ;
# ifdef _UseLocalAccumulators
    Real       *gLocal ;
# endif
    const Real *rC ;
    /* . Function loops. */
    Integer     f, i, ijTx, ijTy, ijTz,  iTx, iTy, iTz, j, m, n, nCFunc ;
    Integer    *Ix, *Iy, *Iz ;
    Real        tI, tIJ ;
    Real       *Cijf ;
    RysQuadrature roots ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ         = ( iBasis == jBasis ) && ( rI == rJ ) ;
    /* . Set pointers. */
    Cijf = &rWork[0] ; g  = &rWork[s3] ; gT = &rWork[2*s3] ;
    Ix   = &iWork[0] ; Iy = &iWork[s3] ; Iz = &iWork[2*s3] ;
    /* . Triple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        /* . Get information about the shell. */
        iamMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        /* . Set the Upper limit for the JSHELL loops. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            /* . Get information about the shell. */
            jamMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            /* . Set the diagonal block flag. */
            isDiagonal = iIsJ && ( iShell == jShell ) ;
            /* . Select the expansion center for the recurrence relations. */
            if ( iamMax >= jamMax )
            {
               iamMaxT = iamMax ;
               jamMaxT = jamMax ;
               dxIJt   = rIJ[0] ;
               dyIJt   = rIJ[1] ;
               dzIJt   = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iamMaxT = jamMax ;
               jamMaxT = iamMax ;
               dxIJt   = - rIJ[0] ;
               dyIJt   = - rIJ[1] ;
               dzIJt   = - rIJ[2] ;
               rC      = rJ ;
            }
            /* . Loop over fitting shells. */
            for ( fShell = 0 ; fShell < fBasis->nShells ; fShell++ )
            {
                /* . Get information about the shell. */
                fammax  = fBasis->shells[fShell].lHigh ;
                nCFuncF = fBasis->shells[fShell].nCBF     ;
                /* . Displacement. */
                xCF = ( rC[0] - rF[0] ) ;
                yCF = ( rC[1] - rF[1] ) ;
                zCF = ( rC[2] - rF[2] ) ;
                /* . Roots. */
/* . +4! */
                nRoots = ( fammax + iamMax + jamMax + 4 ) / 2 + 1 ;
                /* . Strides - sStrides are the same as tStrides. */
                tStrideJ =   fammax + 1 ;
                tStrideI = ( jamMax + 1 ) * tStrideJ ;
                tStrideM = ( iamMax + 1 ) * tStrideI ;
                if ( iamMax >= jamMax ) { tStrideIt = tStrideI ; tStrideJt = tStrideJ ; }
                else                    { tStrideIt = tStrideJ ; tStrideJt = tStrideI ; }
                /* . Index arrays. */
                for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                {
   	            iTx = iBasis->shells[iShell].cbfPowX[i] * tStrideI ;
	            iTy = iBasis->shells[iShell].cbfPowY[i] * tStrideI ;
	            iTz = iBasis->shells[iShell].cbfPowZ[i] * tStrideI ;
                    for ( j = 0 ; j < nCFuncJ ; j++ )
                    {
	                ijTx = jBasis->shells[jShell].cbfPowX[j] * tStrideJ + iTx ;
	                ijTy = jBasis->shells[jShell].cbfPowY[j] * tStrideJ + iTy ;
	                ijTz = jBasis->shells[jShell].cbfPowZ[j] * tStrideJ + iTz ;
                        for ( f = 0 ; f < nCFuncF ; f++, n++ )
                        {
                            Ix[n] = fBasis->shells[fShell].cbfPowX[f] + ijTx ;
                            Iy[n] = fBasis->shells[fShell].cbfPowY[f] + ijTy ;
                            Iz[n] = fBasis->shells[fShell].cbfPowZ[f] + ijTz ;
                        }
                    }
                }
                /* . Function initialization. */
                nCFunc = nCFuncI * nCFuncJ * nCFuncF ;
                for ( n = 0 ; n < nCFunc ; n++ ) g[n] = 0.0e+00 ;
                /* . Triple loop over primitives. */
                for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
                {
                    /* . Get some information for the primitive. */
	            ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	            arri = ai * rIJ2 ;
	            for ( i = 0 ; i < 3 ; i++ ) arI[i] = ai * rI[i] ;
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
	                for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                        /* . Loop over fitting primitives. */
                        for ( fp = 0 ; fp < fBasis->shells[fShell].nPrimitives ; fp++ )
                        {
                            /* . Get some information for the primitive. */
	                    expf = fBasis->shells[fShell].primitives[fp].exponent ;
                            /* . Calculate some factors. */
                            ab    = aa * expf ;
                            aandb = aa + expf ;
                            rho   = ab / aandb ;
                            dnuc  = expfac / ( expf * sqrt ( aandb ) ) ;
                            /* . Calculate the Rys polynomial roots. */
                            c1x  = ( ar[0] - rF[0] ) ;
                            c1y  = ( ar[1] - rF[1] ) ;
                            c1z  = ( ar[2] - rF[2] ) ;
                            RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = - expf * xCF + axac ;
                            c3y  = - expf * yCF + ayac ;
                            c3z  = - expf * zCF + azac ;
                            c4x  =   expf * axac ;
                            c4y  =   expf * ayac ;
                            c4z  =   expf * azac ;
                            /* . Coefficient array. */
                            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                            {
                                tI = dnuc * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                                for ( j = 0 ; j < nCFuncJ ; j++ )
                                {
                                    tIJ = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                                    for ( f = 0 ; f < nCFuncF ; f++, n++ )
                                    {
                                        Cijf[n] = tIJ * fBasis->shells[fShell].primitives[fp].cCBF[f] ;
                                    }
                                }
                            }
# ifdef _UseLocalAccumulators
                            for ( n = 0 ; n < nCFunc ; n++ ) gLocal[n] = 0.0e+00 ;
# endif
                            /* . Loop over Rys roots. */
                            for ( m = 0 ; m < nRoots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expf + u2 ) * fac2 ;
                                xcp00 =   u2 * c1x         * fac ;
                                ycp00 =   u2 * c1y         * fac ;
                                zcp00 =   u2 * c1z         * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                GaussianBasisSubsidiary_f1Cg1  ( iamMax+jamMax+2 ,
                                                                 fammax+2        ,
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
                                                                 fammax+3        , /* . gStrideIJ. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ) ;
                                GaussianBasisSubsidiary_f1Ag1  ( iamMax+jamMax   ,
                                                                 fammax          ,
                                                                 fammax+3        , /* . gStrideIJ. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ,
                                                                 xCF             ,
                                                                 yCF             ,
                                                                 zCF             ,
                                                                 fammax+1        , /* . hStrideIJ. */
                                                                 Hx              ,
                                                                 Hy              ,
                                                                 Hz              ) ;
                                /* . Reset S and T. */
                                for ( i = 0 ; i < tStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = Tx[i] = Ty[i] = Tz[i] = 0.0e+00 ;
                                GaussianBasisSubsidiary_f1Xg2i ( iamMaxT         ,
                                                                 jamMaxT         ,
                                                                 fammax          ,
                                                                 fammax+3        , /* . gStrideIJ. */
                                                                 1               , /* . gStrideF. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ,
                                                                 dxIJt           ,
                                                                 dyIJt           ,
                                                                 dzIJt           ,
                                                                 tStrideIt       ,
                                                                 tStrideJt       ,
                                                                 1               , /* . tStrideF. */
                                                                 Sx              ,
                                                                 Sy              ,
                                                                 Sz              ) ;
                                GaussianBasisSubsidiary_f1Xg2i ( iamMaxT         ,
                                                                 jamMaxT         ,
                                                                 fammax          ,
                                                                 fammax+1        , /* . hStrideIJ. */
                                                                 1               , /* . hStrideF. */
                                                                 Hx              ,
                                                                 Hy              ,
                                                                 Hz              ,
                                                                 dxIJt           ,
                                                                 dyIJt           ,
                                                                 dzIJt           ,
                                                                 tStrideIt       ,
                                                                 tStrideJt       ,
                                                                 1               , /* . tStrideF. */
                                                                 Tx              ,
                                                                 Ty              ,
                                                                 Tz              ) ;
                                /* . Assemble the integrals. */
# ifdef _UseLocalAccumulators
                                for ( n = 0 ; n < nCFunc ; n++ )
                                {
                                    gLocal[n] += ( Tx[Ix[n]] * Sy[Iy[n]] * Sz[Iz[n]] + 
                                                   Sx[Ix[n]] * Ty[Iy[n]] * Sz[Iz[n]] + 
                                                   Sx[Ix[n]] * Sy[Iy[n]] * Tz[Iz[n]] ) ;
                                }
# else
                                for ( n = 0 ; n < nCFunc ; n++ )
                                {
                                    g[n] += Cijf[n] * ( Tx[Ix[n]] * Sy[Iy[n]] * Sz[Iz[n]] + 
                                                        Sx[Ix[n]] * Ty[Iy[n]] * Sz[Iz[n]] + 
                                                        Sx[Ix[n]] * Sy[Iy[n]] * Tz[Iz[n]] ) ;
                                }
# endif
                            } /* . nRoots. */
# ifdef _UseLocalAccumulators
                            for ( n = 0 ; n < nCFunc ; n++ ) g[n] += Cijf[n] * gLocal[n] ;
# endif
                        } /* . fP. */
                    } /* . jP. */
                } /* . iP. */
                /* . Transform and save the integrals. */
                {
                    auto Boolean     skip ;
                    auto Integer     ii, jj, m, n ;
                    auto Cardinal16 *indices16 = block->indices16 ;
                    auto Real       *integrals = block->data      ;
                    values = g  ;
                    work   = gT ;
                    GaussianBasisTransform3 ( nCFuncI                    ,
                                              nCFuncJ                    ,
                                              nCFuncF                    ,
                                              iBasis->shells[iShell].c2s ,
                                              jBasis->shells[jShell].c2s ,
                                              fBasis->shells[fShell].c2s ,
                                              &values                    ,
                                              &work                      ) ;
                    pG = values ;
                    m  = block->count ;
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
                    {
                        ii = iBasis->shells[iShell].nStart + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++ )
                        {
                            skip = isDiagonal && ( j > i ) ;
                            jj   = jBasis->shells[jShell].nStart + j ;
                            for ( f = 0 ; f < fBasis->shells[fShell].nBasis ; f++, n++ )
                            {
                                if ( ! skip )
                                {
                                    auto Integer m3 = 3 * m ;
                                    indices16[m3  ] = ii ;
                                    indices16[m3+1] = jj ;
                                    indices16[m3+2] = fBasis->shells[fShell].nStart + f ;
                                    integrals[m]    = -pG[n] ; /* . -r12 operator. */
                                    m++ ;
                                }
                            }
                        }
                    }
                    block->count = m ;
                }
            } /* . fShell. */
        } /* . jShell. */
    } /* . iShell. */
}
# undef MAXAMP21
# undef MAXAMP23

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-fit anti-Coulomb integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 6 * s3 and real 8 * s3 where s3 = ( maximum shell size )^3. */
# define MAXAMP23 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP3 )
# define MAXAMP25 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP5 )
void GaussianBasisIntegrals_f1Ag2r1 ( const GaussianBasis *iBasis ,
                                      const Real          *rI     ,
                                      const GaussianBasis *jBasis ,
                                      const Real          *rJ     ,
                                      const Real          *rIJ    ,
                                      const Real           rIJ2   ,
                                      const GaussianBasis *fBasis ,
                                      const Real          *rF     ,
                                      const Integer        s3     ,
                                            Integer       *iWork  ,
                                            Real          *rWork  ,
                                            Block         *block  )
{
    Boolean       iIsJ, isDiagonal ;
    Integer       dStrideI, dStrideJ,
                  fammax,          fp, fShell,
                  iamMax, iamMaxT, ip, iShell,
                  jamMax, jamMaxT, jp, jShell,
                  jUpper, nCFuncF, nCFuncI, nCFuncJ, nRoots,
                  tStrideI, tStrideIt, tStrideJ, tStrideJt, tStrideM ;
    Real          aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                  c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dnuc,
                  dxIJt, dyIJt, dzIJt, expfac, expf, fac, fac2, f00, rho, u2,
                  xc00, xCF, xcp00, yc00, yCF, ycp00, zc00, zCF, zcp00 ;
    Real         *pGx, *pGy, *pGz, *pHx, *pHy, *pHz, *values = NULL, *work = NULL ;
    Real          ar[3], arI[3], *gT, *gX, *gY, *gZ, *hX, *hY, *hZ,
                  Gx[MAXAMP25*MAXAMP3], Gy[MAXAMP25*MAXAMP3], Gz[MAXAMP25*MAXAMP3] ,
                  Hx[MAXAMP23*MAXAMP1], Hy[MAXAMP23*MAXAMP1], Hz[MAXAMP23*MAXAMP1] ,
                  Sx  [MAXAMP2*MAXAMP2*MAXAMP1], Sy  [MAXAMP2*MAXAMP2*MAXAMP1], Sz  [MAXAMP2*MAXAMP2*MAXAMP1] ,
                  SxDg[MAXAMP1*MAXAMP1*MAXAMP1], SyDg[MAXAMP1*MAXAMP1*MAXAMP1], SzDg[MAXAMP1*MAXAMP1*MAXAMP1] ,
                  SxDh[MAXAMP1*MAXAMP1*MAXAMP1], SyDh[MAXAMP1*MAXAMP1*MAXAMP1], SzDh[MAXAMP1*MAXAMP1*MAXAMP1] ,
                  Tx  [MAXAMP2*MAXAMP1*MAXAMP1], Ty  [MAXAMP2*MAXAMP1*MAXAMP1], Tz  [MAXAMP2*MAXAMP1*MAXAMP1] ,
                  TxDg[MAXAMP1*MAXAMP1*MAXAMP1], TyDg[MAXAMP1*MAXAMP1*MAXAMP1], TzDg[MAXAMP1*MAXAMP1*MAXAMP1] ,
                  TxDh[MAXAMP1*MAXAMP1*MAXAMP1], TyDh[MAXAMP1*MAXAMP1*MAXAMP1], TzDh[MAXAMP1*MAXAMP1*MAXAMP1] ;
# ifdef _UseLocalAccumulators
    Real         *gXLocal, *gYLocal, *gZLocal,
                 *hXLocal, *hYLocal, *hZLocal ;
# endif
    /* . Function loops. */
    Integer       f, i, iDx, iDy, iDz, ijDx, ijDy, ijDz, ijTx, ijTy, ijTz, iTx, iTy, iTz, j, m, n, nCFunc ;
    Integer      *IDx, *IDy, *IDz ,
                 *ITx, *ITy, *ITz ;
    Real          tI, tIJ ;
    Real         *Cijf ;
    const Real   *rC ;
    RysQuadrature roots ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ         = ( iBasis == jBasis ) && ( rI == rJ ) ;
    /* . Set pointers. */
    Cijf = &rWork[   0] ;
    gX   = &rWork[  s3] ; gY  = &rWork[2*s3] ; gZ  = &rWork[3*s3] ;
    hX   = &rWork[4*s3] ; hY  = &rWork[5*s3] ; hZ  = &rWork[6*s3] ;
    gT   = &rWork[7*s3] ;
    IDx  = &iWork[   0] ; IDy = &iWork[  s3] ; IDz = &iWork[2*s3] ;
    ITx  = &iWork[3*s3] ; ITy = &iWork[4*s3] ; ITz = &iWork[5*s3] ;
    /* . Triple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        /* . Get information about the shell. */
        iamMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        /* . Set the Upper limit for the JSHELL loops. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            /* . Get information about the shell. */
            jamMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            /* . Set the diagonal block flag. */
            isDiagonal = iIsJ && ( iShell == jShell ) ;
            /* . Select the expansion center for the recurrence relations. */
            if ( iamMax >= jamMax )
            {
               iamMaxT = iamMax ;
               jamMaxT = jamMax ;
               dxIJt   = rIJ[0] ;
               dyIJt   = rIJ[1] ;
               dzIJt   = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iamMaxT = jamMax ;
               jamMaxT = iamMax ;
               dxIJt   = - rIJ[0] ;
               dyIJt   = - rIJ[1] ;
               dzIJt   = - rIJ[2] ;
               rC      = rJ ;
            }
            /* . Loop over fitting shells. */
            for ( fShell = 0 ; fShell < fBasis->nShells ; fShell++ )
            {
                /* . Get information about the shell. */
                fammax  = fBasis->shells[fShell].lHigh ;
                nCFuncF = fBasis->shells[fShell].nCBF     ;
                /* . Displacement. */
                xCF = ( rC[0] - rF[0] ) ;
                yCF = ( rC[1] - rF[1] ) ;
                zCF = ( rC[2] - rF[2] ) ;
                /* . Roots. */
/* . +6! */
                nRoots = ( fammax + iamMax + jamMax + 6 ) / 2 + 1 ;
                /* . Strides. */
                dStrideJ =   fammax + 1 ;
                dStrideI = ( jamMax + 1 ) * dStrideJ ;
                tStrideJ =   fammax + 1 ;
                tStrideI = ( jamMax + 2 ) * tStrideJ ;
                tStrideM = ( iamMax + 2 ) * tStrideI ;
                if ( iamMax >= jamMax ) { tStrideIt = tStrideI ; tStrideJt = tStrideJ ; }
                else                    { tStrideIt = tStrideJ ; tStrideJt = tStrideI ; }
                /* . Index arrays. */
                for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                {
   	            iDx = iBasis->shells[iShell].cbfPowX[i] * dStrideI ;
	            iDy = iBasis->shells[iShell].cbfPowY[i] * dStrideI ;
	            iDz = iBasis->shells[iShell].cbfPowZ[i] * dStrideI ;
   	            iTx = iBasis->shells[iShell].cbfPowX[i] * tStrideI ;
	            iTy = iBasis->shells[iShell].cbfPowY[i] * tStrideI ;
	            iTz = iBasis->shells[iShell].cbfPowZ[i] * tStrideI ;
                    for ( j = 0 ; j < nCFuncJ ; j++ )
                    {
	                ijDx = jBasis->shells[jShell].cbfPowX[j] * dStrideJ + iDx ;
	                ijDy = jBasis->shells[jShell].cbfPowY[j] * dStrideJ + iDy ;
	                ijDz = jBasis->shells[jShell].cbfPowZ[j] * dStrideJ + iDz ;
	                ijTx = jBasis->shells[jShell].cbfPowX[j] * tStrideJ + iTx ;
	                ijTy = jBasis->shells[jShell].cbfPowY[j] * tStrideJ + iTy ;
	                ijTz = jBasis->shells[jShell].cbfPowZ[j] * tStrideJ + iTz ;
                        for ( f = 0 ; f < nCFuncF ; f++, n++ )
                        {
                            IDx[n] = fBasis->shells[fShell].cbfPowX[f] + ijDx ;
                            IDy[n] = fBasis->shells[fShell].cbfPowY[f] + ijDy ;
                            IDz[n] = fBasis->shells[fShell].cbfPowZ[f] + ijDz ;
                            ITx[n] = fBasis->shells[fShell].cbfPowX[f] + ijTx ;
                            ITy[n] = fBasis->shells[fShell].cbfPowY[f] + ijTy ;
                            ITz[n] = fBasis->shells[fShell].cbfPowZ[f] + ijTz ;
                        }
                    }
                }
                /* . Function initialization. */
                nCFunc = nCFuncI * nCFuncJ * nCFuncF ;
                for ( n = 0 ; n < nCFunc ; n++ ) gX[n] = gY[n] = gZ[n] = hX[n] = hY[n] = hZ[n] = 0.0e+00 ;
                /* . Triple loop over primitives. */
                for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
                {
                    /* . Get some information for the primitive. */
	            ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	            arri = ai * rIJ2 ;
	            for ( i = 0 ; i < 3 ; i++ ) arI[i] = ai * rI[i] ;
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
	                for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                        /* . Loop over fitting primitives. */
                        for ( fp = 0 ; fp < fBasis->shells[fShell].nPrimitives ; fp++ )
                        {
                            /* . Get some information for the primitive. */
	                    expf = fBasis->shells[fShell].primitives[fp].exponent ;
                            /* . Calculate some factors. */
                            ab    = aa * expf ;
                            aandb = aa + expf ;
                            rho   = ab / aandb ;
                            dnuc  = expfac / ( expf * sqrt ( aandb ) ) ;
                            /* . Calculate the Rys polynomial roots. */
                            c1x = ( ar[0] - rF[0] ) ;
                            c1y = ( ar[1] - rF[1] ) ;
                            c1z = ( ar[2] - rF[2] ) ;
                            RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = - expf * xCF + axac ;
                            c3y  = - expf * yCF + ayac ;
                            c3z  = - expf * zCF + azac ;
                            c4x  =   expf * axac ;
                            c4y  =   expf * ayac ;
                            c4z  =   expf * azac ;
                            /* . Coefficient array. */
                            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                            {
                                tI = dnuc * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                                for ( j = 0 ; j < nCFuncJ ; j++ )
                                {
                                    tIJ = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                                    for ( f = 0 ; f < nCFuncF ; f++, n++ )
                                    {
                                        Cijf[n] = tIJ * fBasis->shells[fShell].primitives[fp].cCBF[f] ;
                                    }
                                }
                            }
# ifdef _UseLocalAccumulators
                            for ( n = 0 ; n < nCFunc ; n++ ) gXLocal[n] = gYLocal[n] = gZLocal[n] = hXLocal[n] = hYLocal[n] = hZLocal[n] = 0.0e+00 ;
# endif
                            /* . Loop over Rys roots. */
                            for ( m = 0 ; m < nRoots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expf + u2 ) * fac2 ;
                                xcp00 =   u2 * c1x         * fac ;
                                ycp00 =   u2 * c1y         * fac ;
                                zcp00 =   u2 * c1z         * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                GaussianBasisSubsidiary_f1Cg1  ( iamMax+jamMax+4 ,
                                                                 fammax+2        ,
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
                                                                 fammax+3        , /* . gStrideIJ. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ) ;
                                GaussianBasisSubsidiary_f1Ag1  ( iamMax+jamMax+2 ,
                                                                 fammax          ,
                                                                 fammax+3        , /* . gStrideIJ. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ,
                                                                 xCF             ,
                                                                 yCF             ,
                                                                 zCF             ,
                                                                 fammax+1        , /* . hStrideIJ. */
                                                                 Hx              ,
                                                                 Hy              ,
                                                                 Hz              ) ;
                                /* . Reset S and T. */
                                for ( i = 0 ; i < tStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = Tx[i] = Ty[i] = Tz[i] = 0.0e+00 ;
                                GaussianBasisSubsidiary_f1Xg2i ( iamMaxT+1       ,
                                                                 jamMaxT+1       ,
                                                                 fammax          ,
                                                                 fammax+3        , /* . gStrideIJ. */
                                                                 1               , /* . gStrideF. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ,
                                                                 dxIJt           ,
                                                                 dyIJt           ,
                                                                 dzIJt           ,
                                                                 tStrideIt       ,
                                                                 tStrideJt       ,
                                                                 1               , /* . sStrideF. */
                                                                 Sx              ,
                                                                 Sy              ,
                                                                 Sz              ) ;
                                GaussianBasisSubsidiary_f1Xg2i ( iamMaxT+1       ,
                                                                 jamMaxT+1       ,
                                                                 fammax          ,
                                                                 fammax+1        , /* . hStrideIJ. */
                                                                 1               , /* . hStrideF. */
                                                                 Hx              ,
                                                                 Hy              ,
                                                                 Hz              ,
                                                                 dxIJt           ,
                                                                 dyIJt           ,
                                                                 dzIJt           ,
                                                                 tStrideIt       ,
                                                                 tStrideJt       ,
                                                                 1               , /* . tStrideF. */
                                                                 Tx              ,
                                                                 Ty              ,
                                                                 Tz              ) ;
                                GaussianBasisSubsidiary_f1Xg2r ( Sx              ,
                                                                 Sy              ,
                                                                 Sz              ,
                                                                 SxDg            ,
                                                                 SyDg            ,
                                                                 SzDg            ,
                                                                 SxDh            ,
                                                                 SyDh            ,
                                                                 SzDh            ,
                                                                 ai              ,
                                                                 aj              ,
                                                                 iamMax          ,
                                                                 jamMax          ,
                                                                 fammax          ,
                                                                 tStrideJ        ,
                                                                 tStrideI        ,
                                                                 dStrideJ        ,
                                                                 dStrideI        ) ;
                                GaussianBasisSubsidiary_f1Xg2r ( Tx              ,
                                                                 Ty              ,
                                                                 Tz              ,
                                                                 TxDg            ,
                                                                 TyDg            ,
                                                                 TzDg            ,
                                                                 TxDh            ,
                                                                 TyDh            ,
                                                                 TzDh            ,
                                                                 ai              ,
                                                                 aj              ,
                                                                 iamMax          ,
                                                                 jamMax          ,
                                                                 fammax          ,
                                                                 tStrideJ        ,
                                                                 tStrideI        ,
                                                                 dStrideJ        ,
                                                                 dStrideI        ) ;
                                /* . Assemble the integrals. */
# ifdef _UseLocalAccumulators
                                for ( n = 0 ; n < nCFunc ; n++ )
                                {
                                    gXLocal[n] += ( TxDg[IDx[n]] * Sy  [ITy[n]] * Sz  [ITz[n]] + 
                                                    SxDg[IDx[n]] * Ty  [ITy[n]] * Sz  [ITz[n]] + 
                                                    SxDg[IDx[n]] * Sy  [ITy[n]] * Tz  [ITz[n]] ) ;
                                    gYLocal[n] += ( Tx  [ITx[n]] * SyDg[IDy[n]] * Sz  [ITz[n]] + 
                                                    Sx  [ITx[n]] * TyDg[IDy[n]] * Sz  [ITz[n]] + 
                                                    Sx  [ITx[n]] * SyDg[IDy[n]] * Tz  [ITz[n]] ) ;
                                    gZLocal[n] += ( Tx  [ITx[n]] * Sy  [ITy[n]] * SzDg[IDz[n]] + 
                                                    Sx  [ITx[n]] * Ty  [ITy[n]] * SzDg[IDz[n]] + 
                                                    Sx  [ITx[n]] * Sy  [ITy[n]] * TzDg[IDz[n]] ) ;
                                    hXLocal[n] += ( TxDh[IDx[n]] * Sy  [ITy[n]] * Sz  [ITz[n]] + 
                                                    SxDh[IDx[n]] * Ty  [ITy[n]] * Sz  [ITz[n]] + 
                                                    SxDh[IDx[n]] * Sy  [ITy[n]] * Tz  [ITz[n]] ) ;
                                    hYLocal[n] += ( Tx  [ITx[n]] * SyDh[IDy[n]] * Sz  [ITz[n]] + 
                                                    Sx  [ITx[n]] * TyDh[IDy[n]] * Sz  [ITz[n]] + 
                                                    Sx  [ITx[n]] * SyDh[IDy[n]] * Tz  [ITz[n]] ) ;
                                    hZLocal[n] += ( Tx  [ITx[n]] * Sy  [ITy[n]] * SzDh[IDz[n]] + 
                                                    Sx  [ITx[n]] * Ty  [ITy[n]] * SzDh[IDz[n]] + 
                                                    Sx  [ITx[n]] * Sy  [ITy[n]] * TzDh[IDz[n]] ) ;
                                }
# else
                                for ( n = 0 ; n < nCFunc ; n++ )
                                {
                                    gX[n] += Cijf[n] * ( TxDg[IDx[n]] * Sy  [ITy[n]] * Sz  [ITz[n]] + 
                                                         SxDg[IDx[n]] * Ty  [ITy[n]] * Sz  [ITz[n]] + 
                                                         SxDg[IDx[n]] * Sy  [ITy[n]] * Tz  [ITz[n]] ) ;
                                    gY[n] += Cijf[n] * ( Tx  [ITx[n]] * SyDg[IDy[n]] * Sz  [ITz[n]] + 
                                                         Sx  [ITx[n]] * TyDg[IDy[n]] * Sz  [ITz[n]] + 
                                                         Sx  [ITx[n]] * SyDg[IDy[n]] * Tz  [ITz[n]] ) ;
                                    gZ[n] += Cijf[n] * ( Tx  [ITx[n]] * Sy  [ITy[n]] * SzDg[IDz[n]] + 
                                                         Sx  [ITx[n]] * Ty  [ITy[n]] * SzDg[IDz[n]] + 
                                                         Sx  [ITx[n]] * Sy  [ITy[n]] * TzDg[IDz[n]] ) ;
                                    hX[n] += Cijf[n] * ( TxDh[IDx[n]] * Sy  [ITy[n]] * Sz  [ITz[n]] + 
                                                         SxDh[IDx[n]] * Ty  [ITy[n]] * Sz  [ITz[n]] + 
                                                         SxDh[IDx[n]] * Sy  [ITy[n]] * Tz  [ITz[n]] ) ;
                                    hY[n] += Cijf[n] * ( Tx  [ITx[n]] * SyDh[IDy[n]] * Sz  [ITz[n]] + 
                                                         Sx  [ITx[n]] * TyDh[IDy[n]] * Sz  [ITz[n]] + 
                                                         Sx  [ITx[n]] * SyDh[IDy[n]] * Tz  [ITz[n]] ) ;
                                    hZ[n] += Cijf[n] * ( Tx  [ITx[n]] * Sy  [ITy[n]] * SzDh[IDz[n]] + 
                                                         Sx  [ITx[n]] * Ty  [ITy[n]] * SzDh[IDz[n]] + 
                                                         Sx  [ITx[n]] * Sy  [ITy[n]] * TzDh[IDz[n]] ) ;
                                }
# endif
                            } /* . nRoots. */
# ifdef _UseLocalAccumulators
                            for ( n = 0 ; n < nCFunc ; n++ )
                            {
                                 gX[n] += Cijf[n] * gXLocal[n] ;
                                 gY[n] += Cijf[n] * gYLocal[n] ;
                                 gZ[n] += Cijf[n] * gZLocal[n] ;
                                 hX[n] += Cijf[n] * hXLocal[n] ;
                                 hY[n] += Cijf[n] * hYLocal[n] ;
                                 hZ[n] += Cijf[n] * hZLocal[n] ;
                            }
# endif
                        } /* . fP. */
                    } /* . jP. */
                } /* . iP. */
                /* . Transform and save the integrals. */
                {
/*                    auto Boolean     skip ;*/
                    auto Integer      ii, jj, m, m3, m6, n ;
                    auto Cardinal16  *indices16 = block->indices16   ;
                    auto Real        *integrals = block->data, scale ;
                    auto RealArray2D *iC2S = iBasis->shells[iShell].c2s ,
                                     *jC2S = jBasis->shells[jShell].c2s ,
                                     *fC2S = fBasis->shells[fShell].c2s ;
                    work   = gT ;
                    values = gX ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pGx = values ;
                    values = gY ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pGy = values ;
                    values = gZ ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pGz = values ;
                    values = hX ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pHx = values ;
                    values = hY ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pHy = values ;
                    values = hZ ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pHz = values ;
                    if ( isDiagonal ) scale = -1.0e+00 ; /* . -r12 operator. */
                    else              scale = -2.0e+00 ;
                    m = block->count ;
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
                    {
                        ii = iBasis->shells[iShell].nStart + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++ )
                        {
/*                            skip = isDiagonal && ( j > i ) ;*/
                            jj   = jBasis->shells[jShell].nStart + j ;
                            for ( f = 0 ; f < fBasis->shells[fShell].nBasis ; f++, m++, n++ )
                            {
/*
                                if ( ! skip )
                                {
*/
                                    m3              = 3 * m ;
                                    m6              = 6 * m ;
                                    indices16[m3  ] = ii ;
                                    indices16[m3+1] = jj ;
                                    indices16[m3+2] = fBasis->shells[fShell].nStart + f ;
                                    integrals[m6  ] = scale * pGx[n] ;
                                    integrals[m6+1] = scale * pGy[n] ;
                                    integrals[m6+2] = scale * pGz[n] ;
                                    integrals[m6+3] = scale * pHx[n] ;
                                    integrals[m6+4] = scale * pHy[n] ;
                                    integrals[m6+5] = scale * pHz[n] ;
                                    /*m++ ;*/
/*
                                }
*/
                            }
                        }
                    }
                    block->count = m ;
                }
            } /* . fShell. */
        } /* . jShell. */
    } /* . iShell. */
}
# undef MAXAMP23
# undef MAXAMP25

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-fit Coulomb integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s3 and real 3 * s3 where s3 = ( maximum shell size )^3. */
# define MAXAMP21 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP1 )
void GaussianBasisIntegrals_f1Cg2i ( const GaussianBasis *iBasis ,
                                     const Real          *rI     ,
                                     const GaussianBasis *jBasis ,
                                     const Real          *rJ     ,
                                     const Real          *rIJ    ,
                                     const Real           rIJ2   ,
                                     const GaussianBasis *fBasis ,
                                     const Real          *rF     ,
                                     const Integer        s3     ,
                                           Integer       *iWork  ,
                                           Real          *rWork  ,
                                           Block         *block  )
{
    Boolean       iIsJ, isDiagonal ;
    Integer       fammax,          fp, fShell,
                  iamMax, iamMaxT, ip, iShell,
                  jamMax, jamMaxT, jp, jShell,
                  jUpper, nCFuncF, nCFuncI, nCFuncJ, nRoots,
                  sStrideI, sStrideIt, sStrideJ, sStrideJt, sStrideM ;
    Real          aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                  c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dnuc, dxIJt, dyIJt, dzIJt,
                  expfac, expf, fac, fac2, f00, rho, u2,
                  xc00, xcp00, yc00, ycp00, zc00, zcp00 ;
    Real         *pG, *values = NULL, *work = NULL ;
    Real          ar[3], arI[3], *g, *gT,
                  Gx[MAXAMP21*MAXAMP1], Gy[MAXAMP21*MAXAMP1], Gz[MAXAMP21*MAXAMP1] ,
                  Sx[MAXAMP1*MAXAMP1*MAXAMP1], Sy[MAXAMP1*MAXAMP1*MAXAMP1], Sz[MAXAMP1*MAXAMP1*MAXAMP1] ;
    const Real   *rC ;
    /* . Function loops. */
    Integer       f, i, ijTx, ijTy, ijTz,  iTx, iTy, iTz, j, m, n, nCFunc ;
    Integer      *Ix, *Iy, *Iz ;
    Real          tI, tIJ ;
    Real         *Cijf ;
    RysQuadrature roots ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ         = ( iBasis == jBasis ) && ( rI == rJ ) ;
    /* . Set pointers. */
    Cijf = &rWork[0] ; g  = &rWork[s3] ; gT = &rWork[2*s3] ;
    Ix   = &iWork[0] ; Iy = &iWork[s3] ; Iz = &iWork[2*s3] ;
    /* . Triple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        /* . Get information about the shell. */
        iamMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF  ;
        /* . Set the Upper limit for the JSHELL loops. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            /* . Get information about the shell. */
            jamMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF  ;
            /* . Set the diagonal block flag. */
            isDiagonal = iIsJ && ( iShell == jShell ) ;
            /* . Select the expansion center for the recurrence relations. */
            if ( iamMax >= jamMax )
            {
               iamMaxT = iamMax ;
               jamMaxT = jamMax ;
               dxIJt   = rIJ[0] ;
               dyIJt   = rIJ[1] ;
               dzIJt   = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iamMaxT = jamMax ;
               jamMaxT = iamMax ;
               dxIJt   = - rIJ[0] ;
               dyIJt   = - rIJ[1] ;
               dzIJt   = - rIJ[2] ;
               rC      = rJ ;
            }
            /* . Loop over fitting shells. */
            for ( fShell = 0 ; fShell < fBasis->nShells ; fShell++ )
            {
                /* . Get information about the shell. */
                fammax  = fBasis->shells[fShell].lHigh ;
                nCFuncF = fBasis->shells[fShell].nCBF  ;
                /* . Roots. */
                nRoots = ( fammax + iamMax + jamMax ) / 2 + 1 ;
                /* . Strides. */
                sStrideJ =   fammax + 1 ;
                sStrideI = ( jamMax + 1 ) * sStrideJ ;
                sStrideM = ( iamMax + 1 ) * sStrideI ;
                if ( iamMax >= jamMax ) { sStrideIt = sStrideI ; sStrideJt = sStrideJ ; }
                else                    { sStrideIt = sStrideJ ; sStrideJt = sStrideI ; }
                /* . Index arrays. */
                for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                {
   	            iTx = iBasis->shells[iShell].cbfPowX[i] * sStrideI ;
	            iTy = iBasis->shells[iShell].cbfPowY[i] * sStrideI ;
	            iTz = iBasis->shells[iShell].cbfPowZ[i] * sStrideI ;
                    for ( j = 0 ; j < nCFuncJ ; j++ )
                    {
	                ijTx = jBasis->shells[jShell].cbfPowX[j] * sStrideJ + iTx ;
	                ijTy = jBasis->shells[jShell].cbfPowY[j] * sStrideJ + iTy ;
	                ijTz = jBasis->shells[jShell].cbfPowZ[j] * sStrideJ + iTz ;
                        for ( f = 0 ; f < nCFuncF ; f++, n++ )
                        {
                            Ix[n] = fBasis->shells[fShell].cbfPowX[f] + ijTx ;
                            Iy[n] = fBasis->shells[fShell].cbfPowY[f] + ijTy ;
                            Iz[n] = fBasis->shells[fShell].cbfPowZ[f] + ijTz ;
                        }
                    }
                }
                /* . Function initialization. */
                nCFunc = nCFuncI * nCFuncJ * nCFuncF ;
                for ( n = 0 ; n < nCFunc ; n++ ) g[n] = 0.0e+00 ;
                /* . Triple loop over primitives. */
                for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
                {
                    /* . Get some information for the primitive. */
	            ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	            arri = ai * rIJ2 ;
	            for ( i = 0 ; i < 3 ; i++ ) arI[i] = ai * rI[i] ;
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
	                for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                        /* . Loop over fitting primitives. */
                        for ( fp = 0 ; fp < fBasis->shells[fShell].nPrimitives ; fp++ )
                        {
                            /* . Get some information for the primitive. */
	                    expf = fBasis->shells[fShell].primitives[fp].exponent ;
                            /* . Calculate some factors. */
                            ab    = aa * expf ;
                            aandb = aa + expf ;
                            rho   = ab / aandb ;
                            dnuc  = expfac / ( expf * sqrt ( aandb ) ) ;
                            /* . Calculate the Rys polynomial roots. */
                            c1x = ( ar[0] - rF[0] ) ;
                            c1y = ( ar[1] - rF[1] ) ;
                            c1z = ( ar[2] - rF[2] ) ;
                            RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c3x  = expf * ( rF[0] - rC[0] ) + axac ;
                            c3y  = expf * ( rF[1] - rC[1] ) + ayac ;
                            c3z  = expf * ( rF[2] - rC[2] ) + azac ;
                            c4x  = expf * axac ;
                            c4y  = expf * ayac ;
                            c4z  = expf * azac ;
                            /* . Coefficient array. */
                            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                            {
                                tI = dnuc * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                                for ( j = 0 ; j < nCFuncJ ; j++ )
                                {
                                    tIJ = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                                    for ( f = 0 ; f < nCFuncF ; f++, n++ )
                                    {
                                        Cijf[n] = tIJ * fBasis->shells[fShell].primitives[fp].cCBF[f] ;
                                    }
                                }
                            }
                            /* . Loop over Rys roots. */
                            for ( m = 0 ; m < nRoots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expf + u2 ) * fac2 ;
                                xcp00 = aa * u2 * c1x * fac ;
                                ycp00 = aa * u2 * c1y * fac ;
                                zcp00 = aa * u2 * c1z * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                GaussianBasisSubsidiary_f1Cg1  ( iamMax+jamMax   ,
                                                                 fammax          ,
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
                                                                 fammax+1        , /* . gStrideIJ. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ) ;
                                /* . Reset S. */
                                for ( i = 0 ; i < sStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = 0.0e+00 ;
                                GaussianBasisSubsidiary_f1Xg2i ( iamMaxT         ,
                                                                 jamMaxT         ,
                                                                 fammax          ,
                                                                 fammax+1        , /* . gStrideIJ. */
                                                                 1               , /* . gStrideF. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ,
                                                                 dxIJt           ,
                                                                 dyIJt           ,
                                                                 dzIJt           ,
                                                                 sStrideIt       ,
                                                                 sStrideJt       ,
                                                                 1               , /* . sStrideF. */
                                                                 Sx              ,
                                                                 Sy              ,
                                                                 Sz              ) ;
                                /* . Assemble the integrals. */
                                for ( n = 0 ; n < nCFunc ; n++ ) g[n] += ( Cijf[n] * Sx[Ix[n]] * Sy[Iy[n]] * Sz[Iz[n]] ) ;
                            } /* . nRoots. */
                        } /* . fP. */
                    } /* . jP. */
                } /* . iP. */
                /* . Transform and save the integrals. */
                {
                    auto Boolean     skip ;
                    auto Integer     ii, jj, m, n ;
                    auto Cardinal16 *indices16 = block->indices16 ;
                    auto Real       *integrals = block->data      ;
                    values = g  ;
                    work   = gT ;
                    GaussianBasisTransform3 ( nCFuncI                    ,
                                              nCFuncJ                    ,
                                              nCFuncF                    ,
                                              iBasis->shells[iShell].c2s ,
                                              jBasis->shells[jShell].c2s ,
                                              fBasis->shells[fShell].c2s ,
                                              &values                    ,
                                              &work                      ) ;
                    pG = values ;
                    m  = block->count ;
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
                    {
                        ii = iBasis->shells[iShell].nStart + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++ )
                        {
                            skip = isDiagonal && ( j > i ) ;
                            jj   = jBasis->shells[jShell].nStart + j ;
                            for ( f = 0 ; f < fBasis->shells[fShell].nBasis ; f++, n++ )
                            {
                                if ( ! skip )
                                {
                                    auto Integer  m3 = 3 * m ;
                                    indices16[m3  ] = ii ;
                                    indices16[m3+1] = jj ;
                                    indices16[m3+2] = fBasis->shells[fShell].nStart + f ;
                                    integrals[m]    = pG[n] ;
                                    m++ ;
                                }
                            }
                        }
                    }
                    block->count = m ;
                }
            } /* . fShell. */
        } /* . jShell. */
    } /* . iShell. */
}
# undef MAXAMP21

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-fit overlap integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 6 * s3 and real 8 * s3 where s3 = ( maximum shell size )^3. */
# define MAXAMP23 ( MAXIMUM_ANGULAR_MOMENTUM + MAXAMP3 )
void GaussianBasisIntegrals_f1Cg2r1 ( const GaussianBasis *iBasis ,
                                      const Real          *rI     ,
                                      const GaussianBasis *jBasis ,
                                      const Real          *rJ     ,
                                      const Real          *rIJ    ,
                                      const Real           rIJ2   ,
                                      const GaussianBasis *fBasis ,
                                      const Real          *rF     ,
                                      const Integer        s3     ,
                                            Integer       *iWork  ,
                                            Real          *rWork  ,
                                            Block         *block  )
{
    Boolean       iIsJ, isDiagonal ;
    Integer       dStrideI, dStrideJ,
                  fammax,          fp, fShell,
                  iamMax, iamMaxT, ip, iShell,
                  jamMax, jamMaxT, jp, jShell,
                  jUpper, nCFuncF, nCFuncI, nCFuncJ, nRoots,
                  tStrideI, tStrideIt, tStrideJ, tStrideJt, tStrideM ;
    Real          aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                  c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dnuc,
                  dxIJt, dyIJt, dzIJt, expfac, expf, fac, fac2, f00, rho, u2,
                  xc00, xCF, xcp00, yc00, yCF, ycp00, zc00, zCF, zcp00 ;
    Real         *pGx, *pGy, *pGz, *pHx, *pHy, *pHz, *values = NULL, *work = NULL ;
    Real          ar[3], arI[3], *gT, *gX, *gY, *gZ, *hX, *hY, *hZ,
                  Gx  [MAXAMP23*MAXAMP1], Gy[MAXAMP23*MAXAMP1], Gz[MAXAMP23*MAXAMP1] ,
                  Sx  [MAXAMP2*MAXAMP2*MAXAMP1], Sy  [MAXAMP2*MAXAMP2*MAXAMP1], Sz  [MAXAMP2*MAXAMP2*MAXAMP1] ,
                  SxDg[MAXAMP1*MAXAMP1*MAXAMP1], SyDg[MAXAMP1*MAXAMP1*MAXAMP1], SzDg[MAXAMP1*MAXAMP1*MAXAMP1] ,
                  SxDh[MAXAMP1*MAXAMP1*MAXAMP1], SyDh[MAXAMP1*MAXAMP1*MAXAMP1], SzDh[MAXAMP1*MAXAMP1*MAXAMP1] ;
    const Real   *rC ;
    /* . Function loops. */
    Integer       f, i, iDx, iDy, iDz, ijDx, ijDy, ijDz, ijTx, ijTy, ijTz, iTx, iTy, iTz, j, m, n, nCFunc ;
    Integer      *IDx, *IDy, *IDz ,
                 *ITx, *ITy, *ITz ;
    Real          tI, tIJ ;
    Real         *Cijf ;
    RysQuadrature roots ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ         = ( iBasis == jBasis ) && ( rI == rJ ) ;
    /* . Set pointers. */
    Cijf = &rWork[   0] ;
    gX   = &rWork[  s3] ; gY  = &rWork[2*s3] ; gZ  = &rWork[3*s3] ;
    hX   = &rWork[4*s3] ; hY  = &rWork[5*s3] ; hZ  = &rWork[6*s3] ;
    gT   = &rWork[7*s3] ;
    IDx  = &iWork[   0] ; IDy = &iWork[  s3] ; IDz = &iWork[2*s3] ;
    ITx  = &iWork[3*s3] ; ITy = &iWork[4*s3] ; ITz = &iWork[5*s3] ;
    /* . Triple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        /* . Get information about the shell. */
        iamMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF  ;
        /* . Set the Upper limit for the JSHELL loops. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            /* . Get information about the shell. */
            jamMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF  ;
            /* . Set the diagonal block flag. */
            isDiagonal = iIsJ && ( iShell == jShell ) ;
            /* . Select the expansion center for the recurrence relations. */
            if ( iamMax >= jamMax )
            {
               iamMaxT = iamMax ;
               jamMaxT = jamMax ;
               dxIJt   = rIJ[0] ;
               dyIJt   = rIJ[1] ;
               dzIJt   = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iamMaxT = jamMax ;
               jamMaxT = iamMax ;
               dxIJt   = - rIJ[0] ;
               dyIJt   = - rIJ[1] ;
               dzIJt   = - rIJ[2] ;
               rC      = rJ ;
            }
            /* . Loop over fitting shells. */
            for ( fShell = 0 ; fShell < fBasis->nShells ; fShell++ )
            {
                /* . Get information about the shell. */
                fammax  = fBasis->shells[fShell].lHigh ;
                nCFuncF = fBasis->shells[fShell].nCBF  ;
                /* . Displacement. */
                xCF = ( rC[0] - rF[0] ) ;
                yCF = ( rC[1] - rF[1] ) ;
                zCF = ( rC[2] - rF[2] ) ;
                /* . Roots. */
                nRoots = ( fammax + iamMax + jamMax + 2 ) / 2 + 1 ;
                /* . Index arrays. */
                /* . Strides. */
                dStrideJ =   fammax + 1 ;
                dStrideI = ( jamMax + 1 ) * dStrideJ ;
                tStrideJ =   fammax + 1 ;
                tStrideI = ( jamMax + 2 ) * tStrideJ ;
                tStrideM = ( iamMax + 2 ) * tStrideI ;
                if ( iamMax >= jamMax ) { tStrideIt = tStrideI ; tStrideJt = tStrideJ ; }
                else                    { tStrideIt = tStrideJ ; tStrideJt = tStrideI ; }
                /* . Index arrays. */
                for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                {
   	            iDx = iBasis->shells[iShell].cbfPowX[i] * dStrideI ;
	            iDy = iBasis->shells[iShell].cbfPowY[i] * dStrideI ;
	            iDz = iBasis->shells[iShell].cbfPowZ[i] * dStrideI ;
   	            iTx = iBasis->shells[iShell].cbfPowX[i] * tStrideI ;
	            iTy = iBasis->shells[iShell].cbfPowY[i] * tStrideI ;
	            iTz = iBasis->shells[iShell].cbfPowZ[i] * tStrideI ;
                    for ( j = 0 ; j < nCFuncJ ; j++ )
                    {
	                ijDx = jBasis->shells[jShell].cbfPowX[j] * dStrideJ + iDx ;
	                ijDy = jBasis->shells[jShell].cbfPowY[j] * dStrideJ + iDy ;
	                ijDz = jBasis->shells[jShell].cbfPowZ[j] * dStrideJ + iDz ;
	                ijTx = jBasis->shells[jShell].cbfPowX[j] * tStrideJ + iTx ;
	                ijTy = jBasis->shells[jShell].cbfPowY[j] * tStrideJ + iTy ;
	                ijTz = jBasis->shells[jShell].cbfPowZ[j] * tStrideJ + iTz ;
                        for ( f = 0 ; f < nCFuncF ; f++, n++ )
                        {
                            IDx[n] = fBasis->shells[fShell].cbfPowX[f] + ijDx ;
                            IDy[n] = fBasis->shells[fShell].cbfPowY[f] + ijDy ;
                            IDz[n] = fBasis->shells[fShell].cbfPowZ[f] + ijDz ;
                            ITx[n] = fBasis->shells[fShell].cbfPowX[f] + ijTx ;
                            ITy[n] = fBasis->shells[fShell].cbfPowY[f] + ijTy ;
                            ITz[n] = fBasis->shells[fShell].cbfPowZ[f] + ijTz ;
                        }
                    }
                }
                /* . Function initialization. */
                nCFunc = nCFuncI * nCFuncJ * nCFuncF ;
                for ( n = 0 ; n < nCFunc ; n++ ) gX[n] = gY[n] = gZ[n] = hX[n] = hY[n] = hZ[n] = 0.0e+00 ;
                /* . Triple loop over primitives. */
                for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
                {
                    /* . Get some information for the primitive. */
	            ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	            arri = ai * rIJ2 ;
	            for ( i = 0 ; i < 3 ; i++ ) arI[i] = ai * rI[i] ;
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
	                for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                        /* . Loop over fitting primitives. */
                        for ( fp = 0 ; fp < fBasis->shells[fShell].nPrimitives ; fp++ )
                        {
                            /* . Get some information for the primitive. */
	                    expf = fBasis->shells[fShell].primitives[fp].exponent ;
                            /* . Calculate some factors. */
                            ab    = aa * expf ;
                            aandb = aa + expf ;
                            rho   = ab / aandb ;
                            dnuc  = expfac / ( expf * sqrt ( aandb ) ) ;
                            /* . Calculate the Rys polynomial roots. */
                            c1x = ( ar[0] - rF[0] ) ;
                            c1y = ( ar[1] - rF[1] ) ;
                            c1z = ( ar[2] - rF[2] ) ;
                            RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = - expf * xCF + axac ;
                            c3y  = - expf * yCF + ayac ;
                            c3z  = - expf * zCF + azac ;
                            c4x  =   expf * axac ;
                            c4y  =   expf * ayac ;
                            c4z  =   expf * azac ;
                            /* . Coefficient array. */
                            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                            {
                                tI = dnuc * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                                for ( j = 0 ; j < nCFuncJ ; j++ )
                                {
                                    tIJ = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                                    for ( f = 0 ; f < nCFuncF ; f++, n++ )
                                    {
                                        Cijf[n] = tIJ * fBasis->shells[fShell].primitives[fp].cCBF[f] ;
                                    }
                                }
                            }
                            /* . Loop over Rys roots. */
                            for ( m = 0 ; m < nRoots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expf + u2 ) * fac2 ;
                                xcp00 =   u2 * c1x         * fac ;
                                ycp00 =   u2 * c1y         * fac ;
                                zcp00 =   u2 * c1z         * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                GaussianBasisSubsidiary_f1Cg1  ( iamMax+jamMax+2 ,
                                                                 fammax          ,
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
                                                                 fammax+1        , /* . gStrideIJ. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ) ;
                                /* . Reset S. */
                                for ( i = 0 ; i < tStrideM ; i++ ) Sx[i] = Sy[i] = Sz[i] = 0.0e+00 ;
                                GaussianBasisSubsidiary_f1Xg2i ( iamMaxT+1       ,
                                                                 jamMaxT+1       ,
                                                                 fammax          ,
                                                                 fammax+1        , /* . gStrideIJ. */
                                                                 1               , /* . gStrideF. */
                                                                 Gx              ,
                                                                 Gy              ,
                                                                 Gz              ,
                                                                 dxIJt           ,
                                                                 dyIJt           ,
                                                                 dzIJt           ,
                                                                 tStrideIt       ,
                                                                 tStrideJt       ,
                                                                 1               , /* . sStrideF. */
                                                                 Sx              ,
                                                                 Sy              ,
                                                                 Sz              ) ;
                                GaussianBasisSubsidiary_f1Xg2r ( Sx              ,
                                                                 Sy              ,
                                                                 Sz              ,
                                                                 SxDg            ,
                                                                 SyDg            ,
                                                                 SzDg            ,
                                                                 SxDh            ,
                                                                 SyDh            ,
                                                                 SzDh            ,
                                                                 ai              ,
                                                                 aj              ,
                                                                 iamMax          ,
                                                                 jamMax          ,
                                                                 fammax          ,
                                                                 tStrideJ        ,
                                                                 tStrideI        ,
                                                                 dStrideJ        ,
                                                                 dStrideI        ) ;
                                /* . Assemble the integrals. */
                                for ( n = 0 ; n < nCFunc ; n++ )
                                {
                                    gX[n] += Cijf[n] * ( SxDg[IDx[n]] * Sy  [ITy[n]] * Sz  [ITz[n]] ) ;
                                    gY[n] += Cijf[n] * ( Sx  [ITx[n]] * SyDg[IDy[n]] * Sz  [ITz[n]] ) ;
                                    gZ[n] += Cijf[n] * ( Sx  [ITx[n]] * Sy  [ITy[n]] * SzDg[IDz[n]] ) ;
                                    hX[n] += Cijf[n] * ( SxDh[IDx[n]] * Sy  [ITy[n]] * Sz  [ITz[n]] ) ;
                                    hY[n] += Cijf[n] * ( Sx  [ITx[n]] * SyDh[IDy[n]] * Sz  [ITz[n]] ) ;
                                    hZ[n] += Cijf[n] * ( Sx  [ITx[n]] * Sy  [ITy[n]] * SzDh[IDz[n]] ) ;
                                }
                            } /* . nRoots. */
                        } /* . fP. */
                    } /* . jP. */
                } /* . iP. */
                /* . Transform and save the integrals. */
                {
/*                    auto Boolean     skip ;*/
                    auto Integer      ii, jj, m, m3, m6, n ;
                    auto Cardinal16  *indices16 = block->indices16   ;
                    auto Real        *integrals = block->data, scale ;
                    auto RealArray2D *iC2S = iBasis->shells[iShell].c2s ,
                                     *jC2S = jBasis->shells[jShell].c2s ,
                                     *fC2S = fBasis->shells[fShell].c2s ;
                    work   = gT ;
                    values = gX ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pGx = values ;
                    values = gY ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pGy = values ;
                    values = gZ ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pGz = values ;
                    values = hX ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pHx = values ;
                    values = hY ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pHy = values ;
                    values = hZ ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pHz = values ;
                    if ( isDiagonal ) scale = 1.0e+00 ;
                    else              scale = 2.0e+00 ;
                    m = block->count ;
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
                    {
                        ii = iBasis->shells[iShell].nStart + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++ )
                        {
/*                            skip = isDiagonal && ( j > i ) ;*/
                            jj   = jBasis->shells[jShell].nStart + j ;
                            for ( f = 0 ; f < fBasis->shells[fShell].nBasis ; f++, m++, n++ )
                            {
/*
                                if ( ! skip )
                                {
*/
                                    m3              = 3 * m ;
                                    m6              = 6 * m ;
                                    indices16[m3  ] = ii ;
                                    indices16[m3+1] = jj ;
                                    indices16[m3+2] = fBasis->shells[fShell].nStart + f ;
                                    integrals[m6  ] = scale * pGx[n] ;
                                    integrals[m6+1] = scale * pGy[n] ;
                                    integrals[m6+2] = scale * pGz[n] ;
                                    integrals[m6+3] = scale * pHx[n] ;
                                    integrals[m6+4] = scale * pHy[n] ;
                                    integrals[m6+5] = scale * pHz[n] ;
                                    /*m++ ;*/
/*
                                }
*/
                            }
                        }
                    }
                    block->count = m ;
                }
            } /* . fShell. */
        } /* . jShell. */
    } /* . iShell. */
}
# undef MAXAMP23

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-fit overlap integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 2 * s3 where s3 = ( maximum shell size )^3. */
void GaussianBasisIntegrals_f1Og2i ( const GaussianBasis *iBasis ,
                                     const Real          *rI     ,
                                     const GaussianBasis *jBasis ,
                                     const Real          *rJ     ,
                                     const GaussianBasis *fBasis ,
                                     const Real          *rF     ,
                                     const Integer        s3     ,
                                           Real          *rWork  ,
                                           Block         *block  )
{
    Boolean iIsJ, isDiagonal ;
    Integer c, dim1, dim2,
            f, fammax, fijx, fijy, fijz, fp, fShell,
            i, iamMax, ip, iShell, ix, iy, iz,
            j, jamMax, jix, jiy, jiz, jp, jShell,
            jUpper, n, nCFuncF, nCFuncI, nCFuncJ ;
    Real    aI, aIJ, aIJF, aJ, aF, eIJ, eIJF, expfac, rIF2, rIJ2, rJF2, ti, tij ;
    Real   *pG, *values = NULL, *work = NULL ;
    Real    cI[3], cIJ[3], cIJF[3], *g, *gT,
            xint[MAXAMP1*MAXAMP1*MAXAMP1], yint[MAXAMP1*MAXAMP1*MAXAMP1], zint[MAXAMP1*MAXAMP1*MAXAMP1] ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ         = ( iBasis == jBasis ) && ( rI == rJ ) ;
    for ( c = 0, rIJ2 = rIF2 = rJF2 = 0.0e+00 ; c < 3 ; c++ )
    {
        rIJ2 += pow ( rI[c] - rJ[c], 2 ) ;
        rIF2 += pow ( rI[c] - rF[c], 2 ) ;
        rJF2 += pow ( rJ[c] - rF[c], 2 ) ;
    }
    /* . Set pointers. */
    g = &rWork[0] ; gT = &rWork[s3] ;
    /* . Triple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        /* . Get information about the shell. */
        iamMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF  ;
        /* . Set the Upper limit for the JSHELL loops. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            /* . Get information about the shell. */
            jamMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF  ;
            /* . Set the diagonal block flag. */
            isDiagonal = iIsJ && ( iShell == jShell ) ;
            /* . Loop over fitting shells. */
            for ( fShell = 0 ; fShell < fBasis->nShells ; fShell++ )
            {
                /* . Get information about the shell. */
                fammax  = fBasis->shells[fShell].lHigh ;
                nCFuncF = fBasis->shells[fShell].nCBF  ;
                /* . Initialize the integral blocks. */
                for ( i = 0 ; i < ( nCFuncF * nCFuncI * nCFuncJ ) ; i++ ) g[i] = 0.0e+00 ;
                /* . Set the array dimensions. */
                dim1 = fammax + 1 ;
                dim2 = dim1 * ( jamMax + 1 ) ;
/*
printf ( "\nShells i %u %d %d %d %d %d %d %d %d | j %u %d %d %d %d %d %d %d %d | f %u %d %d %d %d %d %d %d %d:\n",
 iBasis->isSpherical, iBasis->atomicNumber, iShell, iamMax, nCFuncI, iBasis->shells[iShell].nBasis, iBasis->shells[iShell].nCBF, iBasis->shells[iShell].nStart, iBasis->shells[iShell].nStartC ,
 jBasis->isSpherical, jBasis->atomicNumber, jShell, jamMax, nCFuncJ, jBasis->shells[jShell].nBasis, jBasis->shells[jShell].nCBF, jBasis->shells[jShell].nStart, jBasis->shells[jShell].nStartC ,
 fBasis->isSpherical, fBasis->atomicNumber, fShell, fammax, nCFuncF, fBasis->shells[fShell].nBasis, fBasis->shells[fShell].nCBF, fBasis->shells[fShell].nStart, fBasis->shells[fShell].nStartC ) ;
*/
                /* . Triple loop over primitives. */
                for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
                {
                    /* . Get some information for the primitive. */
	            aI = iBasis->shells[iShell].primitives[ip].exponent ;
	            for ( i = 0 ; i < 3 ; i++ ) cI[i] = aI * rI[i] ;
                    /* . Middle loop over primitives. */
                    for ( jp = 0 ; jp < jBasis->shells[jShell].nPrimitives ; jp++ )
                    {
                        /* . Get some information for the primitive. */
	                aJ  = jBasis->shells[jShell].primitives[jp].exponent ;
                        aIJ = aI + aJ ;
	                eIJ = aI * aJ * rIJ2 ;
	                if ( ( eIJ / aIJ ) > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
	                for ( i = 0 ; i < 3 ; i++ ) cIJ[i] = ( cI[i] + aJ * rJ[i] ) ;
                        /* . Inner loop over primitives. */
                        for ( fp = 0 ; fp < fBasis->shells[fShell].nPrimitives ; fp++ )
                        {
                            /* . Get some information for the primitive. */
	                    aF   = fBasis->shells[fShell].primitives[fp].exponent ;
                            aIJF = aIJ + aF ;
                            eIJF = ( eIJ + aI * aF * rIF2 + aJ * aF * rJF2 ) / aIJF ;
	                    if ( eIJF > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                            expfac = exp ( - eIJF ) ;
	                    for ( i = 0 ; i < 3 ; i++ ) cIJF[i] = ( cIJ[i] + aF * rF[i] ) / aIJF ;
                            /* . Calculate the overlap integrals. */
                            GaussianBasisSubsidiary_f1Og2 ( xint, yint, zint, aIJF, cIJF, rI, rJ, rF, iamMax, jamMax, fammax ) ;
                            /* . Assemble the integrals. */
                            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                            {
   	                        ix = iBasis->shells[iShell].cbfPowX[i] * dim2 ;
	                        iy = iBasis->shells[iShell].cbfPowY[i] * dim2 ;
	                        iz = iBasis->shells[iShell].cbfPowZ[i] * dim2 ;
                                ti = expfac * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                                for ( j = 0 ; j < nCFuncJ ; j++ )
                                {
	                            jix = jBasis->shells[jShell].cbfPowX[j] * dim1 + ix ;
	                            jiy = jBasis->shells[jShell].cbfPowY[j] * dim1 + iy ;
	                            jiz = jBasis->shells[jShell].cbfPowZ[j] * dim1 + iz ;
                                    tij = ti * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                                    for ( f = 0 ; f < nCFuncF ; f++, n++ )
                                    {
                                        fijx  = fBasis->shells[fShell].cbfPowX[f] + jix ;
                                        fijy  = fBasis->shells[fShell].cbfPowY[f] + jiy ;
                                        fijz  = fBasis->shells[fShell].cbfPowZ[f] + jiz ;
                                        g[n] += ( tij * fBasis->shells[fShell].primitives[fp].cCBF[f] * xint[fijx] * yint[fijy] * zint[fijz] ) ;
                                    }
                                }
                            }
                        } /* . fP. */
                    } /* . jP. */
                } /* . iP. */
                /* . Transform and save the integrals. */
                {
                    auto Boolean     skip ;
                    auto Integer     ii, jj, m, n ;
                    auto Cardinal16 *indices16 = block->indices16 ;
                    auto Real       *integrals = block->data      ;
                    values = g  ;
                    work   = gT ;
                    GaussianBasisTransform3 ( nCFuncI                    ,
                                              nCFuncJ                    ,
                                              nCFuncF                    ,
                                              iBasis->shells[iShell].c2s ,
                                              jBasis->shells[jShell].c2s ,
                                              fBasis->shells[fShell].c2s ,
                                              &values                    ,
                                              &work                      ) ;
                    pG = values ;
                    m  = block->count ;
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
                    {
                        ii = iBasis->shells[iShell].nStart + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++ )
                        {
                            skip = isDiagonal && ( j > i ) ;
                            jj   = jBasis->shells[jShell].nStart + j ;
                            for ( f = 0 ; f < fBasis->shells[fShell].nBasis ; f++, n++ )
                            {
                                if ( ! skip )
                                {
                                    auto Integer m3 = 3 * m ;
                                    indices16[m3  ] = ii ;
                                    indices16[m3+1] = jj ;
                                    indices16[m3+2] = fBasis->shells[fShell].nStart + f ;
                                    integrals[m]    = pG[n] ;
                                    m++ ;
                                }
                            }
                        }
                    }
                    block->count = m ;
                }
            } /* . fShell. */
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-fit overlap integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 7 * s3 where s3 = ( maximum shell size )^3. */
void GaussianBasisIntegrals_f1Og2r1 ( const GaussianBasis *iBasis ,
                                      const Real          *rI     ,
                                      const GaussianBasis *jBasis ,
                                      const Real          *rJ     ,
                                      const GaussianBasis *fBasis ,
                                      const Real          *rF     ,
                                      const Integer        s3     ,
                                            Real          *rWork  ,
                                            Block         *block  )
{
    Boolean iIsJ, isDiagonal ;
    Integer c, dim1, ddim2, dim2,
            f, fammax, fijx, fijxd, fijy, fijyd, fijz, fijzd, fp, fShell,
            i, iamMax, ip, iShell, ix, ixd, iy, iyd, iz, izd,
            j, jamMax, jix, jixd, jiy, jiyd, jiz, jizd, jp, jShell,
            jUpper, n, nCFuncF, nCFuncI, nCFuncJ ;
    Real    aI, aIJ, aIJF, aJ, aF, eIJ, eIJF, expfac, rIF2, rIJ2, rJF2, ti, tij, tijf ;
    Real   *pGx, *pGy, *pGz, *pHx, *pHy, *pHz, *values = NULL, *work = NULL ;
    Real    cI[3], cIJ[3], cIJF[3], *gT, *gx, *gy, *gz, *hx, *hy, *hz,
            xidg[MAXAMP1*MAXAMP1*MAXAMP1], yidg[MAXAMP1*MAXAMP1*MAXAMP1], zidg[MAXAMP1*MAXAMP1*MAXAMP1] ,
            xidh[MAXAMP1*MAXAMP1*MAXAMP1], yidh[MAXAMP1*MAXAMP1*MAXAMP1], zidh[MAXAMP1*MAXAMP1*MAXAMP1] ,
            xint[MAXAMP1*MAXAMP2*MAXAMP2], yint[MAXAMP1*MAXAMP2*MAXAMP2], zint[MAXAMP1*MAXAMP2*MAXAMP2] ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ         = ( iBasis == jBasis ) && ( rI == rJ ) ;
    for ( c = 0, rIJ2 = rIF2 = rJF2 = 0.0e+00 ; c < 3 ; c++ )
    {
        rIJ2 += pow ( rI[c] - rJ[c], 2 ) ;
        rIF2 += pow ( rI[c] - rF[c], 2 ) ;
        rJF2 += pow ( rJ[c] - rF[c], 2 ) ;
    }
    /* . Set pointers . */
    gT = &rWork[   0] ;
    gx = &rWork[  s3] ; gy = &rWork[2*s3] ; gz = &rWork[3*s3] ;
    hx = &rWork[4*s3] ; hy = &rWork[5*s3] ; hz = &rWork[6*s3] ;
    /* . Triple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        /* . Get information about the shell. */
        iamMax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF  ;
        /* . Set the Upper limit for the JSHELL loops. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            /* . Get information about the shell. */
            jamMax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF  ;
            /* . Set the diagonal block flag. */
            isDiagonal = iIsJ && ( iShell == jShell ) ;
            /* . Loop over fitting shells. */
            for ( fShell = 0 ; fShell < fBasis->nShells ; fShell++ )
            {
                /* . Get information about the shell. */
                fammax  = fBasis->shells[fShell].lHigh ;
                nCFuncF = fBasis->shells[fShell].nCBF  ;
                /* . Initialize the integral blocks. */
                for ( i = 0 ; i < ( nCFuncF * nCFuncI * nCFuncJ ) ; i++ ) gx[i] = gy[i] = gz[i] = hx[i] = hy[i] = hz[i] = 0.0e+00 ;
                /* . Set the array dimensions. */
                dim1 = fammax + 1 ;
                ddim2 = dim1 * ( jamMax + 1 ) ;
                dim2  = dim1 * ( jamMax + 2 ) ;
                /* . Triple loop over primitives. */
                for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
                {
                    /* . Get some information for the primitive. */
	            aI = iBasis->shells[iShell].primitives[ip].exponent ;
	            for ( i = 0 ; i < 3 ; i++ ) cI[i] = aI * rI[i] ;
                    /* . Inner loop over primitives. */
                    for ( jp = 0 ; jp < jBasis->shells[jShell].nPrimitives ; jp++ )
                    {
                        /* . Get some information for the primitive. */
	                aJ  = jBasis->shells[jShell].primitives[jp].exponent ;
                        aIJ = aI + aJ ;
	                eIJ = aI * aJ * rIJ2 ;
	                if ( ( eIJ / aIJ ) > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
	                for ( i = 0 ; i < 3 ; i++ ) cIJ[i] = ( cI[i] + aJ * rJ[i] ) ;
                        /* . Loop over fitting primitives. */
                        for ( fp = 0 ; fp < fBasis->shells[fShell].nPrimitives ; fp++ )
                        {
                            /* . Get some information for the primitive. */
	                    aF   = fBasis->shells[fShell].primitives[fp].exponent ;
                            aIJF = aIJ + aF ;
                            eIJF = ( eIJ + aI * aF * rIF2 + aJ * aF * rJF2 ) / aIJF ;
	                    if ( eIJF > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                            expfac = exp ( - eIJF ) ;
	                    for ( i = 0 ; i < 3 ; i++ ) cIJF[i] = ( cIJ[i] + aF * rF[i] ) / aIJF ;
                            /* . Calculate the overlap integrals. */
                            GaussianBasisSubsidiary_f1Og2 ( xint, yint, zint, aIJF, cIJF, rI, rJ, rF, iamMax+1, jamMax+1, fammax ) ;
                            GaussianBasisSubsidiary_f1Xg2r ( xint, yint, zint,
                                                              xidg, yidg, zidg,
                                                              xidh, yidh, zidh,
                                                              aI, aJ, iamMax, jamMax, fammax, dim1, dim2, dim1, ddim2 ) ;
                            /* . Assemble the integrals. */
                            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                            {
   	                        ix  = iBasis->shells[iShell].cbfPowX[i] * dim2  ;
	                        iy  = iBasis->shells[iShell].cbfPowY[i] * dim2  ;
	                        iz  = iBasis->shells[iShell].cbfPowZ[i] * dim2  ;
   	                        ixd = iBasis->shells[iShell].cbfPowX[i] * ddim2 ;
	                        iyd = iBasis->shells[iShell].cbfPowY[i] * ddim2 ;
	                        izd = iBasis->shells[iShell].cbfPowZ[i] * ddim2 ;
                                ti  = expfac * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                                for ( j = 0 ; j < nCFuncJ ; j++ )
                                {
	                            jix  = jBasis->shells[jShell].cbfPowX[j] * dim1 + ix  ;
	                            jiy  = jBasis->shells[jShell].cbfPowY[j] * dim1 + iy  ;
	                            jiz  = jBasis->shells[jShell].cbfPowZ[j] * dim1 + iz  ;
	                            jixd = jBasis->shells[jShell].cbfPowX[j] * dim1 + ixd ;
	                            jiyd = jBasis->shells[jShell].cbfPowY[j] * dim1 + iyd ;
	                            jizd = jBasis->shells[jShell].cbfPowZ[j] * dim1 + izd ;
                                    tij  = ti * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                                    for ( f = 0 ; f < nCFuncF ; f++, n++ )
                                    {
                                        fijx   = fBasis->shells[fShell].cbfPowX[f] + jix  ;
                                        fijy   = fBasis->shells[fShell].cbfPowY[f] + jiy  ;
                                        fijz   = fBasis->shells[fShell].cbfPowZ[f] + jiz  ;
                                        fijxd  = fBasis->shells[fShell].cbfPowX[f] + jixd ;
                                        fijyd  = fBasis->shells[fShell].cbfPowY[f] + jiyd ;
                                        fijzd  = fBasis->shells[fShell].cbfPowZ[f] + jizd ;
                                        tijf   = tij * fBasis->shells[fShell].primitives[fp].cCBF[f] ;
                                        gx[n] += tijf * xidg[fijxd] * yint[fijy ] * zint[fijz ] ;
                                        gy[n] += tijf * xint[fijx ] * yidg[fijyd] * zint[fijz ] ;
                                        gz[n] += tijf * xint[fijx ] * yint[fijy ] * zidg[fijzd] ;
                                        hx[n] += tijf * xidh[fijxd] * yint[fijy ] * zint[fijz ] ;
                                        hy[n] += tijf * xint[fijx ] * yidh[fijyd] * zint[fijz ] ;
                                        hz[n] += tijf * xint[fijx ] * yint[fijy ] * zidh[fijzd] ;
                                    }
                                }
                            }
                        } /* . fP. */
                    } /* . jP. */
                } /* . iP. */
                /* . Transform and save the integrals. */
                {
/*                    auto Boolean     skip ;*/
                    auto Integer      ii, jj, m, m3, m6, n ;
                    auto Cardinal16  *indices16 = block->indices16   ;
                    auto Real        *integrals = block->data, scale ;
                    auto RealArray2D *iC2S = iBasis->shells[iShell].c2s ,
                                     *jC2S = jBasis->shells[jShell].c2s ,
                                     *fC2S = fBasis->shells[fShell].c2s ;
                    work   = gT ;
                    values = gx ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pGx = values ;
                    values = gy ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pGy = values ;
                    values = gz ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pGz = values ;
                    values = hx ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pHx = values ;
                    values = hy ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pHy = values ;
                    values = hz ; GaussianBasisTransform3 ( nCFuncI, nCFuncJ, nCFuncF, iC2S, jC2S, fC2S, &values, &work ) ; pHz = values ;
                    if ( isDiagonal ) scale = 1.0e+00 ;
                    else              scale = 2.0e+00 ;
                    m = block->count ;
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
                    {
                        ii = iBasis->shells[iShell].nStart + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++ )
                        {
/*                            skip = isDiagonal && ( j > i ) ;*/
                            jj   = jBasis->shells[jShell].nStart + j ;
                            for ( f = 0 ; f < fBasis->shells[fShell].nBasis ; f++, m++, n++ )
                            {
/*
                                if ( ! skip )
                                {
*/
                                    m3              = 3 * m ;
                                    m6              = 6 * m ;
                                    indices16[m3  ] = ii ;
                                    indices16[m3+1] = jj ;
                                    indices16[m3+2] = fBasis->shells[fShell].nStart + f ;
                                    integrals[m6  ] = scale * pGx[n] ;
                                    integrals[m6+1] = scale * pGy[n] ;
                                    integrals[m6+2] = scale * pGz[n] ;
                                    integrals[m6+3] = scale * pHx[n] ;
                                    integrals[m6+4] = scale * pHy[n] ;
                                    integrals[m6+5] = scale * pHz[n] ;
                                    /*m++ ;*/
/*
                                }
*/
                            }
                        }
                    }
                    block->count = m ;
                }
            } /* . fShell. */
        } /* . jShell. */
    } /* . iShell. */
}
