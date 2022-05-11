/*==================================================================================================================================
! . Integrals - 2 bases, 1 electron, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "Boolean.h"
# include "GaussianBasisIntegrals_f1Xg1.h"
# include "GaussianBasisSubsidiary.h"
# include "GaussianBasisTransform.h"
# include "Integer.h"
# include "RysQuadrature.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Anti-Coulomb integrals.
! . integrals is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 6 * s2 and real 3 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1Ag1i ( const GaussianBasis *iBasis    ,
                                     const Real          *rI        ,
                                     const GaussianBasis *jBasis    ,
                                     const Real          *rJ        ,
                                     const Integer        s2        ,
                                           Integer       *iWork     ,
                                           Real          *rWork     ,
                                           RealArray2D   *integrals )
{
    Boolean       iIsJ ;
    Integer       gStrideI, hStrideI, i, iammax, ip, iShell, iGx, iGy, iGz, iHx, iHy, iHz,
                  j, jammax, jp, jShell, jUpper, m, n, nCFunc, nCFuncI, nCFuncJ, nRoots ;
    Real          aandb, ab, ai, aj, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dfi, dfIJ,
                  fac, fac2, f00, rho, rIJ2, tI, u2, xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real         *pG, *values = NULL, *work = NULL ;
    Integer      *IGx, *IGy, *IGz, *IHx, *IHy, *IHz ;
    Real         *Cij, *g, *gT,
                  Gx [MAXAMP3*MAXAMP3], Gy [MAXAMP3*MAXAMP3], Gz [MAXAMP3*MAXAMP3] ,
                  Hx [MAXAMP1*MAXAMP1], Hy [MAXAMP1*MAXAMP1], Hz [MAXAMP1*MAXAMP1] ;
    RysQuadrature roots ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( integrals, 0.0e+00 ) ;
    /* . Set pointers. */
    Cij = &rWork[   0] ; g   = &rWork[  s2] ; gT  = &rWork[2*s2] ;
    IGx = &iWork[   0] ; IGy = &iWork[  s2] ; IGz = &iWork[2*s2] ;
    IHx = &iWork[3*s2] ; IHy = &iWork[4*s2] ; IHz = &iWork[5*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF  ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].lHigh ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            /* . Roots. */
/* . +4! */
            nRoots   = ( iammax + jammax + 4 ) / 2 + 1 ;
            /* . Strides. */
            gStrideI = jammax + 3 ;
            hStrideI = jammax + 1 ;
            /* . Index arrays. */
            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
            {
   	        iGx = iBasis->shells[iShell].cbfPowX[i] * gStrideI ;
	        iGy = iBasis->shells[iShell].cbfPowY[i] * gStrideI ;
	        iGz = iBasis->shells[iShell].cbfPowZ[i] * gStrideI ;
   	        iHx = iBasis->shells[iShell].cbfPowX[i] * hStrideI ;
	        iHy = iBasis->shells[iShell].cbfPowY[i] * hStrideI ;
	        iHz = iBasis->shells[iShell].cbfPowZ[i] * hStrideI ;
                for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                {
	            IGx[n] = jBasis->shells[jShell].cbfPowX[j] + iGx ;
	            IGy[n] = jBasis->shells[jShell].cbfPowY[j] + iGy ;
	            IGz[n] = jBasis->shells[jShell].cbfPowZ[j] + iGz ;
	            IHx[n] = jBasis->shells[jShell].cbfPowX[j] + iHx ;
	            IHy[n] = jBasis->shells[jShell].cbfPowY[j] + iHy ;
	            IHz[n] = jBasis->shells[jShell].cbfPowZ[j] + iHz ;
                }
            }
            /* . Function initialization. */
            nCFunc = nCFuncI * nCFuncJ ;
            for ( n = 0 ; n < nCFunc ; n++ ) g[n] = 0.0e+00 ;
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai  = iBasis->shells[iShell].primitives[ip].exponent ;
                dfi = PI252 / ai ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nPrimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj = jBasis->shells[jShell].primitives[jp].exponent ;
                    /* . Calculate some factors. */
                    ab    = ai * aj ;
                    aandb = ai + aj ;
                    rho   = ab / aandb ;
                    dfIJ  = dfi / ( aj * sqrt ( aandb ) ) ;
                    /* . Calculate some displacements. */
                    c1x  =   ai * xIJ ;
                    c1y  =   ai * yIJ ;
                    c1z  =   ai * zIJ ;
                    c3x  = - aj * xIJ ;
                    c3y  = - aj * yIJ ;
                    c3z  = - aj * zIJ ;
                    /* . Calculate the rys polynomial roots. */
                    RysQuadrature_Roots ( &roots, nRoots, ( rho * rIJ2 ) ) ;
                    /* . Coefficient array. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
                        tI = dfIJ * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ ) Cij[n] = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                    }
                    /* . Loop over Rys roots. */
                    for ( m = 0 ; m < nRoots ; m++ )
                    {
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] ;
                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( ai + u2 ) * fac2 ;
                        b00   =        u2   * fac2 ;
                        b10   = ( aj + u2 ) * fac2 ;
                        xcp00 = u2 * c1x * fac ;
                        xc00  = u2 * c3x * fac ;
                        ycp00 = u2 * c1y * fac ;
                        yc00  = u2 * c3y * fac ;
                        zcp00 = u2 * c1z * fac ;
                        zc00  = u2 * c3z * fac ;
                        GaussianBasisSubsidiary_f1Cg1 ( iammax+2 ,
                                                        jammax+2 ,
                                                        b00      ,
                                                        b10      ,
                                                        bp01     ,
                                                        f00      ,
                                                        xc00     ,
                                                        xcp00    ,
                                                        yc00     ,
                                                        ycp00    ,
                                                        zc00     ,
                                                        zcp00    ,
                                                        gStrideI ,
                                                        Gx       ,
                                                        Gy       ,
                                                        Gz       ) ;
                        GaussianBasisSubsidiary_f1Ag1 ( iammax   ,
                                                        jammax   ,
                                                        gStrideI ,
                                                        Gx       ,
                                                        Gy       ,
                                                        Gz       ,
                                                        xIJ      ,
                                                        yIJ      ,
                                                        zIJ      ,
                                                        hStrideI ,
                                                        Hx       ,
                                                        Hy       ,
                                                        Hz       ) ;
                        /* . Assemble the integrals. */
                        for ( n = 0 ; n < nCFunc ; n++ )
                        {
                            g[n] += Cij[n] * ( Hx[IHx[n]] * Gy[IGy[n]] * Gz[IGz[n]] + 
                                               Gx[IGx[n]] * Hy[IHy[n]] * Gz[IGz[n]] + 
                                               Gx[IGx[n]] * Gy[IGy[n]] * Hz[IHz[n]] ) ;
                        }
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
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( integrals, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = -pG[n] ; /* . -r12 operator. */
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Anti-Coulomb derivatives.
! . sX, sY and sZ are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 9 * s2 and real 5 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1Ag1r1 ( const GaussianBasis *iBasis ,
                                      const Real          *rI     ,
                                      const GaussianBasis *jBasis ,
                                      const Real          *rJ     ,
                                      const Integer        s2     ,
                                            Integer       *iWork  ,
                                            Real          *rWork  ,
                                            RealArray2D   *sX     ,
                                            RealArray2D   *sY     ,
                                            RealArray2D   *sZ     )
{
    Integer       dStrideI, gStrideI, hStrideI,
                  i, iammax, ip, iShell,
                  iDx, iDy, iDz, iGx, iGy, iGz, iHx, iHy, iHz,
                  j, jammax, jp, jShell, m, n, nCFunc, nCFuncI, nCFuncJ, nRoots ;
    Real          aandb, ab, ai, aj, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dfi, dfIJ,
                  fac, fac2, f00, rho, rIJ2, tI, u2,
                  xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real         *pGx, *pGy, *pGz, *values = NULL, *work = NULL ;
    Integer      *IDx, *IDy, *IDz, *IGx, *IGy, *IGz, *IHx, *IHy, *IHz ;
    Real         *Cij, *gT, *gX, *gY, *gZ,
                  Gx [MAXAMP3*MAXAMP4], Gy [MAXAMP3*MAXAMP4], Gz [MAXAMP3*MAXAMP4] ,
                  GxD[MAXAMP1*MAXAMP1], GyD[MAXAMP1*MAXAMP1], GzD[MAXAMP1*MAXAMP1] ,
                  Hx [MAXAMP1*MAXAMP2], Hy [MAXAMP1*MAXAMP2], Hz [MAXAMP1*MAXAMP2] ,
                  HxD[MAXAMP1*MAXAMP1], HyD[MAXAMP1*MAXAMP1], HzD[MAXAMP1*MAXAMP1] ;
    RealArray2D  *iC2S, *jC2S ;
    RysQuadrature roots ;
    /* . Initialization. */
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( sX, 0.0e+00 ) ;
    RealArray2D_Set ( sY, 0.0e+00 ) ;
    RealArray2D_Set ( sZ, 0.0e+00 ) ;
    /* . Set pointers. */
    Cij = &rWork[   0] ; gT  = &rWork[  s2] ;
    gX  = &rWork[2*s2] ; gY  = &rWork[3*s2] ; gZ  = &rWork[4*s2] ;
    IDx = &iWork[   0] ; IDy = &iWork[  s2] ; IDz = &iWork[2*s2] ;
    IGx = &iWork[3*s2] ; IGy = &iWork[4*s2] ; IGz = &iWork[5*s2] ;
    IHx = &iWork[6*s2] ; IHy = &iWork[7*s2] ; IHz = &iWork[8*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        iC2S    = iBasis->shells[iShell].c2s ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jBasis->nShells ; jShell++ )
        {
            jammax   = jBasis->shells[jShell].lHigh ;
            jC2S     = jBasis->shells[jShell].c2s ;
            nCFuncJ  = jBasis->shells[jShell].nCBF     ;
            /* . Roots. */
            nRoots = ( iammax + jammax + 5 ) / 2 + 1 ; /* . +5! */
            /* . Strides. */
            dStrideI = jammax + 1 ;
            gStrideI = jammax + 3 ;
            hStrideI = jammax + 1 ;
            /* . Index arrays. */
            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
            {
   	        iDx = iBasis->shells[iShell].cbfPowX[i] * dStrideI ;
	        iDy = iBasis->shells[iShell].cbfPowY[i] * dStrideI ;
	        iDz = iBasis->shells[iShell].cbfPowZ[i] * dStrideI ;
   	        iGx = iBasis->shells[iShell].cbfPowX[i] * gStrideI ;
	        iGy = iBasis->shells[iShell].cbfPowY[i] * gStrideI ;
	        iGz = iBasis->shells[iShell].cbfPowZ[i] * gStrideI ;
   	        iHx = iBasis->shells[iShell].cbfPowX[i] * hStrideI ;
	        iHy = iBasis->shells[iShell].cbfPowY[i] * hStrideI ;
	        iHz = iBasis->shells[iShell].cbfPowZ[i] * hStrideI ;
                for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                {
	            IDx[n] = jBasis->shells[jShell].cbfPowX[j] + iDx ;
	            IDy[n] = jBasis->shells[jShell].cbfPowY[j] + iDy ;
	            IDz[n] = jBasis->shells[jShell].cbfPowZ[j] + iDz ;
	            IGx[n] = jBasis->shells[jShell].cbfPowX[j] + iGx ;
	            IGy[n] = jBasis->shells[jShell].cbfPowY[j] + iGy ;
	            IGz[n] = jBasis->shells[jShell].cbfPowZ[j] + iGz ;
	            IHx[n] = jBasis->shells[jShell].cbfPowX[j] + iHx ;
	            IHy[n] = jBasis->shells[jShell].cbfPowY[j] + iHy ;
	            IHz[n] = jBasis->shells[jShell].cbfPowZ[j] + iHz ;
                }
            }
            /* . Initialize the integral blocks. */
            nCFunc = nCFuncI * nCFuncJ ;
            for ( n = 0 ; n < nCFunc ; n++ ) gX[n] = gY[n] = gZ[n] = 0.0e+00 ;
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai  = iBasis->shells[iShell].primitives[ip].exponent ;
                dfi = PI252 / ai ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nPrimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj = jBasis->shells[jShell].primitives[jp].exponent ;
                    /* . Calculate some factors. */
                    ab    = ai * aj ;
                    aandb = ai + aj ;
                    rho   = ab / aandb ;
                    dfIJ  = dfi / ( aj * sqrt ( aandb ) ) ;
                    /* . Calculate some displacements. */
                    c1x  =   ai * xIJ ;
                    c1y  =   ai * yIJ ;
                    c1z  =   ai * zIJ ;
                    c3x  = - aj * xIJ ;
                    c3y  = - aj * yIJ ;
                    c3z  = - aj * zIJ ;
                    /* . Calculate the rys polynomial roots. */
                    RysQuadrature_Roots ( &roots, nRoots, ( rho * rIJ2 ) ) ;
                    /* . Coefficient array. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
                        tI = dfIJ * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ ) Cij[n] = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                    }
                    /* . Loop over the roots and construct the subsidiary integrals. */
                    for ( m = 0 ; m < nRoots ; m++ )
                    {
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] ;
                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( ai + u2 ) * fac2 ;
                        b00   =        u2   * fac2 ;
                        b10   = ( aj + u2 ) * fac2 ;
                        xcp00 = u2 * c1x * fac ;
                        xc00  = u2 * c3x * fac ;
                        ycp00 = u2 * c1y * fac ;
                        yc00  = u2 * c3y * fac ;
                        zcp00 = u2 * c1z * fac ;
                        zc00  = u2 * c3z * fac ;
                        GaussianBasisSubsidiary_f1Cg1  ( iammax+3 , /* . +1. */
                                                         jammax+2 ,
                                                         b00      ,
                                                         b10      ,
                                                         bp01     ,
                                                         f00      ,
                                                         xc00     ,
                                                         xcp00    ,
                                                         yc00     ,
                                                         ycp00    ,
                                                         zc00     ,
                                                         zcp00    ,
                                                         gStrideI ,
                                                         Gx       ,
                                                         Gy       ,
                                                         Gz       ) ;
                        GaussianBasisSubsidiary_f1Ag1  ( iammax+1 , /* . +1. */
                                                         jammax   ,
                                                         gStrideI ,
                                                         Gx       ,
                                                         Gy       ,
                                                         Gz       ,
                                                         xIJ      ,
                                                         yIJ      ,
                                                         zIJ      ,
                                                         hStrideI ,
                                                         Hx       ,
                                                         Hy       ,
                                                         Hz       ) ;
                        GaussianBasisSubsidiary_f1Xg1r ( Gx       ,
                                                         Gy       ,
                                                         Gz       ,
                                                         ai       ,
                                                         iammax   ,
                                                         jammax   ,
                                                         gStrideI ,
                                                         dStrideI ,
                                                         GxD      ,
                                                         GyD      ,
                                                         GzD      ) ;
                        GaussianBasisSubsidiary_f1Xg1r ( Hx       ,
                                                         Hy       ,
                                                         Hz       ,
                                                         ai       ,
                                                         iammax   ,
                                                         jammax   ,
                                                         hStrideI ,
                                                         dStrideI ,
                                                         HxD      ,
                                                         HyD      ,
                                                         HzD      ) ;
                        /* . Assemble the integrals. */
                        for ( n = 0 ; n < nCFunc ; n++ )
                        {
                            gX[n] += Cij[n] * ( HxD[IDx[n]] * Gy [IGy[n]] * Gz [IGz[n]] + 
                                                GxD[IDx[n]] * Hy [IHy[n]] * Gz [IGz[n]] + 
                                                GxD[IDx[n]] * Gy [IGy[n]] * Hz [IHz[n]] ) ;
                            gY[n] += Cij[n] * ( Hx [IHx[n]] * GyD[IDy[n]] * Gz [IGz[n]] + 
                                                Gx [IGx[n]] * HyD[IDy[n]] * Gz [IGz[n]] + 
                                                Gx [IGx[n]] * GyD[IDy[n]] * Hz [IHz[n]] ) ;
                            gZ[n] += Cij[n] * ( Hx [IHx[n]] * Gy [IGy[n]] * GzD[IDz[n]] + 
                                                Gx [IGx[n]] * Hy [IHy[n]] * GzD[IDz[n]] + 
                                                Gx [IGx[n]] * Gy [IGy[n]] * HzD[IDz[n]] ) ;
                        }
                    } /* . nRoots. */
                } /* . jP. */
            } /* . iP. */
            /* . Transform the integrals. */
            work   = gT ;
            values = gX ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pGx = values ;
            values = gY ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pGy = values ;
            values = gZ ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pGz = values ;
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( sX, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = -pGx[n] ; /* . -r12 operator. */
                    Array2D_Item ( sY, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = -pGy[n] ;
                    Array2D_Item ( sZ, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = -pGz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Coulomb integrals.
! . integrals is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s2 and real 3 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1Cg1i ( const GaussianBasis *iBasis    ,
                                     const Real          *rI        ,
                                     const GaussianBasis *jBasis    ,
                                     const Real          *rJ        ,
                                     const Integer        s2        ,
                                           Integer       *iWork     ,
                                           Real          *rWork     ,
                                           RealArray2D   *integrals )
{
    Boolean       iIsJ ;
    Integer       i, iammax, ip, iShell, ix, iy, iz, j, jammax, jdim,
                  jp, jShell, jUpper, m, n, nCFunc, nCFuncI, nCFuncJ, nRoots ;
    Real          aandb, ab, ai, aj, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dfi, dfIJ,
                  fac, fac2, f00, rho, rIJ2, tI, u2, xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real         *pG, *values = NULL, *work = NULL ;
    Integer      *Ix, *Iy, *Iz ;
    Real         *Cij, *g, *gT ,
                  xint[MAXAMP1*MAXAMP1], yint[MAXAMP1*MAXAMP1], zint[MAXAMP1*MAXAMP1] ;
    RysQuadrature roots ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( integrals, 0.0e+00 ) ;
    /* . Set pointers. */
    Cij = &rWork[   0] ; g  = &rWork[s2] ; gT = &rWork[2*s2] ;
    Ix  = &iWork[   0] ; Iy = &iWork[s2] ; Iz = &iWork[2*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF  ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].lHigh ;
            jdim    = jammax + 1 ;
            nCFuncJ = jBasis->shells[jShell].nCBF ;
            /* . Roots. */
            nRoots = ( iammax + jammax ) / 2 + 1 ;
            /* . Index arrays. */
            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
            {
   	        ix = iBasis->shells[iShell].cbfPowX[i] * jdim ;
	        iy = iBasis->shells[iShell].cbfPowY[i] * jdim ;
	        iz = iBasis->shells[iShell].cbfPowZ[i] * jdim ;
                for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                {
	            Ix[n] = jBasis->shells[jShell].cbfPowX[j] + ix ;
	            Iy[n] = jBasis->shells[jShell].cbfPowY[j] + iy ;
	            Iz[n] = jBasis->shells[jShell].cbfPowZ[j] + iz ;
                }
            }
            /* . Function initialization. */
            nCFunc = nCFuncI * nCFuncJ ;
            for ( n = 0 ; n < nCFunc ; n++ ) g[n] = 0.0e+00 ;
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai  = iBasis->shells[iShell].primitives[ip].exponent ;
                dfi = PI252 / ai ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nPrimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj = jBasis->shells[jShell].primitives[jp].exponent ;
                    /* . Calculate some factors. */
                    ab    = ai * aj ;
                    aandb = ai + aj ;
                    rho   = ab / aandb ;
                    dfIJ  = dfi / ( aj * sqrt ( aandb ) ) ;
                    /* . Calculate some displacements. */
                    c1x  = ai * ( rI[0] - rJ[0] ) ;
                    c1y  = ai * ( rI[1] - rJ[1] ) ;
                    c1z  = ai * ( rI[2] - rJ[2] ) ;
                    c3x  = aj * ( rJ[0] - rI[0] ) ;
                    c3y  = aj * ( rJ[1] - rI[1] ) ;
                    c3z  = aj * ( rJ[2] - rI[2] ) ;
                    /* . Calculate the rys polynomial roots. */
                    RysQuadrature_Roots ( &roots, nRoots, ( rho * rIJ2 ) ) ;
                    /* . Coefficient array. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
                        tI = dfIJ * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ ) Cij[n] = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                    }
                    /* . Loop over Rys roots. */
                    for ( m = 0 ; m < nRoots ; m++ )
                    {
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] ;
                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( ai + u2 ) * fac2 ;
                        b00   =        u2   * fac2 ;
                        b10   = ( aj + u2 ) * fac2 ;
                        xcp00 = u2 * c1x * fac ;
                        xc00  = u2 * c3x * fac ;
                        ycp00 = u2 * c1y * fac ;
                        yc00  = u2 * c3y * fac ;
                        zcp00 = u2 * c1z * fac ;
                        zc00  = u2 * c3z * fac ;
                        GaussianBasisSubsidiary_f1Cg1 ( iammax ,
                                                        jammax ,
                                                        b00    ,
                                                        b10    ,
                                                        bp01   ,
                                                        f00    ,
                                                        xc00   ,
                                                        xcp00  ,
                                                        yc00   ,
                                                        ycp00  ,
                                                        zc00   ,
                                                        zcp00  ,
                                                        jdim   ,
                                                        xint   ,
                                                        yint   ,
                                                        zint   ) ;
                        /* . Assemble the integrals. */
                        for ( n = 0 ; n < nCFunc ; n++ ) g[n] += ( Cij[n] * xint[Ix[n]] * yint[Iy[n]] * zint[Iz[n]] ) ;
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
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( integrals, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pG[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Coulomb derivatives.
! . sX, sY and sZ are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s2 and real 5 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1Cg1r1 ( const GaussianBasis *iBasis ,
                                      const Real          *rI     ,
                                      const GaussianBasis *jBasis ,
                                      const Real          *rJ     ,
                                      const Integer        s2     ,
                                            Integer       *iWork  ,
                                            Real          *rWork  ,
                                            RealArray2D   *sX     ,
                                            RealArray2D   *sY     ,
                                            RealArray2D   *sZ     )
{
    Integer       i, iammax, ip, iShell, ix, iy, iz,
                  j, jammax, jdim, jp, jShell,
                  m, n, nCFunc, nCFuncI, nCFuncJ, nRoots ;
    Real          aandb, ab, ai, aj, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dfi, dfIJ,
                  fac, fac2, f00, rho, rIJ2, tI, u2,
                  xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real         *pGx, *pGy, *pGz, *values = NULL, *work = NULL ;
    Integer      *Ix, *Iy, *Iz ;
    Real         *Cij, *gT, *gX, *gY, *gZ,
                  xind[MAXAMP1*MAXAMP1], yind[MAXAMP1*MAXAMP1], zind[MAXAMP1*MAXAMP1] ,
                  xint[MAXAMP1*MAXAMP2], yint[MAXAMP1*MAXAMP2], zint[MAXAMP1*MAXAMP2] ;
    RealArray2D  *iC2S, *jC2S ;
    RysQuadrature roots ;
    /* . Initialization. */
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( sX, 0.0e+00 ) ;
    RealArray2D_Set ( sY, 0.0e+00 ) ;
    RealArray2D_Set ( sZ, 0.0e+00 ) ;
    /* . Set pointers. */
    Cij = &rWork[   0] ; gT = &rWork[  s2] ;
    gX  = &rWork[2*s2] ; gY = &rWork[3*s2] ; gZ = &rWork[4*s2] ;
    Ix  = &iWork[   0] ; Iy = &iWork[  s2] ; Iz = &iWork[2*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        iC2S    = iBasis->shells[iShell].c2s ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jBasis->nShells ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].lHigh ;
            jdim    = jammax + 1 ;
            jC2S    = jBasis->shells[jShell].c2s ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            /* . Roots. */
            nRoots = ( iammax + jammax + 1 ) / 2 + 1 ;
            /* . Index arrays. */
            for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
            {
   	        ix = iBasis->shells[iShell].cbfPowX[i] * jdim ;
	        iy = iBasis->shells[iShell].cbfPowY[i] * jdim ;
	        iz = iBasis->shells[iShell].cbfPowZ[i] * jdim ;
                for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                {
	            Ix[n] = jBasis->shells[jShell].cbfPowX[j] + ix ;
	            Iy[n] = jBasis->shells[jShell].cbfPowY[j] + iy ;
	            Iz[n] = jBasis->shells[jShell].cbfPowZ[j] + iz ;
                }
            }
            /* . Initialize the integral blocks. */
            nCFunc = nCFuncI * nCFuncJ ;
            for ( n = 0 ; n < nCFunc ; n++ ) gX[n] = gY[n] = gZ[n] = 0.0e+00 ;
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai  = iBasis->shells[iShell].primitives[ip].exponent ;
                dfi = PI252 / ai ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nPrimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj = jBasis->shells[jShell].primitives[jp].exponent ;
                    /* . Calculate some factors. */
                    ab    = ai * aj ;
                    aandb = ai + aj ;
                    rho   = ab / aandb ;
                    dfIJ  = dfi / ( aj * sqrt ( aandb ) ) ;
                    /* . Calculate some displacements. */
                    c1x  = ai * ( rI[0] - rJ[0] ) ;
                    c1y  = ai * ( rI[1] - rJ[1] ) ;
                    c1z  = ai * ( rI[2] - rJ[2] ) ;
                    c3x  = aj * ( rJ[0] - rI[0] ) ;
                    c3y  = aj * ( rJ[1] - rI[1] ) ;
                    c3z  = aj * ( rJ[2] - rI[2] ) ;
                    /* . Calculate the rys polynomial roots. */
                    RysQuadrature_Roots ( &roots, nRoots, ( rho * rIJ2 ) ) ;
                    /* . Coefficient array. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
                        tI = dfIJ * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ ) Cij[n] = tI * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
                    }
                    /* . Loop over the roots and construct the subsidiary integrals. */
                    for ( m = 0 ; m < nRoots ; m++ )
                    {
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] ;
                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( ai + u2 ) * fac2 ;
                        b00   =        u2   * fac2 ;
                        b10   = ( aj + u2 ) * fac2 ;
                        xcp00 = u2 * c1x * fac ;
                        xc00  = u2 * c3x * fac ;
                        ycp00 = u2 * c1y * fac ;
                        yc00  = u2 * c3y * fac ;
                        zcp00 = u2 * c1z * fac ;
                        zc00  = u2 * c3z * fac ;
                        GaussianBasisSubsidiary_f1Cg1  ( iammax+1 ,
                                                         jammax   ,
                                                         b00      ,
                                                         b10      ,
                                                         bp01     ,
                                                         f00      ,
                                                         xc00     ,
                                                         xcp00    ,
                                                         yc00     ,
                                                         ycp00    ,
                                                         zc00     ,
                                                         zcp00    ,
                                                         jdim     ,
                                                         xint     ,
                                                         yint     ,
                                                         zint     ) ;
                        GaussianBasisSubsidiary_f1Xg1r ( xint     ,
                                                         yint     ,
                                                         zint     ,
                                                         ai       ,
                                                         iammax   ,
                                                         jammax   ,
                                                         jdim     ,
                                                         jdim     ,
                                                         xind     ,
                                                         yind     ,
                                                         zind     ) ;
                        /* . Assemble the integrals. */
                        for ( n = 0 ; n < nCFunc ; n++ )
                        {
                            gX[n] += ( Cij[n] * xind[Ix[n]] * yint[Iy[n]] * zint[Iz[n]] ) ;
                            gY[n] += ( Cij[n] * xint[Ix[n]] * yind[Iy[n]] * zint[Iz[n]] ) ;
                            gZ[n] += ( Cij[n] * xint[Ix[n]] * yint[Iy[n]] * zind[Iz[n]] ) ;
                        }
                    } /* . nRoots. */
                } /* . jP. */
            } /* . iP. */
            /* . Transform the integrals. */
            work   = gT ;
            values = gX ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pGx = values ;
            values = gY ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pGy = values ;
            values = gZ ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pGz = values ;
             /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( sX, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pGx[n] ;
                    Array2D_Item ( sY, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pGy[n] ;
                    Array2D_Item ( sZ, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pGz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole integrals.
! . dipoleX, dipoleY and dipoleZ are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 4 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1Df1i ( const GaussianBasis *iBasis  ,
                                     const Real          *rI      ,
                                     const GaussianBasis *jBasis  ,
                                     const Real          *rJ      ,
                                     const Real          *center  ,
                                     const Integer        s2      ,
                                           Real          *rWork   ,
                                           RealArray2D   *dipoleX ,
                                           RealArray2D   *dipoleY ,
                                           RealArray2D   *dipoleZ )
{
    Boolean      iIsJ ;
    Integer      i, iammax, ip, iShell,  ix, iy, iz,
                 j, jammax, jdim, jp,  jShell,  jUpper, jxix, jyiy, jziz, n, nCFuncI, nCFuncJ ;
    Real         aa, aainv, ai, aj, arri, expfac, fac, rIJ2, ti, tIJ, xIJ, yIJ, zIJ ;
    Real        *pSx, *pSy, *pSz, *values = NULL, *work = NULL ;
    Real         ar[3], arI[3], *gT, *sx, *sy, *sz ,
                 xo[MAXAMP1*MAXAMP1], yo[MAXAMP1*MAXAMP1], zo[MAXAMP1*MAXAMP1] ,
                 xd[MAXAMP1*MAXAMP1], yd[MAXAMP1*MAXAMP1], zd[MAXAMP1*MAXAMP1] ;
    RealArray2D *iC2S, *jC2S ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( dipoleX, 0.0e+00 ) ;
    RealArray2D_Set ( dipoleY, 0.0e+00 ) ;
    RealArray2D_Set ( dipoleZ, 0.0e+00 ) ;
    /* . Set pointers. */
    gT = &rWork[0] ; sx = &rWork[s2] ; sy = &rWork[2*s2] ; sz = &rWork[3*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        iC2S    = iBasis->shells[iShell].c2s ;
        nCFuncI = iBasis->shells[iShell].nCBF ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].lHigh ;
            jdim    = jammax + 1 ;
            jC2S    = jBasis->shells[jShell].c2s ;
            nCFuncJ = jBasis->shells[jShell].nCBF ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( nCFuncI * nCFuncJ ) ; i++ )
            {
                sx[i] = 0.0e+00 ;
                sy[i] = 0.0e+00 ;
                sz[i] = 0.0e+00 ;
            }
            /* . Outer loop over primitives. */
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
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                    /* . Calculate the subsidiary integrals. */
                    GaussianBasisSubsidiary_f1Og1 ( xo, yo, zo, aa, ar, rI, rJ,         iammax, jammax ) ;
                    GaussianBasisSubsidiary_f1Dg1   ( xd, yd, zd, aa, ar, rI, rJ, center, iammax, jammax ) ;
                    /* . Add in the contributions to the full integrals. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
   	                ix = iBasis->shells[iShell].cbfPowX[i] * jdim ;
	                iy = iBasis->shells[iShell].cbfPowY[i] * jdim ;
	                iz = iBasis->shells[iShell].cbfPowZ[i] * jdim ;
                        ti = expfac * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                        {
	                    jxix = jBasis->shells[jShell].cbfPowX[j] + ix ;
	                    jyiy = jBasis->shells[jShell].cbfPowY[j] + iy ;
	                    jziz = jBasis->shells[jShell].cbfPowZ[j] + iz ;
                            tIJ  = ti * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
    		            sx[n] += tIJ * xd[jxix] * yo[jyiy] * zo[jziz] ;
    		            sy[n] += tIJ * xo[jxix] * yd[jyiy] * zo[jziz] ;
    		            sz[n] += tIJ * xo[jxix] * yo[jyiy] * zd[jziz] ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Transform the integrals. */
            work   = gT ;
            values = sx ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSx = values ;
            values = sy ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSy = values ;
            values = sz ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSz = values ;
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( dipoleX, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSx[n] ;
                    Array2D_Item ( dipoleY, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSy[n] ;
                    Array2D_Item ( dipoleZ, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Kinetic energy and overlap integrals.
! . Kinetic and overlap are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 3 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1KOg1i ( const GaussianBasis *iBasis  ,
                                      const Real          *rI      ,
                                      const GaussianBasis *jBasis  ,
                                      const Real          *rJ      ,
                                      const Integer        s2      ,
                                            Real          *rWork   ,
                                            RealArray2D   *overlap ,
                                            RealArray2D   *kinetic )
{
    Boolean      iIsJ ;
    Integer      i, iammax, ip, iShell,  ixo, ixt, iyo, iyt, izo, izt,
                 j, jammax, jdimo, jdimt, jp,  jShell,  jUpper,
                 jxixo, jyiyo, jzizo, jxixt, jyiyt, jzizt, n, nCFuncI, nCFuncJ ;
    Real         aa, aainv, ai, aj, arri, expfac, fac, rIJ2, ti, tIJ, xIJ, yIJ, zIJ ;
    Real        *pS, *pT, *values = NULL, *work = NULL ;
    Real         ar[3], arI[3], *gT, *s, *t,
                 xo[MAXAMP1*MAXAMP3], yo[MAXAMP1*MAXAMP3], zo[MAXAMP1*MAXAMP3] ,
                 xt[MAXAMP1*MAXAMP1], yt[MAXAMP1*MAXAMP1], zt[MAXAMP1*MAXAMP1] ;
    RealArray2D *iC2S, *jC2S ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( kinetic, 0.0e+00 ) ;
    RealArray2D_Set ( overlap, 0.0e+00 ) ;
    /* . Set pointers. */
    gT = &rWork[0] ; s = &rWork[s2] ; t = &rWork[2*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        iC2S    = iBasis->shells[iShell].c2s ;
        nCFuncI = iBasis->shells[iShell].nCBF ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].lHigh ;
            jdimo   = jammax + 3 ;
            jdimt   = jammax + 1 ;
            jC2S    = jBasis->shells[jShell].c2s ;
            nCFuncJ = jBasis->shells[jShell].nCBF ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( nCFuncI * nCFuncJ ) ; i++ )
            {
                s[i] = 0.0e+00 ;
                t[i] = 0.0e+00 ;
            }
            /* . Outer loop over primitives. */
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
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;

                    /* . Calculate the subsidiary integrals. */
                    GaussianBasisSubsidiary_f1Og1 ( xo, yo, zo, aa, ar, rI, rJ, iammax, jammax + 2 ) ;
                    GaussianBasisSubsidiary_f1Kg1  ( xo, yo, zo, xt, yt, zt, aj, iammax, jammax, jdimo, jdimt ) ;

                    /* . Add in the contributions to the full integrals. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
   	                ixo = iBasis->shells[iShell].cbfPowX[i] * jdimo ;
	                iyo = iBasis->shells[iShell].cbfPowY[i] * jdimo ;
	                izo = iBasis->shells[iShell].cbfPowZ[i] * jdimo ;
   	                ixt = iBasis->shells[iShell].cbfPowX[i] * jdimt ;
	                iyt = iBasis->shells[iShell].cbfPowY[i] * jdimt ;
	                izt = iBasis->shells[iShell].cbfPowZ[i] * jdimt ;
                        ti = expfac * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                        {
	                    jxixo = jBasis->shells[jShell].cbfPowX[j] + ixo ;
	                    jyiyo = jBasis->shells[jShell].cbfPowY[j] + iyo ;
	                    jzizo = jBasis->shells[jShell].cbfPowZ[j] + izo ;
	                    jxixt = jBasis->shells[jShell].cbfPowX[j] + ixt ;
	                    jyiyt = jBasis->shells[jShell].cbfPowY[j] + iyt ;
	                    jzizt = jBasis->shells[jShell].cbfPowZ[j] + izt ;
                            tIJ  = ti * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
    		            s[n] += tIJ * xo[jxixo] * yo[jyiyo] * zo[jzizo] ;
                            t[n] += tIJ * ( xt[jxixt] * yo[jyiyo] * zo[jzizo] +
                                            xo[jxixo] * yt[jyiyt] * zo[jzizo] +
                                            xo[jxixo] * yo[jyiyo] * zt[jzizt] ) ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Transform the integrals. */
            work   = gT ;
            values = s  ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pS = values ;
            values = t  ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pT = values ;
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( overlap, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pS[n] ;
                    Array2D_Item ( kinetic, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pT[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Kinetic energy and overlap derivatives.
! . kineticX, kineticY, kineticZ, overlapX, overlapY and overlapZ are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 7 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1KOg1r1 ( const GaussianBasis *iBasis   ,
                                       const Real          *rI       ,
                                       const GaussianBasis *jBasis   ,
                                       const Real          *rJ       ,
                                       const Integer        s2       ,
                                             Real          *rWork    ,
                                             RealArray2D   *overlapX ,
                                             RealArray2D   *overlapY ,
                                             RealArray2D   *overlapZ ,
                                             RealArray2D   *kineticX ,
                                             RealArray2D   *kineticY ,
                                             RealArray2D   *kineticZ )
{
     Integer     i, iammax, ip, iShell,  ixo, ixt, iyo, iyt, izo, izt,       
                 j, jammax, jdimo, jdimt, jp,  jShell,  jxixo, jyiyo, jzizo, 
                 jxixt, jyiyt, jzizt, n, nCFuncI, nCFuncJ ;                           
     Real        aa, aainv, ai, aj, arri, expfac, fac, rIJ2, ti, tIJ, xIJ, yIJ, zIJ ; 
     Real       *pSx, *pSy, *pSz, *pTx, *pTy, *pTz, *values = NULL, *work = NULL ;  
     Real        ar[3], arI[3], *gT, *sx, *sy, *sz, *tx, *ty, *tz,
                 xo[MAXAMP2*MAXAMP3], yo[MAXAMP2*MAXAMP3], zo[MAXAMP2*MAXAMP3],       
                 xt[MAXAMP1*MAXAMP2], yt[MAXAMP1*MAXAMP2], zt[MAXAMP1*MAXAMP2],       
                 xod[MAXAMP1*MAXAMP3], yod[MAXAMP1*MAXAMP3], zod[MAXAMP1*MAXAMP3],    
                 xtd[MAXAMP1*MAXAMP1], ytd[MAXAMP1*MAXAMP1], ztd[MAXAMP1*MAXAMP1] ;   
    RealArray2D *iC2S, *jC2S ;
    /* . Initialization. */
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( kineticX, 0.0e+00 ) ;
    RealArray2D_Set ( kineticY, 0.0e+00 ) ;
    RealArray2D_Set ( kineticZ, 0.0e+00 ) ;
    RealArray2D_Set ( overlapX, 0.0e+00 ) ;
    RealArray2D_Set ( overlapY, 0.0e+00 ) ;
    RealArray2D_Set ( overlapZ, 0.0e+00 ) ;
    /* . Set pointers. */
    gT = &rWork[   0] ;
    sx = &rWork[  s2] ; sy = &rWork[2*s2] ; sz = &rWork[3*s2] ;
    tx = &rWork[4*s2] ; ty = &rWork[5*s2] ; tz = &rWork[6*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        iC2S    = iBasis->shells[iShell].c2s ;
        nCFuncI = iBasis->shells[iShell].nCBF ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jBasis->nShells ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].lHigh ;
            jdimo   = jammax + 3 ;
            jdimt   = jammax + 1 ;
            jC2S    = jBasis->shells[jShell].c2s ;
            nCFuncJ = jBasis->shells[jShell].nCBF ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( nCFuncI * nCFuncJ ) ; i++ )
            {
                sx[i] = sy[i] = sz[i] = 0.0e+00 ;
                tx[i] = ty[i] = tz[i] = 0.0e+00 ;
            }
            /* . Outer loop over primitives. */
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
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                    /* . Calculate the subsidiary integrals. */
                    GaussianBasisSubsidiary_f1Og1 ( xo, yo, zo, aa, ar, rI, rJ, iammax + 1, jammax + 2 ) ;
                    GaussianBasisSubsidiary_f1Kg1  ( xo, yo, zo, xt, yt, zt, aj, iammax + 1, jammax, jdimo, jdimt ) ;
                    GaussianBasisSubsidiary_f1Xg1r ( xo, yo, zo, ai, iammax, jammax, jdimo, jdimo, xod, yod, zod ) ;
                    GaussianBasisSubsidiary_f1Xg1r ( xt, yt, zt, ai, iammax, jammax, jdimt, jdimt, xtd, ytd, ztd ) ;
                    /* . Add in the contributions to the full integrals. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
   	                ixo = iBasis->shells[iShell].cbfPowX[i] * jdimo ;
	                iyo = iBasis->shells[iShell].cbfPowY[i] * jdimo ;
	                izo = iBasis->shells[iShell].cbfPowZ[i] * jdimo ;
   	                ixt = iBasis->shells[iShell].cbfPowX[i] * jdimt ;
	                iyt = iBasis->shells[iShell].cbfPowY[i] * jdimt ;
	                izt = iBasis->shells[iShell].cbfPowZ[i] * jdimt ;
                        ti  = expfac * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                        {
	                    jxixo = jBasis->shells[jShell].cbfPowX[j] + ixo ;
	                    jyiyo = jBasis->shells[jShell].cbfPowY[j] + iyo ;
	                    jzizo = jBasis->shells[jShell].cbfPowZ[j] + izo ;
	                    jxixt = jBasis->shells[jShell].cbfPowX[j] + ixt ;
	                    jyiyt = jBasis->shells[jShell].cbfPowY[j] + iyt ;
	                    jzizt = jBasis->shells[jShell].cbfPowZ[j] + izt ;
                            tIJ   = ti * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
    		            sx[n] += tIJ * xod[jxixo] * yo[jyiyo] * zo[jzizo] ;
    		            sy[n] += tIJ * xo[jxixo] * yod[jyiyo] * zo[jzizo] ;
    		            sz[n] += tIJ * xo[jxixo] * yo[jyiyo] * zod[jzizo] ;
                            tx[n] += tIJ * ( xtd[jxixt] * yo[jyiyo] * zo[jzizo] +
                                             xod[jxixo] * yt[jyiyt] * zo[jzizo] +
                                             xod[jxixo] * yo[jyiyo] * zt[jzizt] ) ;
                            ty[n] += tIJ * ( xt[jxixt] * yod[jyiyo] * zo[jzizo] +
                                             xo[jxixo] * ytd[jyiyt] * zo[jzizo] +
                                             xo[jxixo] * yod[jyiyo] * zt[jzizt] ) ;
                            tz[n] += tIJ * ( xt[jxixt] * yo[jyiyo] * zod[jzizo] +
                                             xo[jxixo] * yt[jyiyt] * zod[jzizo] +
                                             xo[jxixo] * yo[jyiyo] * ztd[jzizt] ) ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Transform the integrals. */
            work   = gT ;
            values = sx ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSx = values ;
            values = sy ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSy = values ;
            values = sz ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSz = values ;
            values = tx ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pTx = values ;
            values = ty ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pTy = values ;
            values = tz ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pTz = values ;
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( overlapX, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSx[n] ;
                    Array2D_Item ( overlapY, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSy[n] ;
                    Array2D_Item ( overlapZ, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSz[n] ;
                    Array2D_Item ( kineticX, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pTx[n] ;
                    Array2D_Item ( kineticY, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pTy[n] ;
                    Array2D_Item ( kineticZ, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pTz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Overlap integrals.
! . integrals is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 2 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1Og1i ( const GaussianBasis *iBasis    ,
                                     const Real          *rI        ,
                                     const GaussianBasis *jBasis    ,
                                     const Real          *rJ        ,
                                     const Integer        s2        ,
                                           Real          *rWork     ,
                                           RealArray2D   *integrals )
{
    Boolean  iIsJ ;
    Integer  i, iammax, ip, iShell, ix, iy, iz,
             j, jammax, jdim, jp, jShell, jUpper, jxix, jyiy, jziz,
             n, nCFuncI, nCFuncJ ;
    Real     aa, aainv, ai, aj, arri, expfac, fac, rIJ2, xIJ, yIJ, zIJ ;
    Real    *pS, *values = NULL, *work = NULL ;
    Real     ar[3], arI[3], *gT, *s,
             xo[MAXAMP1*MAXAMP1], yo[MAXAMP1*MAXAMP1], zo[MAXAMP1*MAXAMP1] ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( integrals, 0.0e+00 ) ;
    /* . Set pointers. */
    gT = &rWork[0] ; s = &rWork[s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        nCFuncI = iBasis->shells[iShell].nCBF     ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].lHigh ;
            jdim    = jammax + 1 ;
            nCFuncJ = jBasis->shells[jShell].nCBF     ;
            /* . Initialize the integral block. */
            for ( i = 0 ; i < ( nCFuncI * nCFuncJ ) ; i++ ) s[i] = 0.0e+00 ;
            /* . Outer loop over primitives. */
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
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                    /* . Calculate the overlap integrals. */
                    GaussianBasisSubsidiary_f1Og1 ( xo, yo, zo, aa, ar, rI, rJ, iammax, jammax ) ;
                    /* . Add in the contributions to the full integrals. */
                    for ( i = n = 0 ; i < nCFuncI ; i++ )
                    {
   	                ix = iBasis->shells[iShell].cbfPowX[i] * jdim ;
	                iy = iBasis->shells[iShell].cbfPowY[i] * jdim ;
	                iz = iBasis->shells[iShell].cbfPowZ[i] * jdim ;
                        for ( j = 0 ; j < nCFuncJ ; j++ )
                        {
	                    jxix = jBasis->shells[jShell].cbfPowX[j] + ix ;
	                    jyiy = jBasis->shells[jShell].cbfPowY[j] + iy ;
	                    jziz = jBasis->shells[jShell].cbfPowZ[j] + iz ;
    		            s[n] += expfac * iBasis->shells[iShell].primitives[ip].cCBF[i] *
                                             jBasis->shells[jShell].primitives[jp].cCBF[j] * xo[jxix] * yo[jyiy] * zo[jziz] ;
                            n++ ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Transform the integrals. */
            values = s  ;
            work   = gT ;
            GaussianBasisTransform2 ( nCFuncI                    ,
                                      nCFuncJ                    ,
                                      iBasis->shells[iShell].c2s ,
                                      jBasis->shells[jShell].c2s ,
                                      &values                    ,
                                      &work                      ) ;
            pS = values ;
            /* . Put the integrals in the proper place. */
            for ( i = n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( integrals, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pS[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Overlap derivatives.
! . overlapX, overlapY and overlapZ are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 4 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1Og1r1 ( const GaussianBasis *iBasis   ,
                                      const Real          *rI       ,
                                      const GaussianBasis *jBasis   ,
                                      const Real          *rJ       ,
                                      const Integer        s2       ,
                                            Real          *rWork    ,
                                            RealArray2D   *overlapX ,
                                            RealArray2D   *overlapY ,
                                            RealArray2D   *overlapZ )
{
    Integer      i, iammax, ip, iShell, ix, iy, iz,
                 j, jammax, jdim, jp, jShell, jxix, jyiy, jziz, n, nCFuncI, nCFuncJ ;
    Real         aa, aainv, ai, aj, arri, denfac, expfac, fac, rIJ2, xIJ, yIJ, zIJ ;
    Real        *pSx, *pSy, *pSz, *values = NULL, *work = NULL ;
    Real         ar[3], arI[3], *gT, *sx, *sy, *sz ,
                 xd[MAXAMP1*MAXAMP1], yd[MAXAMP1*MAXAMP1], zd[MAXAMP1*MAXAMP1] ,
                 xo[MAXAMP1*MAXAMP2], yo[MAXAMP1*MAXAMP2], zo[MAXAMP1*MAXAMP2] ;
    RealArray2D *iC2S, *jC2S ;
    /* . Initialization. */
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( overlapX, 0.0e+00 ) ;
    RealArray2D_Set ( overlapY, 0.0e+00 ) ;
    RealArray2D_Set ( overlapZ, 0.0e+00 ) ;
    /* . Set pointers. */
    gT = &rWork[ 0] ;
    sx = &rWork[s2] ; sy = &rWork[2*s2] ; sz = &rWork[3*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        iC2S    = iBasis->shells[iShell].c2s ;
        nCFuncI = iBasis->shells[iShell].nCBF ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jBasis->nShells ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].lHigh ;
            jdim    = jammax + 1 ;
            jC2S    = jBasis->shells[jShell].c2s ;
            nCFuncJ = jBasis->shells[jShell].nCBF ;
            /* . Initialize the integral block. */
            for ( i = 0 ; i < ( nCFuncI * nCFuncJ ) ; i++ ) { sx[i] = 0.0e+00 ; sy[i] = 0.0e+00 ; sz[i] = 0.0e+00 ; }
            /* . Outer loop over primitives. */
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
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                    /* . Calculate the overlap integrals and derivatives. */
                    GaussianBasisSubsidiary_f1Og1    ( xo, yo, zo, aa, ar, rI, rJ, iammax+1, jammax ) ;
                    GaussianBasisSubsidiary_f1Xg1r ( xo, yo, zo, ai, iammax, jammax, jdim, jdim, xd, yd, zd ) ;
                    /* . Add in the contributions to the full integrals. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
   	                ix = iBasis->shells[iShell].cbfPowX[i] * jdim ;
	                iy = iBasis->shells[iShell].cbfPowY[i] * jdim ;
	                iz = iBasis->shells[iShell].cbfPowZ[i] * jdim ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                        {
	                    jxix = jBasis->shells[jShell].cbfPowX[j] + ix ;
	                    jyiy = jBasis->shells[jShell].cbfPowY[j] + iy ;
	                    jziz = jBasis->shells[jShell].cbfPowZ[j] + iz ;
                            denfac = expfac * iBasis->shells[iShell].primitives[ip].cCBF[i] * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
    		            sx[n] += denfac * xd[jxix] * yo[jyiy] * zo[jziz] ;
    		            sy[n] += denfac * xo[jxix] * yd[jyiy] * zo[jziz] ;
    		            sz[n] += denfac * xo[jxix] * yo[jyiy] * zd[jziz] ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Transform the integrals. */
            work   = gT ;
            values = sx ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSx = values ;
            values = sy ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSy = values ;
            values = sz ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSz = values ;
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( overlapX, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSx[n] ;
                    Array2D_Item ( overlapY, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSy[n] ;
                    Array2D_Item ( overlapZ, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Quadrupole integrals.
! . The quadrupole arrays are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 7 * s2 where s2 = ( maximum shell size )^2. */
void GaussianBasisIntegrals_f1Qf1i ( const GaussianBasis *iBasis  ,
                                     const Real          *rI      ,
                                     const GaussianBasis *jBasis  ,
                                     const Real          *rJ      ,
                                     const Real          *center  ,
                                     const Integer        s2      ,
                                           Real          *rWork   ,
                                           RealArray2D   *qXX     ,
                                           RealArray2D   *qYY     ,
                                           RealArray2D   *qZZ     ,
                                           RealArray2D   *qXY     ,
                                           RealArray2D   *qXZ     ,
                                           RealArray2D   *qYZ     )
{
    Boolean      iIsJ ;
    Integer      i, iammax, ip, iShell,  ix, iy, iz,
                 j, jammax, jdim, jp,  jShell,  jUpper, jxix, jyiy, jziz, n, nCFuncI, nCFuncJ ;
    Real         aa, aainv, ai, aj, arri, expfac, fac, rIJ2, ti, tIJ, xIJ, yIJ, zIJ ;
    Real        *pSxx, *pSxy, *pSxz, *pSyy, *pSyz, *pSzz, *values = NULL, *work = NULL ;
    Real         ar[3], arI[3], *gT, *sXX, *sXY, *sXZ, *sYY, *sYZ, *sZZ,
                 xo[MAXAMP1*MAXAMP1], yo[MAXAMP1*MAXAMP1], zo[MAXAMP1*MAXAMP1] ,
                 xd[MAXAMP1*MAXAMP1], yd[MAXAMP1*MAXAMP1], zd[MAXAMP1*MAXAMP1] ,
                 xq[MAXAMP1*MAXAMP1], yq[MAXAMP1*MAXAMP1], zq[MAXAMP1*MAXAMP1] ;
    RealArray2D *iC2S, *jC2S ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( qXX, 0.0e+00 ) ;
    RealArray2D_Set ( qYY, 0.0e+00 ) ;
    RealArray2D_Set ( qZZ, 0.0e+00 ) ;
    RealArray2D_Set ( qXY, 0.0e+00 ) ;
    RealArray2D_Set ( qXZ, 0.0e+00 ) ;
    RealArray2D_Set ( qYZ, 0.0e+00 ) ;
    /* . Set pointers. */
    gT  = &rWork[   0] ;
    sXX = &rWork[  s2] ; sXY = &rWork[2*s2] ; sXZ = &rWork[3*s2] ;
    sYY = &rWork[4*s2] ; sYZ = &rWork[5*s2] ; sZZ = &rWork[6*s2] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        iC2S    = iBasis->shells[iShell].c2s ;
        nCFuncI = iBasis->shells[iShell].nCBF ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nShells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].lHigh ;
            jdim    = jammax + 1 ;
            jC2S    = jBasis->shells[jShell].c2s ;
            nCFuncJ = jBasis->shells[jShell].nCBF ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( nCFuncI * nCFuncJ ) ; i++ )
            {
                sXX[i] = 0.0e+00 ;
                sYY[i] = 0.0e+00 ;
                sZZ[i] = 0.0e+00 ;
                sXY[i] = 0.0e+00 ;
                sXZ[i] = 0.0e+00 ;
                sYZ[i] = 0.0e+00 ;
            }
            /* . Outer loop over primitives. */
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
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                    /* . Calculate the subsidiary integrals. */
                    GaussianBasisSubsidiary_f1Og1   ( xo, yo, zo, aa, ar, rI, rJ,         iammax, jammax ) ;
                    GaussianBasisSubsidiary_f1Dg1     ( xd, yd, zd, aa, ar, rI, rJ, center, iammax, jammax ) ;
                    GaussianBasisSubsidiary_f1Qg1 ( xq, yq, zq, aa, ar, rI, rJ, center, iammax, jammax ) ;
                    /* . Add in the contributions to the full integrals. */
                    for ( i = 0, n = 0 ; i < nCFuncI ; i++ )
                    {
   	                ix = iBasis->shells[iShell].cbfPowX[i] * jdim ;
	                iy = iBasis->shells[iShell].cbfPowY[i] * jdim ;
	                iz = iBasis->shells[iShell].cbfPowZ[i] * jdim ;
                        ti = expfac * iBasis->shells[iShell].primitives[ip].cCBF[i] ;
                        for ( j = 0 ; j < nCFuncJ ; j++, n++ )
                        {
	                    jxix = jBasis->shells[jShell].cbfPowX[j] + ix ;
	                    jyiy = jBasis->shells[jShell].cbfPowY[j] + iy ;
	                    jziz = jBasis->shells[jShell].cbfPowZ[j] + iz ;
                            tIJ  = ti * jBasis->shells[jShell].primitives[jp].cCBF[j] ;
    		            sXX[n] += tIJ * xq[jxix] * yo[jyiy] * zo[jziz] ;
    		            sYY[n] += tIJ * xo[jxix] * yq[jyiy] * zo[jziz] ;
    		            sZZ[n] += tIJ * xo[jxix] * yo[jyiy] * zq[jziz] ;
    		            sXY[n] += tIJ * xd[jxix] * yd[jyiy] * zo[jziz] ;
    		            sXZ[n] += tIJ * xd[jxix] * yo[jyiy] * zd[jziz] ;
    		            sYZ[n] += tIJ * xo[jxix] * yd[jyiy] * zd[jziz] ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Transform the integrals. */
            work   = gT  ;
            values = sXX ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSxx = values ;
            values = sYY ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSyy = values ;
            values = sZZ ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSzz = values ;
            values = sXY ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSxy = values ;
            values = sXZ ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSxz = values ;
            values = sYZ ; GaussianBasisTransform2 ( nCFuncI, nCFuncJ, iC2S, jC2S, &values, &work ) ; pSyz = values ;
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nBasis ; j++, n++ )
                {
                    Array2D_Item ( qXX, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSxx[n] ;
                    Array2D_Item ( qYY, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSyy[n] ;
                    Array2D_Item ( qZZ, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSzz[n] ;
                    Array2D_Item ( qXY, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSxy[n] ;
                    Array2D_Item ( qXZ, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSxz[n] ;
                    Array2D_Item ( qYZ, i+iBasis->shells[iShell].nStart, j+jBasis->shells[jShell].nStart ) = pSyz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

