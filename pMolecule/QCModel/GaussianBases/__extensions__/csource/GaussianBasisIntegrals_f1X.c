/*==================================================================================================================================
! . Integrals - 1 basis, 1 electron.
!=================================================================================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "GaussianBasisIntegrals_f1X.h"
# include "GaussianBasisSubsidiary.h"
# include "GaussianBasisTransform.h"
# include "Integer.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the Cartesian multipole integrals for a single basis function.
!
!   All integrals are over the range (-Inf,+Inf) and have the form:
!
!   Sum_p c_p int ( x - Ax )^Nx exp ( - a_p * ( x - Ax )^2 ) dx * ( " y ) * ( " z )
!
!   with extra multipole operators of:
!
!   overlap    : 1
!   dipole     : ( x - Xc ), etc.
!   quadrupole : ( x - Xc )^2 or ( x - Xc ) * ( y - Yc ), etc.
!
!   Only zero or even 1-D polynomial integrals are non-zero. These have the form:
!
!   I_N(a) = (N-1)!!/(2 a)^(N/2) (Pi/a)^(1/2)
!
!   or if N = 2S:
!
!   (2S-1)!!/(2 a)^S (Pi/a)^(1/2)
!
!   (2S-1)!! = (2S-1) * (2S-3) * (2S-5) ... 5 * 3 * 1 (as 2S is even).
!
!   For dipole and quadrupole integrals have:
!
!   ( x - Cx )  :       D_N(a) = I_(N+1) (a) + D * I_N (a)
!   ( x - Cx )^2:       Q_N(a) = D_(N+1) (a) + D * D_N (a)
!                              = I_(N+2) (a) + 2 * D * I_(N+1) (a) + D^2 * I_N (a)
!
!   where D = ( Ax - Cx ).
!
! . All multipole arrays are overwritten by these functions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 4 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Di ( const GaussianBasis *iBasis ,
                                   const Real          *rI     ,
                                   const Real          *rC     ,
                                   const Integer        s1     ,
                                         Real          *rWork  ,
                                         RealArray1D   *dX     ,
                                         RealArray1D   *dY     ,
                                         RealArray1D   *dZ     )
{
    Integer      i, iammax, iP, iShell, iX, iY, iZ, ncfuncI ; 
    Real         aI, cI, xIC, yIC, zIC ;                               
    Real        *pSx, *pSy, *pSz, *values = NULL, *work = NULL ;     
    Real        *gT, *sX, *sY, *sZ,  
                 oO[MAXAMP2], xD[MAXAMP1], yD[MAXAMP1], zD[MAXAMP1] ;  
    RealArray2D *iC2S ;
    /* . Initialization. */
    xIC = rI[0]-rC[0] ;
    yIC = rI[1]-rC[1] ;
    zIC = rI[2]-rC[2] ;
    /* . Set pointers. */
    gT = &rWork[ 0] ;
    sX = &rWork[s1] ; sY = &rWork[2*s1] ; sZ = &rWork[3*s1] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        iC2S    = iBasis->shells[iShell].c2s ;
        ncfuncI = iBasis->shells[iShell].nCBF     ;
        /* . Initialize the integral blocks. */
        for ( i = 0 ; i < ncfuncI ; i++ )
        {
            sX[i] = 0.0e+00 ;
            sY[i] = 0.0e+00 ;
            sZ[i] = 0.0e+00 ;
        }
        /* . Outer loop over primitives. */
        for ( iP = 0 ; iP < iBasis->shells[iShell].nPrimitives ; iP++ )
        {
            /* . Get some information for the primitive. */
	    aI = iBasis->shells[iShell].primitives[iP].exponent ;
            /* . Calculate the subsidiary integrals. */
            _1Overlap    (     oO, aI , iammax+1 ) ;
            _1Derivative ( xD, oO, xIC, iammax   ) ;
            _1Derivative ( yD, oO, yIC, iammax   ) ;
            _1Derivative ( zD, oO, zIC, iammax   ) ;
            /* . Add in the contributions to the full integrals. */
            for ( i = 0 ; i < ncfuncI ; i++ )
            {
   	        iX     = iBasis->shells[iShell].cbfPowX[i] ;
	        iY     = iBasis->shells[iShell].cbfPowY[i] ;
	        iZ     = iBasis->shells[iShell].cbfPowZ[i] ;
                cI     = iBasis->shells[iShell].primitives[iP].cCBF[i] ;
    		sX[i] += cI * xD[iX] * oO[iY] * oO[iZ] ;
    		sY[i] += cI * oO[iX] * yD[iY] * oO[iZ] ;
    		sZ[i] += cI * oO[iX] * oO[iY] * zD[iZ] ;
/* printf ( "1D> %5d %5d %5d %5d %5d %5d | %10.3f %10.3f | %10.3f %10.3f %10.3f | %10.3f %10.3f %10.3f | %10.3f %10.3f %10.3f\n",
             iShell, iP, i, iX, iY, iZ, aI, cI, xIC, yIC, zIC, xD[iX], yD[iY], zD[iZ], oO[iX], oO[iY], oO[iZ] ) ; fflush ( stdout ) ; */
            }
        } /* . iP. */
        /* . Transform the integrals. */
        work   = gT ;
        values = sX ; GaussianBasisTransform1 ( iC2S, &values, &work ) ; pSx = values ;
        values = sY ; GaussianBasisTransform1 ( iC2S, &values, &work ) ; pSy = values ;
        values = sZ ; GaussianBasisTransform1 ( iC2S, &values, &work ) ; pSz = values ;
        /* . Put the integrals in the proper place. */
        for ( i = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
        {
            Array1D_Item ( dX, i+iBasis->shells[iShell].nStart ) = pSx[i] ;
            Array1D_Item ( dY, i+iBasis->shells[iShell].nStart ) = pSy[i] ;
            Array1D_Item ( dZ, i+iBasis->shells[iShell].nStart ) = pSz[i] ;
        }
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Overlap integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 2 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Oi ( const GaussianBasis *iBasis  ,
                                   const Integer        s1      ,
                                         Real          *rWork   ,
                                         RealArray1D   *overlap )
{
    Integer  i, iammax, iP, iShell, iX, iY, iZ, ncfuncI ;
    Real     aI, cI ;
    Real    *pO3, *values = NULL, *work = NULL ;
    Real    *gT, *o3, o1[MAXAMP1] ;
    /* . Set pointers. */
    gT = &rWork[0] ; o3 = &rWork[s1] ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        ncfuncI = iBasis->shells[iShell].nCBF  ;
        /* . Initialize the integral blocks. */
        for ( i = 0 ; i < ncfuncI ; i++ ) o3[i] = 0.0e+00 ;
        /* . Outer loop over primitives. */
        for ( iP = 0 ; iP < iBasis->shells[iShell].nPrimitives ; iP++ )
        {
            /* . Get some information for the primitive. */
	    aI = iBasis->shells[iShell].primitives[iP].exponent ;
            /* . Calculate the subsidiary integrals. */
            _1Overlap ( o1, aI, iammax ) ;
            /* . Add in the contributions to the full integrals. */
            for ( i = 0 ; i < ncfuncI ; i++ )
            {
   	        iX     = iBasis->shells[iShell].cbfPowX[i] ;
	        iY     = iBasis->shells[iShell].cbfPowY[i] ;
	        iZ     = iBasis->shells[iShell].cbfPowZ[i] ;
                cI     = iBasis->shells[iShell].primitives[iP].cCBF[i] ;
    		o3[i] += cI * o1[iX] * o1[iY] * o1[iZ] ;
            }
        } /* . iP. */
        /* . Transform the integrals. */
        values = o3 ;
        work   = gT ;
        GaussianBasisTransform1 ( iBasis->shells[iShell].c2s, &values, &work ) ;
        pO3 = values ;
        /* . Put the integrals in the proper place. */
        for ( i = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
        {
            Array1D_Item ( overlap, i+iBasis->shells[iShell].nStart ) = pO3[i] ;
        }
    } /* . iShell. */
}

/* . Old code. */
/*
void GaussianBasisIntegrals_f1Oi ( const GaussianBasis *iBasis, RealArray1D *overlap )
{
    if ( ( iBasis != NULL ) && ( overlap != NULL ) )
    {
        auto Integer  i, ip, iShell, ix, iy, iz, nCFuncI, t ;
        auto Real     ci, ei, si ;
        for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
        {
            nCFuncI = iBasis->shells[iShell].nCBF     ;
            for ( i = 0 ; i < nCFuncI ; i++ )
            {
   	        ix = iBasis->shells[iShell].cbfPowX[i] ;
	        iy = iBasis->shells[iShell].cbfPowY[i] ;
	        iz = iBasis->shells[iShell].cbfPowZ[i] ;
                if ( IsEven ( ix ) && IsEven ( iy ) && IsEven ( iz ) )
                {
                    si = 0.0e+00 ;
                    for ( ip = 0 ; ip < iBasis->shells[iShell].nPrimitives ; ip++ )
                    {
                        ci = iBasis->shells[iShell].primitives[ip].cCBF[i]  ;
                        ei = iBasis->shells[iShell].primitives[ip].exponent ;
                        for ( t = 1 ; t <= ( ix / 2 ) ; t++ ) ci *= ( Real ) ( 2 * t - 1 ) / ( 2.0e+00 * ei ) ;
                        for ( t = 1 ; t <= ( iy / 2 ) ; t++ ) ci *= ( Real ) ( 2 * t - 1 ) / ( 2.0e+00 * ei ) ;
                        for ( t = 1 ; t <= ( iz / 2 ) ; t++ ) ci *= ( Real ) ( 2 * t - 1 ) / ( 2.0e+00 * ei ) ;
                        si += ci / ( ei * sqrt ( ei ) ) ;
                    }
                    Array1D_Item ( overlap, iBasis->shells[iShell].nStart+i ) = PI32 * si ;
                }
            }
        }
    }
}
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Quadrupole integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 7 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Qi ( const GaussianBasis *iBasis ,
                                   const Real          *rI     ,
                                   const Real          *rC     ,
                                   const Integer        s1     ,
                                         Real          *rWork  ,
                                         RealArray1D   *qXX    ,
                                         RealArray1D   *qYY    ,
                                         RealArray1D   *qZZ    ,
                                         RealArray1D   *qXY    ,
                                         RealArray1D   *qXZ    ,
                                         RealArray1D   *qYZ    )
{
    Integer      i, iammax, iP, iShell, iX, iY, iZ, ncfuncI ;                     
    Real         aI, cI, xIC, yIC, zIC ;                                                   
    Real        *pSxx, *pSyy, *pSzz, *pSxy, *pSxz, *pSyz, *values = NULL, *work = NULL ; 
    Real        *gT, *sXX, *sXY, *sXZ, *sYY, *sYZ, *sZZ,                   
                 oO[MAXAMP3], xD [MAXAMP2], yD [MAXAMP2], zD [MAXAMP2] ,                   
                              xQ [MAXAMP1], yQ [MAXAMP1], zQ [MAXAMP1] ;                   
    RealArray2D *iC2S ;
    /* . Initialization. */
    xIC = rI[0]-rC[0] ;
    yIC = rI[1]-rC[1] ;
    zIC = rI[2]-rC[2] ;
    /* . Set pointers. */
    gT  = &rWork[   0] ;
    sXX = &rWork[  s1] ; sXY = &rWork[2*s1] ; sXZ = &rWork[3*s1] ;
    sYY = &rWork[4*s1] ; sYZ = &rWork[5*s1] ; sZZ = &rWork[6*s1] ;  
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].lHigh ;
        iC2S    = iBasis->shells[iShell].c2s ;
        ncfuncI = iBasis->shells[iShell].nCBF     ;
        /* . Initialize the integral blocks. */
        for ( i = 0 ; i < ncfuncI ; i++ )
        {
            sXX[i] = 0.0e+00 ;
            sYY[i] = 0.0e+00 ;
            sZZ[i] = 0.0e+00 ;
            sXY[i] = 0.0e+00 ;
            sXZ[i] = 0.0e+00 ;
            sYZ[i] = 0.0e+00 ;
        }
        /* . Outer loop over primitives. */
        for ( iP = 0 ; iP < iBasis->shells[iShell].nPrimitives ; iP++ )
        {
            /* . Get some information for the primitive. */
	    aI = iBasis->shells[iShell].primitives[iP].exponent ;
            /* . Calculate the subsidiary integrals. */
            _1Overlap    (     oO, aI , iammax+2 ) ;
            _1Derivative ( xD, oO, xIC, iammax+1 ) ;
            _1Derivative ( yD, oO, yIC, iammax+1 ) ;
            _1Derivative ( zD, oO, zIC, iammax+1 ) ;
            _1Derivative ( xQ, xD, xIC, iammax   ) ;
            _1Derivative ( yQ, yD, yIC, iammax   ) ;
            _1Derivative ( zQ, zD, zIC, iammax   ) ;
            /* . Add in the contributions to the full integrals. */
            for ( i = 0 ; i < ncfuncI ; i++ )
            {
   	        iX      = iBasis->shells[iShell].cbfPowX[i] ;
	        iY      = iBasis->shells[iShell].cbfPowY[i] ;
	        iZ      = iBasis->shells[iShell].cbfPowZ[i] ;
                cI      = iBasis->shells[iShell].primitives[iP].cCBF[i] ;
    		sXX[i] += cI * xQ[iX] * oO[iY] * oO[iZ] ;
    		sYY[i] += cI * oO[iX] * yQ[iY] * oO[iZ] ;
    		sZZ[i] += cI * oO[iX] * oO[iY] * zQ[iZ] ;
    		sXY[i] += cI * xD[iX] * yD[iY] * oO[iZ] ;
    		sXZ[i] += cI * xD[iX] * oO[iY] * zD[iZ] ;
    		sYZ[i] += cI * oO[iX] * yD[iY] * zD[iZ] ;
            }
        } /* . iP. */
        /* . Transform the integrals. */
        work   = gT  ;
        values = sXX ; GaussianBasisTransform1 ( iC2S, &values, &work ) ; pSxx = values ;
        values = sYY ; GaussianBasisTransform1 ( iC2S, &values, &work ) ; pSyy = values ;
        values = sZZ ; GaussianBasisTransform1 ( iC2S, &values, &work ) ; pSzz = values ;
        values = sXY ; GaussianBasisTransform1 ( iC2S, &values, &work ) ; pSxy = values ;
        values = sXZ ; GaussianBasisTransform1 ( iC2S, &values, &work ) ; pSxz = values ;
        values = sYZ ; GaussianBasisTransform1 ( iC2S, &values, &work ) ; pSyz = values ;
        /* . Put the integrals in the proper place. */
        for ( i = 0 ; i < iBasis->shells[iShell].nBasis ; i++ )
        {
            Array1D_Item ( qXX, i+iBasis->shells[iShell].nStart ) = pSxx[i] ;
            Array1D_Item ( qYY, i+iBasis->shells[iShell].nStart ) = pSyy[i] ;
            Array1D_Item ( qZZ, i+iBasis->shells[iShell].nStart ) = pSzz[i] ;
            Array1D_Item ( qXY, i+iBasis->shells[iShell].nStart ) = pSxy[i] ;
            Array1D_Item ( qXZ, i+iBasis->shells[iShell].nStart ) = pSxz[i] ;
            Array1D_Item ( qYZ, i+iBasis->shells[iShell].nStart ) = pSyz[i] ;
        }
    } /* . iShell. */
}
