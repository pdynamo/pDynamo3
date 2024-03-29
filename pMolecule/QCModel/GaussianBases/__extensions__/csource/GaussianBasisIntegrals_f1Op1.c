/*==================================================================================================================================
! . Integrals - 1 basis, 0 electrons, 1 nucleus/point.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "GaussianBasisIntegrals_f1Op1.h"
# include "GaussianBasisTransform.h"

/*
! . Notes:
!
!   - The order of derivatives is independent of the Cartesian basis function order.
!   - The shape of the output matrices is N * G.
!
!   - The output arrays should be appropriately initialized before entry to all of these functions.
!
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at the given points.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 3 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Op1i ( const GaussianBasis *iBasis ,
                                     const Real          *rI     ,
                                     const Coordinates3  *rG     ,
                                     const Integer        s1     ,
                                           Real          *rWork  ,
                                           RealArray2D   *f      )
{
    if ( ( iBasis != NULL ) && 
         ( rI     != NULL ) && 
         ( rG     != NULL ) && 
         ( f      != NULL ) )
    {
        Integer  g, i, ip, iShell, iStart, ix, iy, iz, nCFuncI ;
        Real     dx, dy, dz, eip, r2 ;
        Real    *pG0, *values = NULL, *work = NULL ;     
        Real    *e0, *g0, *gT, x0[MAXAMP1], y0[MAXAMP1], z0[MAXAMP1] ;
        /* . Set pointers. */
        e0 = &rWork[0] ; g0 = &rWork[s1] ; gT = &rWork[2*s1] ;
        /* . Loop over points. */
        for ( g = 0 ; g < Coordinates3_Rows ( rG ) ; g++ )
        {
            /* . Get the distance squared between atom centers. */
            dx = Coordinates3_Item ( rG, g, 0 ) - rI[0] ;
            dy = Coordinates3_Item ( rG, g, 1 ) - rI[1] ;
            dz = Coordinates3_Item ( rG, g, 2 ) - rI[2] ;
            r2 = dx * dx + dy * dy + dz * dz ;
            /* . Check for negligible contributions. */
            /* if ( r2 >= iBasis->maximumRange ) return ; */
            /* . Form the angular functions. */
            x0[0] = 1.0e+00 ; y0[0] = 1.0e+00 ; z0[0] = 1.0e+00 ;
            for ( i = 1 ; i <= iBasis->lHigh ; i++ )
            {
                x0[i] = dx * x0[i-1] ;
                y0[i] = dy * y0[i-1] ;
                z0[i] = dz * z0[i-1] ;
            }
            /* . Loop over shells. */
            for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
            {
                /* . Get information about the shell. */
                nCFuncI = iBasis->shells[iShell].nCBF     ;
                /* . Check for negligible contributions. */
                /* if ( r2 >= iBasis->shells[iShell].maximumRange ) return ; */
                /* . Form the exponential factors. */
                for ( i = 0 ; i < nCFuncI ; i++ ) { e0[i] = 0.0e+00 ; }
                for ( ip = 0; ip < iBasis->shells[iShell].nPrimitives ; ip ++ )
                {
                    eip = exp ( - iBasis->shells[iShell].primitives[ip].exponent * r2 ) ;
                    for ( i = 0 ; i < nCFuncI ; i++ ) { e0[i] += iBasis->shells[iShell].primitives[ip].cCBF[i] * eip ; }
                }
                /* . Form the Cartesian function values. */
                for ( i = 0 ; i < nCFuncI ; i++ )
                {
                    ix = iBasis->shells[iShell].cbfPowX[i] ;
                    iy = iBasis->shells[iShell].cbfPowY[i] ;
                    iz = iBasis->shells[iShell].cbfPowZ[i] ;
                    g0[i] = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                }
                /* . Transform the integrals. */
                values = g0 ;
                work   = gT ;
                GaussianBasisTransform1 ( iBasis->shells[iShell].c2s, &values, &work ) ;
                pG0 = values ;
                /* . Put the values in the proper place. */
                iStart = iBasis->shells[iShell].nStart ;
                for ( i = 0 ; i < iBasis->shells[iShell].nBasis ; i++ ) { Array2D_Item ( f, i+iStart, g ) = pG0[i] ; }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at the given points and their first derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 10 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Op1ir1 ( const GaussianBasis *iBasis ,
                                       const Real          *rI     ,
                                       const Coordinates3  *rG     ,
                                       const Integer        s1     ,
                                             Real          *rWork  ,
                                             RealArray2D   *f      ,
                                             RealArray2D   *fX     ,
                                             RealArray2D   *fY     ,
                                             RealArray2D   *fZ     )
{
    if ( ( iBasis != NULL ) && 
         ( rI     != NULL ) && 
         ( rG     != NULL ) && 
         ( f      != NULL ) && 
         ( fX     != NULL ) && 
         ( fY     != NULL ) && 
         ( fZ     != NULL ) )
    {
        Integer  g, i, ip, iShell, iStart, ix, iy, iz, nCFuncI, nFuncI ;
        Real     dx, dy, dz, eip, r2 ;
        Real    *pG01, *values = NULL, *work = NULL ;     
        Real    *e0, *e1, *g01, *gT,
                 x0[MAXAMP2], y0[MAXAMP2], z0[MAXAMP2] ,
                 x1[MAXAMP2], y1[MAXAMP2], z1[MAXAMP2] ;
        /* . Set pointers. */
        e0 = &rWork[0] ; e1 = &rWork[s1] ; g01 = &rWork[2*s1] ; gT = &rWork[6*s1] ; /* . Note 4 each for g01 and gT! */
        /* . Loop over points. */
        for ( g = 0 ; g < Coordinates3_Rows ( rG ) ; g++ )
        {
            /* . Get the distance squared between atom centers. */
            dx = Coordinates3_Item ( rG, g, 0 ) - rI[0] ;
            dy = Coordinates3_Item ( rG, g, 1 ) - rI[1] ;
            dz = Coordinates3_Item ( rG, g, 2 ) - rI[2] ;
            r2 = dx * dx + dy * dy + dz * dz ;
            /* . Check for negligible contributions. */
            /* if ( r2 >= iBasis->maximumRange ) return ; */
            /* . Form the angular functions. */
            x0[0] = 1.0e+00 ; y0[0] = 1.0e+00 ; z0[0] = 1.0e+00 ;
            x1[0] = 0.0e+00 ; y1[0] = 0.0e+00 ; z1[0] = 0.0e+00 ;
            for ( i = 1 ; i <= iBasis->lHigh + 1; i++ )
            {
                x0[i] = dx * x0[i-1] ;
                y0[i] = dy * y0[i-1] ;
                z0[i] = dz * z0[i-1] ;
                x1[i] = ( Real ) ( i ) * x0[i-1] ;
                y1[i] = ( Real ) ( i ) * y0[i-1] ;
                z1[i] = ( Real ) ( i ) * z0[i-1] ;
            }
            /* . Outer loop over shells. */
            for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
            {
                /* . Get information about the shell. */
                nCFuncI = iBasis->shells[iShell].nCBF     ;
                /* . Check for negligible contributions. */
                /* if ( r2 >= iBasis->shells[iShell].maximumRange ) return ; */
                /* . Form the exponential factors. */
                for ( i = 0 ; i < nCFuncI ; i++ ) { e0[i] = 0.0e+00 ; e1[i] = 0.0e+00 ; }
                for ( ip = 0; ip < iBasis->shells[iShell].nPrimitives ; ip ++ )
                {
                    eip = exp ( - iBasis->shells[iShell].primitives[ip].exponent * r2 ) ;
                    for ( i = 0 ; i < nCFuncI ; i++ )
                    {
                        e0[i] += iBasis->shells[iShell].primitives[ip].cCBF[i]  * eip ;
                        e1[i] += iBasis->shells[iShell].primitives[ip].cCBF[i]  *
                                 iBasis->shells[iShell].primitives[ip].exponent * eip ;
                    }
                }
                for ( i = 0 ; i < nCFuncI ; i++ ) { e1[i] *= -2.0e+00 ; }
                /* . Form the Cartesian function values. */
                for ( i = 0 ; i < nCFuncI ; i++ )
                {
                    ix = iBasis->shells[iShell].cbfPowX[i] ;
                    iy = iBasis->shells[iShell].cbfPowY[i] ;
                    iz = iBasis->shells[iShell].cbfPowZ[i] ;
                    g01[i]           = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                    g01[i+  nCFuncI] = ( x1[ix] * e0[i] + x0[ix+1] * e1[i] ) * y0[iy] * z0[iz] ;
                    g01[i+2*nCFuncI] = ( y1[iy] * e0[i] + y0[iy+1] * e1[i] ) * x0[ix] * z0[iz] ;
                    g01[i+3*nCFuncI] = ( z1[iz] * e0[i] + z0[iz+1] * e1[i] ) * x0[ix] * y0[iy] ;
                }
                /* . Transform the integrals. */
                values = g01 ;
                work   = gT  ;
                GaussianBasisTransform1M ( nCFuncI, 4, iBasis->shells[iShell].c2s, True, &values, &work ) ;
                pG01 = values ;
                /* . Put the values in the proper place. */
                iStart = iBasis->shells[iShell].nStart ;
                nFuncI = iBasis->shells[iShell].nBasis ;
                for ( i = 0 ; i < nFuncI ; i++ )
                {
                    Array2D_Item ( f , i+iStart, g ) = pG01[i] ;
                    Array2D_Item ( fX, i+iStart, g ) = pG01[i+   nFuncI] ;
                    Array2D_Item ( fY, i+iStart, g ) = pG01[i+ 2*nFuncI] ;
                    Array2D_Item ( fZ, i+iStart, g ) = pG01[i+ 3*nFuncI] ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at the given points and their first and second derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 23 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Op1ir12 ( const GaussianBasis *iBasis ,
                                        const Real          *rI     ,
                                        const Coordinates3  *rG     ,
                                        const Integer        s1     ,
                                              Real          *rWork  ,
                                              RealArray2D   *f      ,
                                              RealArray2D   *fX     ,
                                              RealArray2D   *fY     ,
                                              RealArray2D   *fZ     ,
                                              RealArray2D   *fXX    ,
                                              RealArray2D   *fXY    ,
                                              RealArray2D   *fXZ    ,
                                              RealArray2D   *fYY    ,
                                              RealArray2D   *fYZ    ,
                                              RealArray2D   *fZZ    )
{
    if ( ( iBasis != NULL ) && 
         ( rI     != NULL ) && 
         ( rG     != NULL ) && 
         ( f      != NULL ) && 
         ( fX     != NULL ) && 
         ( fY     != NULL ) && 
         ( fZ     != NULL ) && 
         ( fXX    != NULL ) && 
         ( fXY    != NULL ) && 
         ( fXZ    != NULL ) && 
         ( fYY    != NULL ) && 
         ( fYZ    != NULL ) && 
         ( fZZ    != NULL ) )
    {
        Integer  g, i, ip, iShell, iStart, ix, iy, iz, nCFuncI, nFuncI ;
        Real     dx, dy, dz, e, eip, ee, r2 ;
        Real    *pG012, *values = NULL, *work = NULL ;     
        Real    *e0, *e1, *e2, *g012, *gT,
                 x0[MAXAMP3], y0[MAXAMP3], z0[MAXAMP3] ,
                 x1[MAXAMP3], y1[MAXAMP3], z1[MAXAMP3] ,
                 x2[MAXAMP3], y2[MAXAMP3], z2[MAXAMP3] ;
        /* . Set pointers. */
        e0 = &rWork[0] ; e1 = &rWork[s1] ; e2 = &rWork[2*s1] ; g012 = &rWork[3*s1] ; gT = &rWork[13*s1] ; /* . Note 10 each for g012 and gT! */
        /* . Loop over points. */
        for ( g = 0 ; g < Coordinates3_Rows ( rG ) ; g++ )
        {
            /* . Get the distance squared between atom centers. */
            dx = Coordinates3_Item ( rG, g, 0 ) - rI[0] ;
            dy = Coordinates3_Item ( rG, g, 1 ) - rI[1] ;
            dz = Coordinates3_Item ( rG, g, 2 ) - rI[2] ;
            r2 = dx * dx + dy * dy + dz * dz ;
            /* . Check for negligible contributions. */
            /* if ( r2 >= iBasis->maximumRange ) return ; */
            /* . Form the angular functions. */
            x0[0] = 1.0e+00 ; y0[0] = 1.0e+00 ; z0[0] = 1.0e+00 ;
            x1[0] = 0.0e+00 ; y1[0] = 0.0e+00 ; z1[0] = 0.0e+00 ;
            x2[0] = 0.0e+00 ; y2[0] = 0.0e+00 ; z2[0] = 0.0e+00 ;
            for ( i = 1 ; i <= iBasis->lHigh + 2; i++ )
            {
                x0[i] = dx * x0[i-1] ;
                y0[i] = dy * y0[i-1] ;
                z0[i] = dz * z0[i-1] ;
                x1[i] = ( Real ) ( i ) * x0[i-1] ;
                y1[i] = ( Real ) ( i ) * y0[i-1] ;
                z1[i] = ( Real ) ( i ) * z0[i-1] ;
                x2[i] = ( Real ) ( i ) * x1[i-1] ;
                y2[i] = ( Real ) ( i ) * y1[i-1] ;
                z2[i] = ( Real ) ( i ) * z1[i-1] ;
            }
            /* . Outer loop over shells. */
            for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
            {
                /* . Get information about the shell. */
                nCFuncI = iBasis->shells[iShell].nCBF     ;
                /* . Check for negligible contributions. */
                /* if ( r2 >= iBasis->shells[iShell].maximumRange ) return ; */
                /* . Form the exponential factors. */
                for ( i = 0 ; i < nCFuncI ; i++ ) { e0[i] = 0.0e+00 ; e1[i] = 0.0e+00 ; e2[i] = 0.0e+00 ; }
                for ( ip = 0; ip < iBasis->shells[iShell].nPrimitives ; ip ++ )
                {
                    e   = iBasis->shells[iShell].primitives[ip].exponent ;
                    ee  = e * e ;
                    eip = exp ( - e * r2 ) ;
                    for ( i = 0 ; i < nCFuncI ; i++ )
                    {
                        e0[i] += iBasis->shells[iShell].primitives[ip].cCBF[i]       * eip ;
                        e1[i] += iBasis->shells[iShell].primitives[ip].cCBF[i]  * e  * eip ;
                        e2[i] += iBasis->shells[iShell].primitives[ip].cCBF[i]  * ee * eip ;
                    }
                }
                for ( i = 0 ; i < nCFuncI ; i++ ) { e1[i] *= -2.0e+00 ; e2[i] *= 4.0e+00 ; }
                /* . Form the Cartesian function values. */
                for ( i = 0 ; i < nCFuncI ; i++ )
                {
                    ix = iBasis->shells[iShell].cbfPowX[i] ;
                    iy = iBasis->shells[iShell].cbfPowY[i] ;
                    iz = iBasis->shells[iShell].cbfPowZ[i] ;
                    g012[i]          = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                    g012[i+  nCFuncI] = ( x1[ix] * e0[i] + x0[ix+1] * e1[i] ) * y0[iy] * z0[iz] ;
                    g012[i+2*nCFuncI] = ( y1[iy] * e0[i] + y0[iy+1] * e1[i] ) * x0[ix] * z0[iz] ;
                    g012[i+3*nCFuncI] = ( z1[iz] * e0[i] + z0[iz+1] * e1[i] ) * x0[ix] * y0[iy] ;
                    g012[i+4*nCFuncI] = ( x2[ix] * e0[i] + ( dx * x1[ix] + x1[ix+1] ) * e1[i] + x0[ix+2] * e2[i] ) * y0[iy] * z0[iz] ;
                    g012[i+5*nCFuncI] = ( x1[ix]*y1[iy]*e0[i] + ( x1[ix]*y0[iy+1] + x0[ix+1]*y1[iy] )*e1[i] + x0[ix+1]*y0[iy+1]*e2[i] ) * z0[iz] ;
                    g012[i+6*nCFuncI] = ( x1[ix]*z1[iz]*e0[i] + ( x1[ix]*z0[iz+1] + x0[ix+1]*z1[iz] )*e1[i] + x0[ix+1]*z0[iz+1]*e2[i] ) * y0[iy] ;
                    g012[i+7*nCFuncI] = ( y2[iy] * e0[i] + ( dy * y1[iy] + y1[iy+1] ) * e1[i] + y0[iy+2] * e2[i] ) * x0[ix] * z0[iz] ;
                    g012[i+8*nCFuncI] = ( y1[iy]*z1[iz]*e0[i] + ( y1[iy]*z0[iz+1] + y0[iy+1]*z1[iz] )*e1[i] + y0[iy+1]*z0[iz+1]*e2[i] ) * x0[ix] ;
                    g012[i+9*nCFuncI] = ( z2[iz] * e0[i] + ( dz * z1[iz] + z1[iz+1] ) * e1[i] + z0[iz+2] * e2[i] ) * x0[ix] * y0[iy] ;
                }
                /* . Transform the integrals. */
                values = g012 ;
                work   = gT   ;
                GaussianBasisTransform1M ( nCFuncI, 10, iBasis->shells[iShell].c2s, True, &values, &work ) ;
                pG012 = values ;
                /* . Put the values in the proper place. */
                iStart = iBasis->shells[iShell].nStart ;
                nFuncI = iBasis->shells[iShell].nBasis ;
                for ( i = 0 ; i < nFuncI ; i++ )
                {
                    Array2D_Item ( f  , i+iStart, g ) = pG012[i] ;
                    Array2D_Item ( fX , i+iStart, g ) = pG012[i+   nFuncI] ;
                    Array2D_Item ( fY , i+iStart, g ) = pG012[i+ 2*nFuncI] ;
                    Array2D_Item ( fZ , i+iStart, g ) = pG012[i+ 3*nFuncI] ;
                    Array2D_Item ( fXX, i+iStart, g ) = pG012[i+ 4*nFuncI] ;
                    Array2D_Item ( fXY, i+iStart, g ) = pG012[i+ 5*nFuncI] ;
                    Array2D_Item ( fXZ, i+iStart, g ) = pG012[i+ 6*nFuncI] ;
                    Array2D_Item ( fYY, i+iStart, g ) = pG012[i+ 7*nFuncI] ;
                    Array2D_Item ( fYZ, i+iStart, g ) = pG012[i+ 8*nFuncI] ;
                    Array2D_Item ( fZZ, i+iStart, g ) = pG012[i+ 9*nFuncI] ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at the given points and their first, second and third derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: real 44 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Op1ir123 ( const GaussianBasis *iBasis ,
                                         const Real          *rI     ,
                                         const Coordinates3  *rG     ,
                                         const Integer        s1     ,
                                               Real          *rWork  ,
                                               RealArray2D   *f      ,
                                               RealArray2D   *fX     ,
                                               RealArray2D   *fY     ,
                                               RealArray2D   *fZ     ,
                                               RealArray2D   *fXX    ,
                                               RealArray2D   *fXY    ,
                                               RealArray2D   *fXZ    ,
                                               RealArray2D   *fYY    ,
                                               RealArray2D   *fYZ    ,
                                               RealArray2D   *fZZ    ,
                                               RealArray2D   *fXXX   ,
                                               RealArray2D   *fXXY   ,
                                               RealArray2D   *fXXZ   ,
                                               RealArray2D   *fXYY   ,
                                               RealArray2D   *fXYZ   ,
                                               RealArray2D   *fXZZ   ,
                                               RealArray2D   *fYYY   ,
                                               RealArray2D   *fYYZ   ,
                                               RealArray2D   *fYZZ   ,
                                               RealArray2D   *fZZZ   )
{
    if ( ( iBasis != NULL ) && 
         ( rI     != NULL ) && 
         ( rG     != NULL ) && 
         ( f      != NULL ) && 
         ( fX     != NULL ) && 
         ( fY     != NULL ) && 
         ( fZ     != NULL ) && 
         ( fXX    != NULL ) && 
         ( fXY    != NULL ) && 
         ( fXZ    != NULL ) && 
         ( fYY    != NULL ) && 
         ( fYZ    != NULL ) && 
         ( fZZ    != NULL ) && 
         ( fXXX   != NULL ) && 
         ( fXXY   != NULL ) && 
         ( fXXZ   != NULL ) && 
         ( fXYY   != NULL ) && 
         ( fXYZ   != NULL ) && 
         ( fXZZ   != NULL ) && 
         ( fYYY   != NULL ) && 
         ( fYYZ   != NULL ) && 
         ( fYZZ   != NULL ) && 
         ( fZZZ   != NULL ) )
    {
        Integer  g, i, ip, iShell, iStart, ix, iy, iz, nCFuncI, nFuncI ;
        Real     dx, dy, dz, e, eip, ee, eee, r2 ;
        Real    *pG013, *values = NULL, *work = NULL ;     
        Real    *e0, *e1, *e2, *e3, *g013, *gT,
                 x0[MAXAMP4], y0[MAXAMP4], z0[MAXAMP4] ,
                 x1[MAXAMP4], y1[MAXAMP4], z1[MAXAMP4] ,
                 x2[MAXAMP4], y2[MAXAMP4], z2[MAXAMP4] ,
                 x3[MAXAMP4], y3[MAXAMP4], z3[MAXAMP4] ;
        /* . Set pointers. */
        e0   = &rWork[   0] ; e1 = &rWork[   s1] ; e2 = &rWork[2*s1] ; e3 = &rWork[3*s1] ;
        g013 = &rWork[4*s1] ; gT = &rWork[24*s1] ; /* . Note 20 each for g013 and gT! */
        /* . Loop over points. */
        for ( g = 0 ; g < Coordinates3_Rows ( rG ) ; g++ )
        {
            /* . Get the distance squared between atom centers. */
            dx = Coordinates3_Item ( rG, g, 0 ) - rI[0] ;
            dy = Coordinates3_Item ( rG, g, 1 ) - rI[1] ;
            dz = Coordinates3_Item ( rG, g, 2 ) - rI[2] ;
            r2 = dx * dx + dy * dy + dz * dz ;
            /* . Check for negligible contributions. */
            /* if ( r2 >= iBasis->maximumRange ) return ; */
            /* . Form the angular functions. */
            x0[0] = 1.0e+00 ; y0[0] = 1.0e+00 ; z0[0] = 1.0e+00 ;
            x1[0] = 0.0e+00 ; y1[0] = 0.0e+00 ; z1[0] = 0.0e+00 ;
            x2[0] = 0.0e+00 ; y2[0] = 0.0e+00 ; z2[0] = 0.0e+00 ;
            x3[0] = 0.0e+00 ; y3[0] = 0.0e+00 ; z3[0] = 0.0e+00 ;
            for ( i = 1 ; i <= iBasis->lHigh + 3 ; i++ )
            {
                x0[i] = dx * x0[i-1] ;
                y0[i] = dy * y0[i-1] ;
                z0[i] = dz * z0[i-1] ;
                x1[i] = ( Real ) ( i ) * x0[i-1] ;
                y1[i] = ( Real ) ( i ) * y0[i-1] ;
                z1[i] = ( Real ) ( i ) * z0[i-1] ;
                x2[i] = ( Real ) ( i ) * x1[i-1] ;
                y2[i] = ( Real ) ( i ) * y1[i-1] ;
                z2[i] = ( Real ) ( i ) * z1[i-1] ;
                x3[i] = ( Real ) ( i ) * x2[i-1] ;
                y3[i] = ( Real ) ( i ) * y2[i-1] ;
                z3[i] = ( Real ) ( i ) * z2[i-1] ;
            }
            /* . Outer loop over shells. */
            for ( iShell = 0 ; iShell < iBasis->nShells ; iShell++ )
            {
                /* . Get information about the shell. */
                nCFuncI = iBasis->shells[iShell].nCBF     ;
                /* . Check for negligible contributions. */
                /* if ( r2 >= iBasis->shells[iShell].maximumRange ) return ; */
                /* . Form the exponential factors. */
                for ( i = 0 ; i < nCFuncI ; i++ ) { e0[i] = 0.0e+00 ; e1[i] = 0.0e+00 ; e2[i] = 0.0e+00 ; e3[i] = 0.0e+00 ; }
                for ( ip = 0; ip < iBasis->shells[iShell].nPrimitives ; ip ++ )
                {
                    e   = iBasis->shells[iShell].primitives[ip].exponent ;
                    ee  = e  * e ;
                    eee = ee * e ;
                    eip = exp ( - e * r2 ) ;
                    for ( i = 0 ; i < nCFuncI ; i++ )
                    {
                        e0[i] += iBasis->shells[iShell].primitives[ip].cCBF[i]        * eip ;
                        e1[i] += iBasis->shells[iShell].primitives[ip].cCBF[i]  * e   * eip ;
                        e2[i] += iBasis->shells[iShell].primitives[ip].cCBF[i]  * ee  * eip ;
                        e3[i] += iBasis->shells[iShell].primitives[ip].cCBF[i]  * eee * eip ;
                    }
                }
                for ( i = 0 ; i < nCFuncI ; i++ ) { e1[i] *= -2.0e+00 ; e2[i] *= 4.0e+00 ; e3[i] *= -8.0e+00 ; }
                /* . Form the Cartesian function values. */
                for ( i = 0 ; i < nCFuncI ; i++ )
                {
                    ix = iBasis->shells[iShell].cbfPowX[i] ;
                    iy = iBasis->shells[iShell].cbfPowY[i] ;
                    iz = iBasis->shells[iShell].cbfPowZ[i] ;
                    g013[i]          = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                    g013[i+   nCFuncI] = ( x1[ix] * e0[i] + x0[ix+1] * e1[i] ) * y0[iy] * z0[iz] ;
                    g013[i+ 2*nCFuncI] = ( y1[iy] * e0[i] + y0[iy+1] * e1[i] ) * x0[ix] * z0[iz] ;
                    g013[i+ 3*nCFuncI] = ( z1[iz] * e0[i] + z0[iz+1] * e1[i] ) * x0[ix] * y0[iy] ;
                    g013[i+ 4*nCFuncI] = ( x2[ix] * e0[i] + ( dx * x1[ix] + x1[ix+1] ) * e1[i] + x0[ix+2] * e2[i] ) * y0[iy] * z0[iz] ;
                    g013[i+ 5*nCFuncI] = ( x1[ix]*y1[iy]*e0[i] + ( x1[ix]*y0[iy+1] + x0[ix+1]*y1[iy] )*e1[i] + x0[ix+1]*y0[iy+1]*e2[i] ) * z0[iz] ;
                    g013[i+ 6*nCFuncI] = ( x1[ix]*z1[iz]*e0[i] + ( x1[ix]*z0[iz+1] + x0[ix+1]*z1[iz] )*e1[i] + x0[ix+1]*z0[iz+1]*e2[i] ) * y0[iy] ;
                    g013[i+ 7*nCFuncI] = ( y2[iy] * e0[i] + ( dy * y1[iy] + y1[iy+1] ) * e1[i] + y0[iy+2] * e2[i] ) * x0[ix] * z0[iz] ;
                    g013[i+ 8*nCFuncI] = ( y1[iy]*z1[iz]*e0[i] + ( y1[iy]*z0[iz+1] + y0[iy+1]*z1[iz] )*e1[i] + y0[iy+1]*z0[iz+1]*e2[i] ) * x0[ix] ;
                    g013[i+ 9*nCFuncI] = ( z2[iz] * e0[i] + ( dz * z1[iz] + z1[iz+1] ) * e1[i] + z0[iz+2] * e2[i] ) * x0[ix] * y0[iy] ;
                    g013[i+10*nCFuncI] = ( x3[ix] * e0[i] + ( x1[ix] + 2.0e+00*dx*x2[ix] + x2[ix+1] )*e1[i] + ( dx*dx*x1[ix] + dx*x1[ix+1] + x1[ix+2] )*e2[i] + x0[ix+3] * e3[i] ) * y0[iy] * z0[iz] ;
                    g013[i+11*nCFuncI] = ( x2[ix]*y1[iy]*e0[i] + ( ( dx*x1[ix] + x1[ix+1] )*y1[iy]   + x2[ix]*y0[iy+1] )*e1[i] +
                                                          ( ( dx*x1[ix] + x1[ix+1] )*y0[iy+1] + x0[ix+2]*y1[iy] )*e2[i] + x0[ix+2]*y0[iy+1]*e3[i] ) * z0[iz] ;
                    g013[i+12*nCFuncI] = ( x2[ix]*z1[iz]*e0[i] + ( ( dx*x1[ix] + x1[ix+1] )*z1[iz]   + x2[ix]*z0[iz+1] )*e1[i] +
                                                          ( ( dx*x1[ix] + x1[ix+1] )*z0[iz+1] + x0[ix+2]*z1[iz] )*e2[i] + x0[ix+2]*z0[iz+1]*e3[i] ) * y0[iy] ;
                    g013[i+13*nCFuncI] = ( y2[iy]*x1[ix]*e0[i] + ( ( dy*y1[iy] + y1[iy+1] )*x1[ix]   + y2[iy]*x0[ix+1] )*e1[i] +
                                                          ( ( dy*y1[iy] + y1[iy+1] )*x0[ix+1] + y0[iy+2]*x1[ix] )*e2[i] + y0[iy+2]*x0[ix+1]*e3[i] ) * z0[iz] ;
                    g013[i+14*nCFuncI] = ( x1[ix]*y1[iy]*z1[iz]*e0[i] + ( x1[ix]*y1[iy  ]*z0[iz+1] + x1[ix  ]*y0[iy+1]*z1[iz  ] + x0[ix+1]*y1[iy  ]*z1[iz] )*e1[i] +
                                                                 ( x1[ix]*y0[iy+1]*z0[iz+1] + x0[ix+1]*y1[iy  ]*z0[iz+1] + x0[ix+1]*y0[iy+1]*z1[iz] )*e2[i] + x0[ix+1]*y0[ix+1]*z0[iz+1]*e3[i] ) ;
                    g013[i+15*nCFuncI] = ( z2[iz]*x1[ix]*e0[i] + ( ( dz*z1[iz] + z1[iz+1] )*x1[ix]   + z2[iz]*x0[ix+1] )*e1[i] +
                                                          ( ( dz*z1[iz] + z1[iz+1] )*x0[ix+1] + z0[iz+2]*x1[ix] )*e2[i] + z0[iz+2]*x0[ix+1]*e3[i] ) * y0[iy] ;
                    g013[i+16*nCFuncI] = ( y3[iy] * e0[i] + ( y1[iy] + 2.0e+00*dy*y2[iy] + y2[iy+1] )*e1[i] + ( dy*dy*y1[iy] + dy*y1[iy+1] + y1[iy+2] )*e2[i] + y0[iy+3] * e3[i] ) * x0[ix] * z0[iz] ;
                    g013[i+17*nCFuncI] = ( y2[iy]*z1[iz]*e0[i] + ( ( dy*y1[iy] + y1[iy+1] )*z1[iz]   + y2[iy]*z0[iz+1] )*e1[i] +
                                                          ( ( dy*y1[iy] + y1[iy+1] )*z0[iz+1] + y0[iy+2]*z1[iz] )*e2[i] + y0[iy+2]*z0[iz+1]*e3[i] ) * x0[ix] ;
                    g013[i+18*nCFuncI] = ( z2[iz]*y1[iy]*e0[i] + ( ( dz*z1[iz] + z1[iz+1] )*y1[iy]   + z2[iz]*y0[iy+1] )*e1[i] +
                                                          ( ( dz*z1[iz] + z1[iz+1] )*y0[iy+1] + z0[iz+2]*y1[iy] )*e2[i] + z0[iz+2]*y0[iy+1]*e3[i] ) * x0[ix] ;
                    g013[i+19*nCFuncI] = ( z3[iz] * e0[i] + ( z1[iz] + 2.0e+00*dz*z2[iz] + z2[iz+1] )*e1[i] + ( dz*dz*z1[iz] + dz*z1[iz+1] + z1[iz+2] )*e2[i] + z0[iz+3] * e3[i] ) * x0[ix] * y0[iy] ;
                }
                /* . Transform the integrals. */
                values = g013 ;
                work   = gT   ;
                GaussianBasisTransform1M ( nCFuncI, 20, iBasis->shells[iShell].c2s, True, &values, &work ) ;
                pG013 = values ;
                /* . Put the values in the proper place. */
                iStart = iBasis->shells[iShell].nStart ;
                nFuncI = iBasis->shells[iShell].nBasis ;
                for ( i = 0 ; i < nFuncI ; i++ )
                {
                    Array2D_Item ( f   , i+iStart, g ) = pG013[i] ;
                    Array2D_Item ( fX  , i+iStart, g ) = pG013[i+   nFuncI] ;
                    Array2D_Item ( fY  , i+iStart, g ) = pG013[i+ 2*nFuncI] ;
                    Array2D_Item ( fZ  , i+iStart, g ) = pG013[i+ 3*nFuncI] ;
                    Array2D_Item ( fXX , i+iStart, g ) = pG013[i+ 4*nFuncI] ;
                    Array2D_Item ( fXY , i+iStart, g ) = pG013[i+ 5*nFuncI] ;
                    Array2D_Item ( fXZ , i+iStart, g ) = pG013[i+ 6*nFuncI] ;
                    Array2D_Item ( fYY , i+iStart, g ) = pG013[i+ 7*nFuncI] ;
                    Array2D_Item ( fYZ , i+iStart, g ) = pG013[i+ 8*nFuncI] ;
                    Array2D_Item ( fZZ , i+iStart, g ) = pG013[i+ 9*nFuncI] ;
                    Array2D_Item ( fXXX, i+iStart, g ) = pG013[i+10*nFuncI] ;
                    Array2D_Item ( fXXY, i+iStart, g ) = pG013[i+11*nFuncI] ;
                    Array2D_Item ( fXXZ, i+iStart, g ) = pG013[i+12*nFuncI] ;
                    Array2D_Item ( fXYY, i+iStart, g ) = pG013[i+13*nFuncI] ;
                    Array2D_Item ( fXYZ, i+iStart, g ) = pG013[i+14*nFuncI] ;
                    Array2D_Item ( fXZZ, i+iStart, g ) = pG013[i+15*nFuncI] ;
                    Array2D_Item ( fYYY, i+iStart, g ) = pG013[i+16*nFuncI] ;
                    Array2D_Item ( fYYZ, i+iStart, g ) = pG013[i+17*nFuncI] ;
                    Array2D_Item ( fYZZ, i+iStart, g ) = pG013[i+18*nFuncI] ;
                    Array2D_Item ( fZZZ, i+iStart, g ) = pG013[i+19*nFuncI] ;
                }
            }
        }
    }
}
