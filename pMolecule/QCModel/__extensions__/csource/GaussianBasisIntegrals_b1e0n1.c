/*==================================================================================================================================
! . Integrals - 1 basis, 0 electrons, 1 nucleus/point.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "GaussianBasisIntegrals_b1e0n1.h"

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
void GaussianBasisIntegrals_Grid ( const GaussianBasis *iBasis ,
                                   const Real          *rI     ,
                                   const Coordinates3  *rG     ,
                                         RealArray2D   *f      )
{
    if ( ( iBasis != NULL ) && 
         ( rI     != NULL ) && 
         ( rG     != NULL ) && 
         ( f      != NULL ) )
    {
        Integer  g, i, icbfind, ip, iShell, iStart, ix, iy, iz, ncfunci ;
        Real     dx, dy, dz, eip, r2 ;
        Real     e0[MAXCBF], g0[MAXCBF], x0[MAXAMP1], y0[MAXAMP1], z0[MAXAMP1] ;
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
            for ( i = 1 ; i <= iBasis->maximum_angularmomentum ; i++ )
            {
                x0[i] = dx * x0[i-1] ;
                y0[i] = dy * y0[i-1] ;
                z0[i] = dz * z0[i-1] ;
            }
            /* . Outer loop over shells. */
            for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
            {
                /* . Get information about the shell. */
                icbfind = iBasis->shells[iShell].type->cbfindex ;
                ncfunci = iBasis->shells[iShell].type->ncbf     ;
                /* . Check for negligible contributions. */
                /* if ( r2 >= iBasis->shells[iShell].maximumRange ) return ; */
                /* . Form the exponential factors. */
                for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] = 0.0e+00 ; }
                for ( ip = 0; ip < iBasis->shells[iShell].nprimitives ; ip ++ )
                {
                    eip = exp ( - iBasis->shells[iShell].primitives[ip].exponent * r2 ) ;
                    for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] += iBasis->shells[iShell].primitives[ip].ccbf[i] * eip ; }
                }
                /* . Form the Cartesian function values. */
                for ( i = 0 ; i < ncfunci ; i++ )
                {
                    ix = CBFPOWX[i+icbfind] ;
                    iy = CBFPOWY[i+icbfind] ;
                    iz = CBFPOWZ[i+icbfind] ;
                    g0[i] = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                }
                /* . Put the values in the proper place. */
                iStart = iBasis->shells[iShell].nstartw ;
                for ( i = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ ) { Array2D_Item ( f, i+iStart, g ) = g0[i] ; }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at the given points and their first derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_GridD ( const GaussianBasis *iBasis ,
                                    const Real          *rI     ,
                                    const Coordinates3  *rG     ,
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
        Integer  g, i, icbfind, ip, iShell, iStart, ix, iy, iz, ncfunci ;
        Real     dx, dy, dz, eip, r2 ;
        Real     e0[MAXCBF], e1[MAXCBF], g01[4*MAXCBF], x0[MAXAMP2], y0[MAXAMP2], z0[MAXAMP2] ,
                                                        x1[MAXAMP2], y1[MAXAMP2], z1[MAXAMP2] ;
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
            for ( i = 1 ; i <= iBasis->maximum_angularmomentum + 1; i++ )
            {
                x0[i] = dx * x0[i-1] ;
                y0[i] = dy * y0[i-1] ;
                z0[i] = dz * z0[i-1] ;
                x1[i] = ( Real ) ( i ) * x0[i-1] ;
                y1[i] = ( Real ) ( i ) * y0[i-1] ;
                z1[i] = ( Real ) ( i ) * z0[i-1] ;
            }
            /* . Outer loop over shells. */
            for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
            {
                /* . Get information about the shell. */
                icbfind = iBasis->shells[iShell].type->cbfindex ;
                ncfunci = iBasis->shells[iShell].type->ncbf     ;
                /* . Check for negligible contributions. */
                /* if ( r2 >= iBasis->shells[iShell].maximumRange ) return ; */
                /* . Form the exponential factors. */
                for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] = 0.0e+00 ; e1[i] = 0.0e+00 ; }
                for ( ip = 0; ip < iBasis->shells[iShell].nprimitives ; ip ++ )
                {
                    eip = exp ( - iBasis->shells[iShell].primitives[ip].exponent * r2 ) ;
                    for ( i = 0 ; i < ncfunci ; i++ )
                    {
                        e0[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * eip ;
                        e1[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  *
                                 iBasis->shells[iShell].primitives[ip].exponent * eip ;
                    }
                }
                for ( i = 0 ; i < ncfunci ; i++ ) { e1[i] *= -2.0e+00 ; }
                /* . Form the Cartesian function values. */
                for ( i = 0 ; i < ncfunci ; i++ )
                {
                    ix = CBFPOWX[i+icbfind] ;
                    iy = CBFPOWY[i+icbfind] ;
                    iz = CBFPOWZ[i+icbfind] ;
                    g01[i]          = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                    g01[i+  MAXCBF] = ( x1[ix] * e0[i] + x0[ix+1] * e1[i] ) * y0[iy] * z0[iz] ;
                    g01[i+2*MAXCBF] = ( y1[iy] * e0[i] + y0[iy+1] * e1[i] ) * x0[ix] * z0[iz] ;
                    g01[i+3*MAXCBF] = ( z1[iz] * e0[i] + z0[iz+1] * e1[i] ) * x0[ix] * y0[iy] ;
                }
                /* . Put the values in the proper place. */
                iStart = iBasis->shells[iShell].nstartw ;
                for ( i = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
                {
                    Array2D_Item ( f , i+iStart, g ) = g01[i] ;
                    Array2D_Item ( fX, i+iStart, g ) = g01[i+   MAXCBF] ;
                    Array2D_Item ( fY, i+iStart, g ) = g01[i+ 2*MAXCBF] ;
                    Array2D_Item ( fZ, i+iStart, g ) = g01[i+ 3*MAXCBF] ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at the given points and their first and second derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_GridD2 ( const GaussianBasis *iBasis ,
                                     const Real          *rI     ,
                                     const Coordinates3  *rG     ,
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
        Integer  g, i, icbfind, ip, iShell, iStart, ix, iy, iz, ncfunci ;
        Real     dx, dy, dz, e, eip, ee, r2 ;
        Real     e0[MAXCBF], e1[MAXCBF], e2[MAXCBF], g012[10*MAXCBF], x0[MAXAMP3], y0[MAXAMP3], z0[MAXAMP3] ,
                                                                      x1[MAXAMP3], y1[MAXAMP3], z1[MAXAMP3] ,
                                                                      x2[MAXAMP3], y2[MAXAMP3], z2[MAXAMP3] ;
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
            for ( i = 1 ; i <= iBasis->maximum_angularmomentum + 2; i++ )
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
            for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
            {
                /* . Get information about the shell. */
                icbfind = iBasis->shells[iShell].type->cbfindex ;
                ncfunci = iBasis->shells[iShell].type->ncbf     ;
                /* . Check for negligible contributions. */
                /* if ( r2 >= iBasis->shells[iShell].maximumRange ) return ; */
                /* . Form the exponential factors. */
                for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] = 0.0e+00 ; e1[i] = 0.0e+00 ; e2[i] = 0.0e+00 ; }
                for ( ip = 0; ip < iBasis->shells[iShell].nprimitives ; ip ++ )
                {
                    e   = iBasis->shells[iShell].primitives[ip].exponent ;
                    ee  = e * e ;
                    eip = exp ( - e * r2 ) ;
                    for ( i = 0 ; i < ncfunci ; i++ )
                    {
                        e0[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]       * eip ;
                        e1[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * e  * eip ;
                        e2[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * ee * eip ;
                    }
                }
                for ( i = 0 ; i < ncfunci ; i++ ) { e1[i] *= -2.0e+00 ; e2[i] *= 4.0e+00 ; }
                /* . Form the Cartesian function values. */
                for ( i = 0 ; i < ncfunci ; i++ )
                {
                    ix = CBFPOWX[i+icbfind] ;
                    iy = CBFPOWY[i+icbfind] ;
                    iz = CBFPOWZ[i+icbfind] ;
                    g012[i]          = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                    g012[i+  MAXCBF] = ( x1[ix] * e0[i] + x0[ix+1] * e1[i] ) * y0[iy] * z0[iz] ;
                    g012[i+2*MAXCBF] = ( y1[iy] * e0[i] + y0[iy+1] * e1[i] ) * x0[ix] * z0[iz] ;
                    g012[i+3*MAXCBF] = ( z1[iz] * e0[i] + z0[iz+1] * e1[i] ) * x0[ix] * y0[iy] ;
                    g012[i+4*MAXCBF] = ( x2[ix] * e0[i] + ( dx * x1[ix] + x1[ix+1] ) * e1[i] + x0[ix+2] * e2[i] ) * y0[iy] * z0[iz] ;
                    g012[i+5*MAXCBF] = ( x1[ix]*y1[iy]*e0[i] + ( x1[ix]*y0[iy+1] + x0[ix+1]*y1[iy] )*e1[i] + x0[ix+1]*y0[iy+1]*e2[i] ) * z0[iz] ;
                    g012[i+6*MAXCBF] = ( x1[ix]*z1[iz]*e0[i] + ( x1[ix]*z0[iz+1] + x0[ix+1]*z1[iz] )*e1[i] + x0[ix+1]*z0[iz+1]*e2[i] ) * y0[iy] ;
                    g012[i+7*MAXCBF] = ( y2[iy] * e0[i] + ( dy * y1[iy] + y1[iy+1] ) * e1[i] + y0[iy+2] * e2[i] ) * x0[ix] * z0[iz] ;
                    g012[i+8*MAXCBF] = ( y1[iy]*z1[iz]*e0[i] + ( y1[iy]*z0[iz+1] + y0[iy+1]*z1[iz] )*e1[i] + y0[iy+1]*z0[iz+1]*e2[i] ) * x0[ix] ;
                    g012[i+9*MAXCBF] = ( z2[iz] * e0[i] + ( dz * z1[iz] + z1[iz+1] ) * e1[i] + z0[iz+2] * e2[i] ) * x0[ix] * y0[iy] ;
                }
                /* . Put the values in the proper place. */
                iStart = iBasis->shells[iShell].nstartw ;
                for ( i = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
                {
                    Array2D_Item ( f  , i+iStart, g ) = g012[i] ;
                    Array2D_Item ( fX , i+iStart, g ) = g012[i+   MAXCBF] ;
                    Array2D_Item ( fY , i+iStart, g ) = g012[i+ 2*MAXCBF] ;
                    Array2D_Item ( fZ , i+iStart, g ) = g012[i+ 3*MAXCBF] ;
                    Array2D_Item ( fXX, i+iStart, g ) = g012[i+ 4*MAXCBF] ;
                    Array2D_Item ( fXY, i+iStart, g ) = g012[i+ 5*MAXCBF] ;
                    Array2D_Item ( fXZ, i+iStart, g ) = g012[i+ 6*MAXCBF] ;
                    Array2D_Item ( fYY, i+iStart, g ) = g012[i+ 7*MAXCBF] ;
                    Array2D_Item ( fYZ, i+iStart, g ) = g012[i+ 8*MAXCBF] ;
                    Array2D_Item ( fZZ, i+iStart, g ) = g012[i+ 9*MAXCBF] ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of the basis functions at the given points and their first, second and third derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_GridD3 ( const GaussianBasis *iBasis ,
                                     const Real          *rI     ,
                                     const Coordinates3  *rG     ,
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
        Integer  g, i, icbfind, ip, iShell, iStart, ix, iy, iz, ncfunci ;
        Real     dx, dy, dz, e, eip, ee, eee, r2 ;
        Real     e0[MAXCBF], e1[MAXCBF], e2[MAXCBF], e3[MAXCBF], g013[20*MAXCBF], x0[MAXAMP4], y0[MAXAMP4], z0[MAXAMP4] ,
                                                                                  x1[MAXAMP4], y1[MAXAMP4], z1[MAXAMP4] ,
                                                                                  x2[MAXAMP4], y2[MAXAMP4], z2[MAXAMP4] ,
                                                                                  x3[MAXAMP4], y3[MAXAMP4], z3[MAXAMP4] ;
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
            for ( i = 1 ; i <= iBasis->maximum_angularmomentum + 3 ; i++ )
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
            for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
            {
                /* . Get information about the shell. */
                icbfind = iBasis->shells[iShell].type->cbfindex ;
                ncfunci = iBasis->shells[iShell].type->ncbf     ;
                /* . Check for negligible contributions. */
                /* if ( r2 >= iBasis->shells[iShell].maximumRange ) return ; */
                /* . Form the exponential factors. */
                for ( i = 0 ; i < ncfunci ; i++ ) { e0[i] = 0.0e+00 ; e1[i] = 0.0e+00 ; e2[i] = 0.0e+00 ; e3[i] = 0.0e+00 ; }
                for ( ip = 0; ip < iBasis->shells[iShell].nprimitives ; ip ++ )
                {
                    e   = iBasis->shells[iShell].primitives[ip].exponent ;
                    ee  = e  * e ;
                    eee = ee * e ;
                    eip = exp ( - e * r2 ) ;
                    for ( i = 0 ; i < ncfunci ; i++ )
                    {
                        e0[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]        * eip ;
                        e1[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * e   * eip ;
                        e2[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * ee  * eip ;
                        e3[i] += iBasis->shells[iShell].primitives[ip].ccbf[i]  * eee * eip ;
                    }
                }
                for ( i = 0 ; i < ncfunci ; i++ ) { e1[i] *= -2.0e+00 ; e2[i] *= 4.0e+00 ; e3[i] *= -8.0e+00 ; }
                /* . Form the Cartesian function values. */
                for ( i = 0 ; i < ncfunci ; i++ )
                {
                    ix = CBFPOWX[i+icbfind] ;
                    iy = CBFPOWY[i+icbfind] ;
                    iz = CBFPOWZ[i+icbfind] ;
                    g013[i]          = x0[ix] * y0[iy] * z0[iz] * e0[i] ;
                    g013[i+   MAXCBF] = ( x1[ix] * e0[i] + x0[ix+1] * e1[i] ) * y0[iy] * z0[iz] ;
                    g013[i+ 2*MAXCBF] = ( y1[iy] * e0[i] + y0[iy+1] * e1[i] ) * x0[ix] * z0[iz] ;
                    g013[i+ 3*MAXCBF] = ( z1[iz] * e0[i] + z0[iz+1] * e1[i] ) * x0[ix] * y0[iy] ;
                    g013[i+ 4*MAXCBF] = ( x2[ix] * e0[i] + ( dx * x1[ix] + x1[ix+1] ) * e1[i] + x0[ix+2] * e2[i] ) * y0[iy] * z0[iz] ;
                    g013[i+ 5*MAXCBF] = ( x1[ix]*y1[iy]*e0[i] + ( x1[ix]*y0[iy+1] + x0[ix+1]*y1[iy] )*e1[i] + x0[ix+1]*y0[iy+1]*e2[i] ) * z0[iz] ;
                    g013[i+ 6*MAXCBF] = ( x1[ix]*z1[iz]*e0[i] + ( x1[ix]*z0[iz+1] + x0[ix+1]*z1[iz] )*e1[i] + x0[ix+1]*z0[iz+1]*e2[i] ) * y0[iy] ;
                    g013[i+ 7*MAXCBF] = ( y2[iy] * e0[i] + ( dy * y1[iy] + y1[iy+1] ) * e1[i] + y0[iy+2] * e2[i] ) * x0[ix] * z0[iz] ;
                    g013[i+ 8*MAXCBF] = ( y1[iy]*z1[iz]*e0[i] + ( y1[iy]*z0[iz+1] + y0[iy+1]*z1[iz] )*e1[i] + y0[iy+1]*z0[iz+1]*e2[i] ) * x0[ix] ;
                    g013[i+ 9*MAXCBF] = ( z2[iz] * e0[i] + ( dz * z1[iz] + z1[iz+1] ) * e1[i] + z0[iz+2] * e2[i] ) * x0[ix] * y0[iy] ;
                    g013[i+10*MAXCBF] = ( x3[ix] * e0[i] + ( x1[ix] + 2.0e+00*dx*x2[ix] + x2[ix+1] )*e1[i] + ( dx*dx*x1[ix] + dx*x1[ix+1] + x1[ix+2] )*e2[i] + x0[ix+3] * e3[i] ) * y0[iy] * z0[iz] ;
                    g013[i+11*MAXCBF] = ( x2[ix]*y1[iy]*e0[i] + ( ( dx*x1[ix] + x1[ix+1] )*y1[iy]   + x2[ix]*y0[iy+1] )*e1[i] +
                                                                ( ( dx*x1[ix] + x1[ix+1] )*y0[iy+1] + x0[ix+2]*y1[iy] )*e2[i] + x0[ix+2]*y0[iy+1]*e3[i] ) * z0[iz] ;
                    g013[i+12*MAXCBF] = ( x2[ix]*z1[iz]*e0[i] + ( ( dx*x1[ix] + x1[ix+1] )*z1[iz]   + x2[ix]*z0[iz+1] )*e1[i] +
                                                                ( ( dx*x1[ix] + x1[ix+1] )*z0[iz+1] + x0[ix+2]*z1[iz] )*e2[i] + x0[ix+2]*z0[iz+1]*e3[i] ) * y0[iy] ;
                    g013[i+13*MAXCBF] = ( y2[iy]*x1[ix]*e0[i] + ( ( dy*y1[iy] + y1[iy+1] )*x1[ix]   + y2[iy]*x0[ix+1] )*e1[i] +
                                                                ( ( dy*y1[iy] + y1[iy+1] )*x0[ix+1] + y0[iy+2]*x1[ix] )*e2[i] + y0[iy+2]*x0[ix+1]*e3[i] ) * z0[iz] ;
                    g013[i+14*MAXCBF] = ( x1[ix]*y1[iy]*z1[iz]*e0[i] + ( x1[ix]*y1[iy  ]*z0[iz+1] + x1[ix  ]*y0[iy+1]*z1[iz  ] + x0[ix+1]*y1[iy  ]*z1[iz] )*e1[i] +
                                                                       ( x1[ix]*y0[iy+1]*z0[iz+1] + x0[ix+1]*y1[iy  ]*z0[iz+1] + x0[ix+1]*y0[iy+1]*z1[iz] )*e2[i] + x0[ix+1]*y0[ix+1]*z0[iz+1]*e3[i] ) ;
                    g013[i+15*MAXCBF] = ( z2[iz]*x1[ix]*e0[i] + ( ( dz*z1[iz] + z1[iz+1] )*x1[ix]   + z2[iz]*x0[ix+1] )*e1[i] +
                                                                ( ( dz*z1[iz] + z1[iz+1] )*x0[ix+1] + z0[iz+2]*x1[ix] )*e2[i] + z0[iz+2]*x0[ix+1]*e3[i] ) * y0[iy] ;
                    g013[i+16*MAXCBF] = ( y3[iy] * e0[i] + ( y1[iy] + 2.0e+00*dy*y2[iy] + y2[iy+1] )*e1[i] + ( dy*dy*y1[iy] + dy*y1[iy+1] + y1[iy+2] )*e2[i] + y0[iy+3] * e3[i] ) * x0[ix] * z0[iz] ;
                    g013[i+17*MAXCBF] = ( y2[iy]*z1[iz]*e0[i] + ( ( dy*y1[iy] + y1[iy+1] )*z1[iz]   + y2[iy]*z0[iz+1] )*e1[i] +
                                                                ( ( dy*y1[iy] + y1[iy+1] )*z0[iz+1] + y0[iy+2]*z1[iz] )*e2[i] + y0[iy+2]*z0[iz+1]*e3[i] ) * x0[ix] ;
                    g013[i+18*MAXCBF] = ( z2[iz]*y1[iy]*e0[i] + ( ( dz*z1[iz] + z1[iz+1] )*y1[iy]   + z2[iz]*y0[iy+1] )*e1[i] +
                                                                ( ( dz*z1[iz] + z1[iz+1] )*y0[iy+1] + z0[iz+2]*y1[iy] )*e2[i] + z0[iz+2]*y0[iy+1]*e3[i] ) * x0[ix] ;
                    g013[i+19*MAXCBF] = ( z3[iz] * e0[i] + ( z1[iz] + 2.0e+00*dz*z2[iz] + z2[iz+1] )*e1[i] + ( dz*dz*z1[iz] + dz*z1[iz+1] + z1[iz+2] )*e2[i] + z0[iz+3] * e3[i] ) * x0[ix] * y0[iy] ;
                }
                /* . Put the values in the proper place. */
                iStart = iBasis->shells[iShell].nstartw ;
                for ( i = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
                {
                    Array2D_Item ( f   , i+iStart, g ) = g013[i] ;
                    Array2D_Item ( fX  , i+iStart, g ) = g013[i+   MAXCBF] ;
                    Array2D_Item ( fY  , i+iStart, g ) = g013[i+ 2*MAXCBF] ;
                    Array2D_Item ( fZ  , i+iStart, g ) = g013[i+ 3*MAXCBF] ;
                    Array2D_Item ( fXX , i+iStart, g ) = g013[i+ 4*MAXCBF] ;
                    Array2D_Item ( fXY , i+iStart, g ) = g013[i+ 5*MAXCBF] ;
                    Array2D_Item ( fXZ , i+iStart, g ) = g013[i+ 6*MAXCBF] ;
                    Array2D_Item ( fYY , i+iStart, g ) = g013[i+ 7*MAXCBF] ;
                    Array2D_Item ( fYZ , i+iStart, g ) = g013[i+ 8*MAXCBF] ;
                    Array2D_Item ( fZZ , i+iStart, g ) = g013[i+ 9*MAXCBF] ;
                    Array2D_Item ( fXXX, i+iStart, g ) = g013[i+10*MAXCBF] ;
                    Array2D_Item ( fXXY, i+iStart, g ) = g013[i+11*MAXCBF] ;
                    Array2D_Item ( fXXZ, i+iStart, g ) = g013[i+12*MAXCBF] ;
                    Array2D_Item ( fXYY, i+iStart, g ) = g013[i+13*MAXCBF] ;
                    Array2D_Item ( fXYZ, i+iStart, g ) = g013[i+14*MAXCBF] ;
                    Array2D_Item ( fXZZ, i+iStart, g ) = g013[i+15*MAXCBF] ;
                    Array2D_Item ( fYYY, i+iStart, g ) = g013[i+16*MAXCBF] ;
                    Array2D_Item ( fYYZ, i+iStart, g ) = g013[i+17*MAXCBF] ;
                    Array2D_Item ( fYZZ, i+iStart, g ) = g013[i+18*MAXCBF] ;
                    Array2D_Item ( fZZZ, i+iStart, g ) = g013[i+19*MAXCBF] ;
                }
            }
        }
    }
}
