/*==================================================================================================================================
! . Gaussian basis set rotation.
!=================================================================================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "GaussianBasisRotate.h"
# include "IntegerArray1D.h"
# include "IntegerArrayND.h"
# include "RealArray1D.h"

/*
! . The basis rotation matrices are not orthogonal.
! . To get the inverse rotation matrix, transpose the original 3x3 rotation and reconstruct the transformation. I.e.:

    r = 3x3 rotation matrix ;
    R = MakeRotation ( l, r   ) ;
    S = MakeRotation ( l, r^T ) ;
    R * S => I.

! . For orbitals       , have phi_i' = c_i^T * Tc * bfs, so c_i' = Tc^T * c_i.
! . For o-representaion, have c_i = X * a_i and a_i = Y^T * c_i, so a_i' = ( X^T * T * Y )^T a_i.

*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Simple real factorial - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real Factorial ( const Integer n )
{
    Real f = 1.0e+00 ;
    if ( n > 1 )
    {
        auto Integer i ;
        for ( i = 2 ; i <= n ; i++ ) { f *= ( Real ) i ; }
    }
    return f ;
 }
 
/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a rotation matrix, Tc, from l = 0, lMaximum given a 3x3 input rotation matrix, R.
! . Tc is appropriate for Cartesian basis sets.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_MakeLRotations ( const Integer      L      ,
                                    const Matrix33    *R      ,
                                          RealArray2D *Tc     ,
                                          Status      *status )
{
    if ( ( L  >= 0    ) &&
         ( R  != NULL ) &&
         ( Tc != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer d ;
        d = ( ( L + 1 ) * ( L + 2 ) * ( L + 3 ) ) / 6 ; /* . Minimum dimension of Tc. */
        if ( ( View2D_Columns ( Tc ) == d ) && ( View2D_Rows ( Tc ) == d ) )
        {
            auto Integer         shape[3] ;
            auto IntegerArray1D *indicesL ;
            auto IntegerArray2D *keys     ;
            auto IntegerArrayND *indicesK ;
            auto RealArray1D    *factors  ;
            auto RealArray2D    *pXYZ     ;
            shape[0] = shape[1] = shape[2] = L+1 ;
            factors  = RealArray1D_AllocateWithExtent     ( d  ,        status ) ; /* . Multinomial factors. */
            indicesK = IntegerArrayND_AllocateWithShape   ( 3  , shape, status ) ; /* . The n-values of the xyz combinations. */
            indicesL = IntegerArray1D_AllocateWithExtent  ( L+2,        status ) ; /* . Number of xyz combinations per l. */
            keys     = IntegerArray2D_AllocateWithExtents ( d  , 3    , status ) ; /* . x, y, z for each n-value. */
            pXYZ     = RealArray2D_AllocateWithExtents    ( d  , 3    , status ) ; /* . The x, y and z polynomial coefficients. */
            if ( Status_IsOK ( status ) )
            {
                auto Integer c, iX, iY, iZ, jX, jY, jZ, kX, kY, kZ, l, n, nX, nY, nZ, r, x, y, z ;
                auto Real    lF, mF, vX, vY, vZ, xF, yF, zF ;
                /* . Fill work arrays. */
                n = 0 ;
                for ( l = 0 ; l <= L ; l++ )
                {
                    lF = Factorial ( l ) ;
                    Array1D_Item ( indicesL, l ) = n ;
                    for ( z = 0 ; z <= l ; z++ )
                    {
                        zF = Factorial ( z ) ;
                        for ( y = 0 ; y <= ( l - z ) ; y++ )
                        {
                            x   = l - y - z ;
                            xF  = Factorial ( x ) ;
                            yF  = Factorial ( y ) ;
                            Array1D_Item   ( factors , n       ) = ( lF / ( xF * yF * zF ) ) ;
                            ArrayND_Item3D ( indicesK, x, y, z ) = n ;
                            Array2D_Item   ( keys    , n, 0    ) = x ;
                            Array2D_Item   ( keys    , n, 1    ) = y ;
                            Array2D_Item   ( keys    , n, 2    ) = z ;
                            n += 1 ;
                        }
                    }
                }
                Array1D_Item ( indicesL, L+1 ) = n ;
                /* . Create polynomials in x, y and z up to the required order, L. */
                /* . x'^l = ( Rxx x + Rxy y + Rxz z )^l, etc. */
                for ( r = 0 ; r < d ; r++ )
                {
                    mF = Array1D_Item ( factors, r ) ;
                    x  = Array2D_Item ( keys, r, 0 ) ;
                    y  = Array2D_Item ( keys, r, 1 ) ;
                    z  = Array2D_Item ( keys, r, 2 ) ;
                    Array2D_Item ( pXYZ, r, 0 ) = mF * pow ( Array2D_Item ( R, 0, 0 ), x ) * pow ( Array2D_Item ( R, 0, 1 ), y ) * pow ( Array2D_Item ( R, 0, 2 ), z ) ;
                    Array2D_Item ( pXYZ, r, 1 ) = mF * pow ( Array2D_Item ( R, 1, 0 ), x ) * pow ( Array2D_Item ( R, 1, 1 ), y ) * pow ( Array2D_Item ( R, 1, 2 ), z ) ;
                    Array2D_Item ( pXYZ, r, 2 ) = mF * pow ( Array2D_Item ( R, 2, 0 ), x ) * pow ( Array2D_Item ( R, 2, 1 ), y ) * pow ( Array2D_Item ( R, 2, 2 ), z ) ;
                }
                /* . Create transformations for each basis function. */
                RealArray2D_Set ( Tc, 0.0e+00 ) ;
                for ( r = 0 ; r < d ; r++ )
                {
                    x = Array2D_Item ( keys, r, 0 ) ;
                    y = Array2D_Item ( keys, r, 1 ) ;
                    z = Array2D_Item ( keys, r, 2 ) ;
                    for ( nX = Array1D_Item ( indicesL, x ) ; nX < Array1D_Item ( indicesL, x+1 ) ; nX++ )
                    {
                        iX = Array2D_Item ( keys, nX, 0 ) ;
                        jX = Array2D_Item ( keys, nX, 1 ) ;
                        kX = Array2D_Item ( keys, nX, 2 ) ;
                        vX = Array2D_Item ( pXYZ, nX, 0 ) ;
                        for ( nY = Array1D_Item ( indicesL, y ) ; nY < Array1D_Item ( indicesL, y+1 ) ; nY++ )
                        {
                            iY = Array2D_Item ( keys, nY, 0 ) ;
                            jY = Array2D_Item ( keys, nY, 1 ) ;
                            kY = Array2D_Item ( keys, nY, 2 ) ;
                            vY = Array2D_Item ( pXYZ, nY, 1 ) ;
                            for ( nZ = Array1D_Item ( indicesL, z ) ; nZ < Array1D_Item ( indicesL, z+1 ) ; nZ++ )
                            {
                                iZ = Array2D_Item ( keys, nZ, 0 ) ;
                                jZ = Array2D_Item ( keys, nZ, 1 ) ;
                                kZ = Array2D_Item ( keys, nZ, 2 ) ;
                                vZ = Array2D_Item ( pXYZ, nZ, 2 ) ;
                                c  = ArrayND_Item3D ( indicesK, iX+iY+iZ, jX+jY+jZ, kX+kY+kZ ) ;
                                Array2D_Item ( Tc, r, c ) += ( vX * vY * vZ ) ;
                            }
                        }
                    }
                }
            }
            /* . Finish up. */
            IntegerArray1D_Deallocate ( &indicesL ) ;
            IntegerArray2D_Deallocate ( &keys     ) ;
            IntegerArrayND_Deallocate ( &indicesK ) ;
            RealArray1D_Deallocate    ( &factors  ) ;
            RealArray2D_Deallocate    ( &pXYZ     ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate a rotation matrix, T, for a basis given an input rotation matrix, Tc, generated by GaussianBasis_MakeLRotations.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_MakeRotationMatrix (       GaussianBasis *self   ,
                                        const RealArray2D   *Tc     ,
                                        const Boolean        doC2O  ,
                                              RealArray2D   *T      ,
                                              Status        *status )
{
    if ( ( self != NULL ) &&
         ( Tc   != NULL ) &&
         ( T    != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean isOK ;
        auto Integer L, n ;
        L = self->maximum_angularmomentum ;
        n = ( ( L + 1 ) * ( L + 2 ) * ( L + 3 ) ) / 6 ; /* . Minimum dimension of Tc. */
        isOK = ( ( View2D_Columns ( Tc ) >= n ) && ( View2D_Rows ( Tc ) >= n ) ) ;
        if ( doC2O ) n = View2D_Columns ( self->c2o ) ;
        else         n = self->nbasisw ;
        isOK = isOK && ( ( View2D_Columns ( T ) == n ) && ( View2D_Rows ( T ) == n ) ) ;
        if ( isOK )
        {
            auto RealArray2D *Tl = T, *Tt = NULL ;
            if ( doC2O )
            {
                Tl = RealArray2D_AllocateWithExtents ( self->nbasisw, self->nbasisw, status ) ;
                Tt = RealArray2D_AllocateWithExtents ( n            , self->nbasisw, status ) ;
            }
            if ( Tl != NULL )
            {
                auto Integer     iShell, lMin, lStart, start ;
                auto RealArray2D cView, lView ;
                RealArray2D_Set ( Tl, 0.0e+00 ) ;
                for ( iShell = 0 ; iShell < self->nshells ; iShell++ )
                {
                    lMin   = self->shells[iShell].type->angularmomentum_low  ;
                    lStart = ( lMin * ( lMin + 1 ) * ( lMin + 2 ) ) / 6 ;
                    start  = self->shells[iShell].nstartw ;
                    n      = self->shells[iShell].nbasisw ;
                    RealArray2D_View   ( Tc, lStart, lStart, n, n, 1, 1, False, &cView, status ) ;
                    RealArray2D_View   ( Tl,  start,  start, n, n, 1, 1, False, &lView, status ) ;
                    RealArray2D_CopyTo ( &cView, &lView, status ) ;
                }
            }
            if ( doC2O )
            {
                RealArray2D_MatrixMultiply ( True , False, 1.0e+00, self->c2o, Tl, 0.0e+00, Tt, status ) ;
                RealArray2D_MatrixMultiply ( False, False, 1.0e+00, Tt, self->o2c, 0.0e+00, T , status ) ;
                RealArray2D_Deallocate ( &Tl ) ;
                RealArray2D_Deallocate ( &Tt ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}
