/*==================================================================================================================================
! . Integrals - 1 basis, 1 electron.
!=================================================================================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "GaussianBasisIntegrals_b1e1n0.h"
# include "Integer.h"
# include "NumericalMacros.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the self overlap integrals for a basis.
! . All zero or even 1-D polynomials are non-zero in Cartesians.
! . 1-D integrals for x^n (n even) are (n-1)!!/(2 a)^(n/2) (Pi/a)^(1/2).
! . For spherical harmonics all zero except for s-functions.
! . SelfOverlap should be appropriately initialized before entry to this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_SelfOverlap ( const GaussianBasis *iBasis, RealArray1D *selfOverlap )
{
    if ( ( iBasis != NULL ) && ( selfOverlap != NULL ) )
    {
        auto Integer  i, icbfind, ip, ishell, ix, iy, iz, ncfunci, t ;
        auto Real     ci, ei, si ;
        /* . Loop over shells. */
        for ( ishell = 0 ; ishell < iBasis->nshells ; ishell++ )
        {
            /* . Get information about the shell. */
            icbfind = iBasis->shells[ishell].type->cbfindex ;
            ncfunci = iBasis->shells[ishell].type->ncbf     ;
            /* . Loop over functions in the shell. */
            for ( i = 0 ; i < ncfunci ; i++ )
            {
   	        ix = CBFPOWX[i+icbfind] ;
	        iy = CBFPOWY[i+icbfind] ;
	        iz = CBFPOWZ[i+icbfind] ;
                /* . Zero or even polynomials only. */
                if ( IsEven ( ix ) && IsEven ( iy ) && IsEven ( iz ) )
                {
                    si = 0.0e+00 ;
                    for ( ip = 0 ; ip < iBasis->shells[ishell].nprimitives ; ip++ )
                    {
                        ci = iBasis->shells[ishell].primitives[ip].ccbf[i]  ;
                        ei = iBasis->shells[ishell].primitives[ip].exponent ;
                        for ( t = 1 ; t <= ( ix / 2 ) ; t++ ) ci *= ( Real ) ( 2 * t - 1 ) / ( 2.0e+00 * ei ) ;
                        for ( t = 1 ; t <= ( iy / 2 ) ; t++ ) ci *= ( Real ) ( 2 * t - 1 ) / ( 2.0e+00 * ei ) ;
                        for ( t = 1 ; t <= ( iz / 2 ) ; t++ ) ci *= ( Real ) ( 2 * t - 1 ) / ( 2.0e+00 * ei ) ;
                        si += ci / ( ei * sqrt ( ei ) ) ;
                    }
                    Array1D_Item ( selfOverlap, iBasis->shells[ishell].nstartw+i ) = PI32 * si ;
                }
            }
        }
    }
}
