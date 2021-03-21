/*==================================================================================================================================
! . Procedures for calculating QC dispersion interactions using the DFT-D2 model of S. Grimme (JCC 27, 1787-1799, 2006).
!=================================================================================================================================*/

# include <math.h>

# include "Boolean.h"
# include "Integer.h"
# include "QCDispersionDFTD2.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Energy and gradients.
! . All in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _LogTolerance 27.63102111592855 /* . Equivalent to - ln ( 10^(-12) ). */
Real QCDispersionDFTD2_Energy ( const Real          s6           ,
                                const Real          sR           ,
                                const Real          dR           ,
                                const RealArray1D  *sqrtC6       ,
                                const RealArray1D  *r0           ,
                                const Coordinates3 *coordinates3 ,
                                      Coordinates3 *gradients3   )
{
    auto Real energy = 0.0e+00 ;
    if ( ( sqrtC6 != NULL ) && ( r0 != NULL ) && ( coordinates3 != NULL ) )
    {
        auto Boolean  doGradients = ( gradients3 != NULL ) ;
        auto Integer  i, j ;
        auto Real     cI, cIJ, damp, dampF, dF, dX, dY, dZ, eLocal, expArg, gX, gY, gZ, r, rI, rIJ, r2, r6 ;
        for ( i = 1 ; i < View1D_Extent ( sqrtC6 ) ; i++ )
        {
            cI = Array1D_Item ( sqrtC6, i ) * s6 ;
            rI = Array1D_Item ( r0    , i )      ;
            for ( j = 0 ; j < i ; j++ )
            {
	        Coordinates3_DifferenceRow ( coordinates3, i, j, dX, dY, dZ ) ; /* rI - rJ */
                cIJ    = cI * Array1D_Item ( sqrtC6, j ) ;
                rIJ    = sR * ( rI + Array1D_Item ( r0, j ) ) ;
                r2     = dX * dX + dY * dY + dZ * dZ ;
                r6     = r2 * r2 * r2 ;
                r      = sqrt ( r2 ) ;
                expArg = dR * ( r / rIJ - 1.0e+00 ) ;
                     if ( expArg >  _LogTolerance ) { damp = 1.0e+00 ; dampF = 0.0e+00 ; }
                else if ( expArg < -_LogTolerance ) { damp = 0.0e+00 ; dampF = 0.0e+00 ; }
                else { dampF = exp ( -expArg ) ; damp = 1.0e+00 / ( 1.0e+00 + dampF ) ; }
                eLocal  = cIJ * damp / r6 ;
                energy -= eLocal ;
                if ( doGradients )
                {
                    dF = ( 6.0e+00 / r2 - ( dR * damp * dampF ) / ( r * rIJ ) ) * eLocal ;
                    gX = dF * dX ;
                    gY = dF * dY ;
                    gZ = dF * dZ ;
                    Coordinates3_IncrementRow ( gradients3, i, gX, gY, gZ ) ;
                    Coordinates3_DecrementRow ( gradients3, j, gX, gY, gZ ) ;
                }
            }
        }
    }
    return energy ;
}
# undef _LogTolerance
