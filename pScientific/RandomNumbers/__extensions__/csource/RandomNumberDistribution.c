/*==================================================================================================================================
! . This module implements various random number distributions.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "RandomNumberDistribution.h"

/*==================================================================================================================================
! . Gaussian/Normal distribution.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Box-Mueller method.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real RandomNumberDistribution_GaussianBoxMueller ( RandomNumberGenerator *rng, const Real mu, const Real sigma )
{
    if ( rng->hasGaussian )
    {
        rng->hasGaussian = False ;
        return ( sigma * rng->gaussian + mu ) ;
    }
    else
    {
        Real f, r2, x, y ;
        do
        {
            x  = -1.0 + 2.0 * RandomNumberGenerator_NextRealOpen ( rng ) ;
            y  = -1.0 + 2.0 * RandomNumberGenerator_NextRealOpen ( rng ) ;
            r2 = x * x + y * y ;
        }
        while ( r2 > 1.0 || r2 == 0.0 ) ;
        f = sqrt ( -2.0 * log ( r2 ) / r2 ) ;
        rng->gaussian    = f * x ;
        rng->hasGaussian = True  ;
        return ( sigma * f * y + mu ) ;
    }
}
