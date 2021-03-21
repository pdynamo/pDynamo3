# ifndef _RANDOMNUMBERDISTRIBUTIONS
# define _RANDOMNUMBERDISTRIBUTIONS

# include "RandomNumberGenerator.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Gaussian/Normal distribution. */
extern Real RandomNumberDistribution_GaussianBoxMueller ( RandomNumberGenerator *rng, const Real mu, const Real sigma ) ;

# endif
