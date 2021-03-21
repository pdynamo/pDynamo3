# ifndef _GAUSSIANBASISNORMALIZE
# define _GAUSSIANBASISNORMALIZE

# include "Boolean.h"
# include "GaussianBasis.h"
# include "Real.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real GaussianBasis_Normalize (       GaussianBasis *self               ,
                                      const Boolean        checkNormalization ,
                                            Status        *status             ) ;
# endif
