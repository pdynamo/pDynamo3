# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_B1E1N0
# define _GAUSSIANBASISCONTAINERINTEGRALS_B1E1N0

# include "GaussianBasisContainer.h"
# include "IntegerArray1D.h"
# include "RealArray1D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_SelfOverlap ( const GaussianBasisContainer *self         ,
                                                          const IntegerArray1D         *basisIndices ,
                                                                RealArray1D            *selfOverlap  ,
                                                                Status                 *status       ) ;
# endif
