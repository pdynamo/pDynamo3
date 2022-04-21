# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_F1XG2
# define _GAUSSIANBASISCONTAINERINTEGRALS_F1XG2

# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "IntegerArray1D.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_f1Xg2i  ( const GaussianBasisContainer *self         ,
                                                      const GaussianBasisContainer *other        ,
                                                      const Coordinates3           *coordinates3 ,
                                                      const GaussianBasisOperator   operator     ,
                                                            BlockStorage           *fitIntegrals ,
                                                            Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1Xg2R1 ( const GaussianBasisContainer *self         ,
                                                      const GaussianBasisContainer *other        ,
                                                      const Coordinates3           *coordinates3 ,
                                                      const SymmetricMatrix        *density      ,
                                                      const RealArray1D            *xVector      ,
                                                      const GaussianBasisOperator   operator     ,
                                                            Coordinates3           *gradients3   ,
                                                            Status                 *status       ) ;
# endif
