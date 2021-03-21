# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_B3E1N0
# define _GAUSSIANBASISCONTAINERINTEGRALS_B3E1N0

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
extern void GaussianBasisContainerIntegrals_ElectronFit  ( const GaussianBasisContainer *self         ,
                                                           const IntegerArray1D         *selfIndices  ,
                                                           const GaussianBasisContainer *other        ,
                                                           const IntegerArray1D         *otherIndices ,
                                                           const Coordinates3           *coordinates3 ,
                                                                 BlockStorage           *fitIntegrals ,
                                                                 Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_ElectronFitD ( const GaussianBasisContainer *self         ,
                                                           const IntegerArray1D         *selfIndices  ,
                                                           const GaussianBasisContainer *other        ,
                                                           const IntegerArray1D         *otherIndices ,
                                                           const Coordinates3           *coordinates3 ,
                                                           const SymmetricMatrix        *sDensity     ,
                                                           const RealArray1D            *oDensity     ,
                                                                 Coordinates3           *gradients3   ,
                                                                 Status                 *status       ) ;
# endif
