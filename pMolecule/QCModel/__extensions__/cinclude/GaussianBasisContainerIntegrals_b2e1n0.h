# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_B2E1N0
# define _GAUSSIANBASISCONTAINERINTEGRALS_B2E1N0

# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "IntegerArray1D.h"
# include "RealArray1D.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_2Coulomb         ( const GaussianBasisContainer *self         ,
                                                               const IntegerArray1D         *basisIndices ,
                                                               const Coordinates3           *coordinates3 ,
                                                                     SymmetricMatrix        *integrals    ,
                                                                     Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_2CoulombD        ( const GaussianBasisContainer *self         ,
                                                               const IntegerArray1D         *basisIndices ,
                                                               const Coordinates3           *coordinates3 ,
                                                               const RealArray1D            *fPotential   ,
                                                               const RealArray1D            *wVector      ,
                                                                     Coordinates3           *gradients3   ,
                                                                     Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_2Overlap         ( const GaussianBasisContainer *self         ,
                                                               const IntegerArray1D         *basisIndices ,
                                                               const Coordinates3           *coordinates3 ,
                                                                     SymmetricMatrix        *integrals    ,
                                                                     Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_Dipole           ( const GaussianBasisContainer *self         ,
                                                               const IntegerArray1D         *basisIndices ,
                                                               const Coordinates3           *coordinates3 ,
                                                                     Vector3                *center       ,
                                                                     SymmetricMatrix        *dipoleX      ,
                                                                     SymmetricMatrix        *dipoleY      ,
                                                                     SymmetricMatrix        *dipoleZ      ,
                                                                     Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_Kinetic2Overlap  ( const GaussianBasisContainer *self         ,
                                                               const IntegerArray1D         *basisIndices ,
                                                               const Coordinates3           *coordinates3 ,
                                                                     SymmetricMatrix        *kinetic      ,
                                                                     SymmetricMatrix        *overlap      ,
                                                                     Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_Kinetic2OverlapD ( const GaussianBasisContainer *self         ,
                                                               const IntegerArray1D         *basisIndices ,
                                                               const Coordinates3           *coordinates3 ,
                                                               const SymmetricMatrix        *kDensity     ,
                                                               const SymmetricMatrix        *oDensity     ,
                                                                     Coordinates3           *gradients3   ,
                                                                     Status                 *status       ) ;
# endif
