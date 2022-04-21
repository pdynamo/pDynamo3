# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_F1XG1
# define _GAUSSIANBASISCONTAINERINTEGRALS_F1XG1

# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "IntegerArray1D.h"
# include "RealArray1D.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_f1Af1i   ( const GaussianBasisContainer *self         ,
                                                       const Coordinates3           *coordinates3 ,
                                                             SymmetricMatrix        *integrals    ,
                                                             Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1Cf1i   ( const GaussianBasisContainer *self         ,
                                                       const Coordinates3           *coordinates3 ,
                                                             SymmetricMatrix        *integrals    ,
                                                             Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1Df1i   ( const GaussianBasisContainer *self         ,
                                                       const Coordinates3           *coordinates3 ,
                                                             Vector3                *center       ,
                                                             SymmetricMatrix        *dipoleX      ,
                                                             SymmetricMatrix        *dipoleY      ,
                                                             SymmetricMatrix        *dipoleZ      ,
                                                             Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1KOf1i  ( const GaussianBasisContainer *self         ,
                                                       const Coordinates3           *coordinates3 ,
                                                             SymmetricMatrix        *kinetic      ,
                                                             SymmetricMatrix        *overlap      ,
                                                             Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1KOf1R1 ( const GaussianBasisContainer *self         ,
                                                       const Coordinates3           *coordinates3 ,
                                                       const SymmetricMatrix        *kDensity     ,
                                                       const SymmetricMatrix        *oDensity     ,
                                                             Coordinates3           *gradients3   ,
                                                             Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1Of1i   ( const GaussianBasisContainer *self         ,
                                                       const Coordinates3           *coordinates3 ,
                                                             SymmetricMatrix        *integrals    ,
                                                             Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1Qf1i   ( const GaussianBasisContainer *self         ,
                                                       const Coordinates3           *coordinates3 ,
                                                             Vector3                *center       ,
                                                             SymmetricMatrix        *qXX          ,
                                                             SymmetricMatrix        *qYY          ,
                                                             SymmetricMatrix        *qZZ          ,
                                                             SymmetricMatrix        *qXY          ,
                                                             SymmetricMatrix        *qXZ          ,
                                                             SymmetricMatrix        *qYZ          ,
                                                             Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1Xf1R1  ( const GaussianBasisContainer *self         ,
                                                       const Coordinates3           *coordinates3 ,
                                                       const RealArray1D            *aVector      ,
                                                       const RealArray1D            *xVector      ,
                                                       const GaussianBasisOperator   operator     ,
                                                             Coordinates3           *gradients3   ,
                                                             Status                 *status       ) ;
# endif
