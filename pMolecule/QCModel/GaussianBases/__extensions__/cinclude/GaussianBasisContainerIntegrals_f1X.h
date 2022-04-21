# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_F1X
# define _GAUSSIANBASISCONTAINERINTEGRALS_F1X

# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "IntegerArray1D.h"
# include "RealArray1D.h"
# include "Status.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_f1Di ( const GaussianBasisContainer *self         ,
                                                   const Coordinates3           *coordinates3 ,
                                                         Vector3                *center       ,
                                                         RealArray1D            *dipoleX      ,
                                                         RealArray1D            *dipoleY      ,
                                                         RealArray1D            *dipoleZ      ,
                                                         Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1Oi ( const GaussianBasisContainer *self         ,
                                                         RealArray1D            *overlap      ,
                                                         Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1Qi ( const GaussianBasisContainer *self         ,
                                                   const Coordinates3           *coordinates3 ,
                                                         Vector3                *center       ,
                                                         RealArray1D            *qXX          ,
                                                         RealArray1D            *qYY          ,
                                                         RealArray1D            *qZZ          ,
                                                         RealArray1D            *qXY          ,
                                                         RealArray1D            *qXZ          ,
                                                         RealArray1D            *qYZ          ,
                                                         Status                 *status       ) ;
# endif
