# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_F1CP1
# define _GAUSSIANBASISCONTAINERINTEGRALS_F1CP1

# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "IntegerArray1D.h"
# include "RealArray1D.h"
# include "Selection.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_f1Cm1R1 ( const GaussianBasisContainer *self            ,
                                                      const RealArray1D            *charges         ,
                                                      const RealArray1D            *widthsE         ,
                                                      const RealArray1D            *widthsN         ,
                                                      const Coordinates3           *coordinates3    ,
                                                      const Coordinates3           *coordinates3G   ,
                                                            Selection              *selectionG      ,
                                                      const RealArray1D            *fitCoefficients ,
                                                            Coordinates3           *gradients3      ,
                                                            Coordinates3           *gradients3G     ,
                                                            Status                 *status          ) ;
extern void GaussianBasisContainerIntegrals_f1Cm1V  ( const GaussianBasisContainer *self            ,
                                                      const RealArray1D            *charges         ,
                                                      const RealArray1D            *widthsE         ,
                                                      const RealArray1D            *widthsN         ,
                                                      const Coordinates3           *coordinates3    ,
                                                      const Coordinates3           *coordinates3G   ,
                                                            Selection              *selectionG      ,
                                                            RealArray1D            *integrals       ,
                                                            Status                 *status          ) ;
extern void GaussianBasisContainerIntegrals_f1Cp1V  ( const GaussianBasisContainer *self            ,
                                                      const RealArray1D            *widthsE         ,
                                                      const RealArray1D            *widthsN         ,
                                                      const Coordinates3           *coordinates3    ,
                                                      const Coordinates3           *coordinates3G   ,
                                                            Selection              *selectionG      ,
                                                      const RealArray1D            *fitCoefficients ,
                                                            RealArray1D            *potentials      ,
                                                            Status                 *status          ) ;
# endif
