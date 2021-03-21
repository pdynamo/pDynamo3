# ifndef _GAUSSIANBASISINTEGRALS_B2E1N1
# define _GAUSSIANBASISINTEGRALS_B2E1N1

# include "Coordinates3.h"
# include "GaussianBasis.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Selection.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_ElectronNuclear           ( const GaussianBasis *self       ,
                                                               const Real          *rS         ,
                                                               const GaussianBasis *other      ,
                                                               const Real          *rO         ,
                                                               const RealArray1D   *charges    ,
                                                               const RealArray1D   *widthsE    ,
                                                               const RealArray1D   *widthsN    ,
                                                               const Coordinates3  *rN         ,
                                                               const Selection     *selectionN ,
                                                                     RealArray2D   *integrals  ) ;
extern void GaussianBasisIntegrals_ElectronNuclearD          ( const GaussianBasis *self       ,
                                                               const Real          *rS         ,
                                                               const GaussianBasis *other      ,
                                                               const Real          *rO         ,
                                                               const RealArray1D   *charges    ,
                                                               const RealArray1D   *widthsE    ,
                                                               const RealArray1D   *widthsN    ,
                                                               const Coordinates3  *rN         ,
                                                               const Selection     *selectionN ,
                                                               const RealArray2D   *dOneIJ     ,
                                                                     Real          *dRi        ,
                                                                     Real          *dRj        ,
                                                                     Coordinates3  *gN         ) ;
extern void GaussianBasisIntegrals_ElectronNuclearPotentials ( const GaussianBasis *self       ,
                                                               const Real          *rS         ,
                                                               const GaussianBasis *other      ,
                                                               const Real          *rO         ,
                                                               const RealArray1D   *widthsE    ,
                                                               const RealArray1D   *widthsN    ,
                                                               const Coordinates3  *rN         ,
                                                               const Selection     *selectionN ,
                                                               const RealArray2D   *dOneIJ     ,
                                                                     RealArray1D   *potentials ) ;
# endif
